/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Pipal project is distributed under the MIT License.                                       *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_PIPAL_ITERATE_HH
#define INCLUDE_PIPAL_ITERATE_HH

// Pipal includes
#include "Pipal/Types.hh"
#include "Pipal/Input.hh"
#include "Pipal/Parameter.hh"
#include "Pipal/Counter.hh"

namespace Pipal
{

  // Constructor
  void buildIterate(Iterate & z, Parameter & p, Input & i, Counter & c)
  {
    // Initialize quantities
    z.rho_ = p.rho_init;
    z.f = 0.0;
    z.fu = 0.0;
    z.g.setZero(i.nV);
    z.r1.setZero(i.nE);
    z.r2.setZero(i.nE);
    z.cE.setZero(i.nE);
    z.JE.setZero(i.nE, i.nV); // SparseMatrix
    z.JEnnz = 0;
    z.lE.setZero(i.nE);
    z.s1.setZero(i.nI);
    z.s2.setZero(i.nI);
    z.cI.setZero(i.nI);
    z.JI.setZero(i.nI, i.nV); // SparseMatrix
    z.JInnz = 0;
    z.lI.setConstant(i.nI, 0.5);
    z.H.setZero(i.nV, i.nV); // SparseMatrix
    z.Hnnz = 0;
    z.v = 0.0;
    z.vu = 0.0;
    z.v0 = 1.0;
    z.phi = 0.0;
    z.Annz = 0;
    z.shift = 0.0;
    z.b.setZero(i.nA);
    z.kkt.setZero(3);
    z.kkt_.setConstant(INFTY, 1, p.opt_err_mem);
    z.err = 0;
    z.fs = 1.0;
    z.cEs.setOnes(i.nE);
    z.cEu.setZero(i.nE);
    z.cIs.setOnes(i.nI);
    z.cIu.setZero(i.nI);
    z.A.setZero(i.nA, i.nA);
    z.shift22 = 0.0;
    z.v_ = 0.0;
    z.cut_ = false;

    z.x     = i.x0;
    z.rho   = p.rho_init;
    z.mu    = p.mu_init;
    z.lE.setZero(i.nE);
    z.lI.setConstant(i.nI, 0.5);
    z.err   = 0;
    evalScalings(z, p, i, c);
    evalFunctions(z, i, c);
    evalGradients(z, i, c);
    evalDependent(z, p, i);
    z.v0    = 1.0;
    evalInfeasibility(z, i);
    z.v0    = z.v;
    evalInfeasibility(z, i);
    z.v_    = z.v;
    evalHessian(z, i, c);
    z.Hnnz  = nnz(z.H);
    z.JEnnz = nnz(z.JE);
    z.JInnz = nnz(z.JI);
    initNewtonMatrix(z, i);
    evalNewtonMatrix(z, p, i, c);
  }

  // Termination checker
  Integer checkTermination(Iterate & z, Parameter & p, Input & i, Counter & c)
  {
    // Update termination based on optimality error of nonlinear optimization problem
    if (z.kkt(1) <= p.opt_err_tol && z.v <= p.opt_err_tol) {return 1;}

    // Update termination based on optimality error of feasibility problem
    if (z.kkt(0) <= p.opt_err_tol && z.v > p.opt_err_tol) {return 2;}

    // Update termination based on iteration count
    if (c.k >= p.iter_max) {return 3;}

    // Update termination based on invalid bounds
    if (i.vi > 0) {return 4;}

    // Update termination based on function evaluation error
    if (z.err > 0) {return 5;}

    return 0;
  }

  // Dependent quantity evaluator
  void evalDependent(Iterate & z, Parameter & p, Input & i)
  {
    // Evaluate quantities dependent on penalty and interior-point parameters
    evalSlacks(z, p, i);
    evalMerit(z, i);
    evalKKTErrors(z, i);
  }

  // Function evaluator
  void evalFunctions(Iterate & z, Input & i, Counter & c)
  {
    // Evaluate x in original space
    Vector x_orig;
    evalXOriginal(z, i, x_orig);

    // Initialize/Reset evaluation flag
    z.err = 0;

    // Increment function evaluation counter
    incrementFunctionCount(c);

    // Try AMPL functions evaluation
    Vector c_orig;
    try
    {
      // Evaluate AMPL functions
      i.f_fun(x_orig, z.f);
      i.c_fun(x_orig, c_orig);
    }
    catch (...)
    {
      // Set evaluation flag, default values, and return
      z.err = 1;
      z.f   = QUIET_NAN;
      z.cE.setConstant(i.nE, QUIET_NAN);
      z.cI.setConstant(i.nI, QUIET_NAN);
      z.fu  = QUIET_NAN;
      z.cEu.setConstant(i.nE, QUIET_NAN);
      z.cIu.setConstant(i.nI, QUIET_NAN);
      return;
    }

    // Set equality constraint values
    if (i.nE > 0) {z.cE = c_orig(i.I6) - i.b6;}

    // Initialize inequality constraint values
    if (i.nI > 0) {z.cI.setZero(i.nI);}

    // Set inequality constraint values
    if (i.n3 > 0) {
      z.cI(Eigen::seq(0, i.n3-1)) = i.l3 - z.x(Eigen::seq(i.n1, i.n1+i.n3-1));
    }
    if (i.n4 > 0) {
      z.cI(Eigen::seq(i.n3, i.n3+i.n4-1)) = -i.u4 + z.x(Eigen::seq(i.n1+i.n3, i.n1+i.n3+i.n4-1));
    }
    if (i.n5 > 0) {
      z.cI(Eigen::seq(i.n3+i.n4, i.n3+i.n4+i.n5-1)) = i.l5 - z.x(Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1));
      z.cI(Eigen::seq(i.n3+i.n4+i.n5, i.n3+i.n4+i.n5+i.n5-1)) = -i.u5 + z.x(Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1));
    }
    if (i.n7 > 0) {
      z.cI(Eigen::seq(i.n3+i.n4+i.n5+i.n5, i.n3+i.n4+i.n5+i.n5+i.n7-1)) = i.l7 - c_orig(i.I7);
    }
    if (i.n8 > 0) {
      z.cI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8-1)) = -i.u8 + c_orig(i.I8);
    }
    if (i.n9 > 0) {
      z.cI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9-1)) =  i.l9 - c_orig(i.I9);
      z.cI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9+i.n9-1)) = -i.u9 + c_orig(i.I9);
    }

    // Store unscaled quantities
    z.fu = z.f;
    if (i.nE > 0) {z.cEu = z.cE;}
    if (i.nI > 0) {z.cIu = z.cI;}

    // Scale quantities
    z.f *= z.fs;
    if (i.nE > 0) {z.cE = (z.cEs.array()*z.cE.array()).matrix();}
    if (i.nI > 0) {z.cI = (z.cIs.array()*z.cI.array()).matrix();}

  }

  // Gradient evaluator
  void evalGradients(Iterate & z, Input & i, Counter & c)
  {
    // Evaluate x in original space
    Vector x_orig;
    evalXOriginal(z, i, x_orig);

    // Initialize/Reset evaluation flag
    z.err = 0;

    // Increment gradient evaluation counter
    incrementGradientCount(c);

    // Try AMPL gradients evaluation
    z.g.resize(i.nV);
    z.JE.resize(i.nE, i.nV);
    z.JI.resize(i.nI, i.nV);
    Vector g_orig;
    Matrix J_orig;
    try
    {
      // Evaluate AMPL gradients
      i.g_fun(x_orig, g_orig);
      i.J_fun(x_orig, J_orig);
    }
    catch (...)
    {
      // Set evaluation flag, default values, and return
      z.err = 1;
      return;
    }

    // Set objective gradient
    z.g << g_orig(i.I1), g_orig(i.I3), g_orig(i.I4), g_orig(i.I5);

    // Set equality constraint Jacobian
    if (i.nE > 0) {
      z.JE << J_orig(i.I6, i.I1), J_orig(i.I6, i.I3), J_orig(i.I6, i.I4), J_orig(i.I6, i.I5);
    }

    // Initialize inequality constraint Jacobian
    if (i.nI > 0) {z.JI.resize(i.nI, i.nV);}

    // Set inequality constraint Jacobian
    if (i.n3 > 0) {
      z.JI(Eigen::seq(0, i.n3-1), Eigen::seq(i.n1, i.n1+i.n3-1)) = -SparseMatrix::Identity(i.n3, i.n3);
    }
    if (i.n4 > 0) {
      z.JI(Eigen::seq(i.n3, i.n3+i.n4-1), Eigen::seq(i.n1+i.n3, i.n1+i.n3+i.n4-1)) = SparseMatrix::Identity(i.n4, i.n4);
    }
    if (i.n5 > 0) {
      z.JI(Eigen::seq(i.n3+i.n4, i.n3+i.n4+i.n5-1), Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1)) = -SparseMatrix::Identity(i.n5, i.n5);
      z.JI(Eigen::seq(i.n3+i.n4+i.n5, i.n3+i.n4+i.n5+i.n5-1), Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1)) = SparseMatrix::Identity(i.n5, i.n5);
    }
    if (i.n7 > 0) {
      z.JI(Eigen::seq(i.n3+i.n4+i.n5+i.n5, i.n3+i.n4+i.n5+i.n5+i.n7-1), Eigen::seq(0, i.n1+i.n3+i.n4+i.n5-1)) <<
        -J_orig(i.I7, i.I1), J_orig(i.I7,i.I3), J_orig(i.I7,i.I4), J_orig(i.I7,i.I5);
    }
    if (i.n8 > 0) {
      z.JI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8-1), Eigen::seq(0, i.n1+i.n3+i.n4+i.n5-1)) <<
        J_orig(i.I8,i.I1), J_orig(i.I8,i.I3), J_orig(i.I8,i.I4), J_orig(i.I8,i.I5);
    }
    if (i.n9 > 0) {
      z.JI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9-1), Eigen::seq(0, i.n1+i.n3+i.n4+i.n5-1)) <<
        -J_orig(i.I9,i.I1), J_orig(i.I9,i.I3), J_orig(i.I9,i.I4), J_orig(i.I9,i.I5);
      z.JI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9+i.n9-1), Eigen::seq(0, i.n1+i.n3+i.n4+i.n5-1)) <<
        J_orig(i.I9,i.I1), J_orig(i.I9,i.I3), J_orig(i.I9,i.I4), J_orig(i.I9,i.I5);
    }

    // Scale objective gradient
    z.g *= z.fs;

    // Scale constraint Jacobians
    if (i.nE > 0) {z.JE = z.cEs.asDiagonal() * z.JE;} // SPARSE
    if (i.nI > 0) {z.JI = z.cIs.asDiagonal() * z.JI;} // SPARSE
  }

  // Hessian evaluator
  void evalHessian(Iterate & z, Input & i, Counter & c)
  {
    // Evaluate lambda in original space
    Vector l_orig, x_orig;
    evalLambdaOriginal(z, i, l_orig);
    evalXOriginal(z, i, x_orig);

    // Initialize/Reset evaluation flag
    z.err = 0;

    // Increment Hessian evaluation counter
    incrementHessianCount(c);

    // Try AMPL Hessian evaluation
    Matrix H_orig;
    try
    {
      // Evaluate H_orig
      i.H_fun(x_orig, l_orig, H_orig);
    }
    catch (...)
    {
      // Set evaluation flag, default values, and return
      z.err = 1;
      z.H.resize(i.nV, i.nV);
      return;
    }

    // Set Hessian of the Lagrangian
    z.H <<
      H_orig(i.I1, i.I1), H_orig(i.I1, i.I3), H_orig(i.I1, i.I4), H_orig(i.I1, i.I5),
      H_orig(i.I3, i.I1), H_orig(i.I3, i.I3), H_orig(i.I3, i.I4), H_orig(i.I3, i.I5),
      H_orig(i.I4, i.I1), H_orig(i.I4, i.I3), H_orig(i.I4, i.I4), H_orig(i.I4, i.I5),
      H_orig(i.I5, i.I1), H_orig(i.I5, i.I3), H_orig(i.I5, i.I4), H_orig(i.I5, i.I5);

    // Rescale H
    z.H *= z.rho*z.fs;
  }

  // Infeasibility evaluator
  void evalInfeasibility(Iterate & z, Input & i)
  {
    // Evaluate scaled and unscaled feasibility violations
    z.v  = evalViolation(i, z.cE, z.cI) / std::max(1.0, z.v0);
    z.vu = evalViolation(i, z.cEu, z.cIu);
  }

  // KKT error evaluator
  Real evalKKTError(Iterate & z, Input & i, Real const rho, Real const mu)
  {
    // Initialize optimality vector
    Vector kkt(i.nV+2*i.nE+2*i.nI);
    kkt.setZero();

    // Set gradient of penalty objective
    kkt(Eigen::seq(0, i.nV-1)) = rho*z.g;

    // Set gradient of Lagrangian for constraints
    if (i.nE > 0) {kkt(Eigen::seq(0, i.nV-1)) = kkt(Eigen::seq(0, i.nV-1)) + (z.lE.transpose()*z.JE).transpose();} // OPTMIZE
    if (i.nI > 0) {kkt(Eigen::seq(0, i.nV-1)) = kkt(Eigen::seq(0, i.nV-1)) + (z.lI.transpose()*z.JI).transpose();} // OPTMIZE

    // Set complementarity for constraint slacks
    if (i.nE > 0) {
      kkt(Eigen::seq(i.nV, i.nV+2*i.nE-1)) <<
        (z.r1.array() * (1.0 + z.lE.array()).array()) - mu,
        (z.r2.array() * (1.0 - z.lE.array()).array()) - mu;
    }
    if (i.nI > 0) {
      kkt(Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1)) <<
        (z.s1.array() * (0.0 + z.lI.array()).array()) - mu,
        (z.s2.array() * (1.0 - z.lI.array()).array()) - mu;
    }

    // Scale complementarity
    if (rho > 0) {kkt = (1.0 / std::max(1.0, (rho*z.g).template lpNorm<Eigen::Infinity>()))*kkt;}

    // Evaluate optimality error
    return kkt.template lpNorm<Eigen::Infinity>();
  }

  // KKT errors evaluator
  void evalKKTErrors(Iterate & z, Input & i)
  {
    // Loop to compute optimality errors
    z.kkt(0) = evalKKTError(z, i, 0, 0);
    z.kkt(1) = evalKKTError(z, i, z.rho, 0);
    z.kkt(2) = evalKKTError(z, i, z.rho, z.mu);
  }

  // Evaluator of lambda in original space
  void evalLambdaOriginal(Iterate & z, Input & i, Vector & l)
  {
    // Initialize multipliers in original space
    l.setZero(i.nE+i.n7+i.n8+i.n9);

    // Scale equality constraint multipliers
    Vector lE;
    if (i.nE > 0) {lE = (z.lE.array()*(z.cEs/(z.rho*z.fs)).array()).matrix();}

    // Set equality constraint multipliers in original space
    if (i.nE > 0) {l(i.I6) = lE;}

    // Scale inequality constraint multipliers
    Vector lI;
    if (i.n7+i.n8+i.n9 > 0) {lI = (z.lI.array()*(z.cIs/(z.rho*z.fs)).array()).matrix();}

    // Set inequality constraint multipliers in original space
    if (i.n7 > 0) {
      l(i.I7) = -lI(Eigen::seq(i.n3+i.n4+i.n5+i.n5, i.n3+i.n4+i.n5+i.n5+i.n7-1));
    }
    if (i.n8 > 0) {
      l(i.I8) = lI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8-1));
    }
    if (i.n9 > 0) {
      l(i.I9) = lI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9+i.n9-1))
                -lI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9-1));
    }
  }

  // Matrices evaluator
  void evalMatrices(Iterate & z, Parameter & p, Input & i, Counter & c)
  {
    // Evaluate Hessian and Newton matrices
    evalHessian(z, i, c);
    evalNewtonMatrix(z, p, i, c);
  }

  // Merit evaluator
  void evalMerit(Iterate & z, Input & i)
  {
    // Initialize merit for objective
    z.phi = z.rho*z.f;

    // Update merit for slacks
    if (i.nE > 0) {
      Vector r_all(z.r1.size() + z.r2.size()); r_all << z.r1, z.r2;
      z.phi -= z.mu * r_all.array().log().sum() + r_all.sum();
    }
    if (i.nI > 0) {
      Eigen::VectorXd s_all(z.s1.size() + z.s2.size()); s_all << z.s1, z.s2;
      z.phi -= z.mu * s_all.array().log().sum() + z.s2.sum();
    }
  }

  // Newton matrix evaluator
  void evalNewtonMatrix(Iterate & z, Parameter & p, Input & i, Counter & c)
  {
    // Check for equality constraints
    if (i.nE > 0)
    {
      // Set diagonal terms
      for (Integer j{0}; j < i.nE; ++j) {
        z.A(i.nV+j, i.nV+j)           = (1.0 + z.lE(j))/z.r1(j);
        z.A(i.nV+i.nE+j, i.nV+i.nE+j) = (1.0 - z.lE(j))/z.r2(j);
      }

      // Set constraint Jacobian
      z.A(Eigen::seq(i.nV+2*i.nE+2*i.nI, i.nV+3*i.nE+2*i.nI-1), Eigen::seq(0, i.nV-1)) = z.JE;
    }

    // Check for inequality constraints
    if (i.nI > 0)
    {
      // Set diagonal terms
      for (Integer j{0}; j < i.nI; ++j) {
        z.A(i.nV+2*i.nE+j, i.nV+2*i.nE+j)           = (0.0 + z.lI(j))/z.s1(j);
        z.A(i.nV+2*i.nE+i.nI+j, i.nV+2*i.nE+i.nI+j) = (1.0 - z.lI(j))/z.s2(j);
      }

      // Set constraint Jacobian
      z.A(Eigen::seq(i.nV+3*i.nE+2*i.nI, i.nV+3*i.nE+3*i.nI-1), Eigen::seq(0, i.nV-1)) = z.JI;
    }

    // Set minimum potential shift
    Real min_shift{std::max(p.shift_min, p.shift_factor1*z.shift)};

    // Initialize Hessian modification
    if (z.cut_ == 1) {z.shift = std::min(p.shift_max,min_shift/p.shift_factor2);}
    else {z.shift = 0;}

    // Initialize inertia correction loop
    bool done{false};
    z.shift22 = 0;

    // Loop until inertia is correct
    while (!done && z.shift < p.shift_max)
    {
      // Set Hessian of Lagrangian
      z.A(Eigen::seq(0, i.nV-1), Eigen::seq(0, i.nV-1)) = z.H+z.shift*SparseMatrix::Identity(i.nV, i.nV);

      // Set diagonal terms
      for (Integer j{0}; j < i.nE; ++j) {
        z.A(i.nV+2*i.nE+2*i.nI+j, i.nV+2*i.nE+2*i.nI+j) = -z.shift22;
      }

      // Set diagonal terms
      for (Integer j{0}; j < i.nI; ++j) {
        z.A(i.nV+3*i.nE+2*i.nI+j, i.nV+3*i.nE+2*i.nI+j) = -z.shift22;
      }

      // Set number of nonzeros in (upper triangle of) Newton matrix
      z.Annz = nnz<MatrixView::TRIL>(z.A);

      // Factor primal-dual matrix
      z.ldlt.compute(z.A);

      // Approximate number of negative pivots (inertia)
      Integer neig{static_cast<Integer>((z.ldlt.vectorD().array() < 0.0).count())};

      // Increment factorization counter
      incrementFactorizationCount(c);

      // Set number of nonnegative eigenvalues
      Integer peig{i.nA - neig};

      // Check inertia
      if (peig < i.nV+2*i.nE+2*i.nI) {z.shift = std::max(min_shift,z.shift/p.shift_factor2);}
      else if (neig < i.nE+i.nI && z.shift22 == 0) {z.shift22 = p.shift_min;}
      else {done = true;}
    }

    // Update Hessian
    z.H = z.H+z.shift*SparseMatrix::Identity(i.nV, i.nV);
  }

  // Newton right-hand side evaluator
  void evalNewtonRhs(Iterate & z, Input & i)
  {
    // Initialize right-hand side vector
    z.b.setZero(i.nA);

    // Set gradient of objective
    z.b(Eigen::seq(0, i.nV-1)) = z.rho*z.g;

    // Set gradient of Lagrangian for constraints
    if (i.nE > 0) {z.b(Eigen::seq(0, i.nV-1)) = z.b(Eigen::seq(0, i.nV-1)) + (z.lE.transpose()*z.JE).transpose();}
    if (i.nI > 0) {z.b(Eigen::seq(0, i.nV-1)) = z.b(Eigen::seq(0, i.nV-1)) + (z.lI.transpose()*z.JI).transpose();}

    // Set complementarity for constraint slacks
    if (i.nE > 0) {
      // compute element-wise complementarity terms with safe element-wise division
      Vector tmpE1((1.0 + z.lE.array() - z.mu * z.r1.cwiseInverse().array()).matrix());
      Vector tmpE2((1.0 - z.lE.array() - z.mu * z.r2.cwiseInverse().array()).matrix());
      z.b(Eigen::seq(i.nV, i.nV+2*i.nE-1)) << tmpE1, tmpE2;
    }
    if (i.nI > 0) {
      // compute element-wise complementarity terms with safe element-wise division
      Vector tmpI1((z.lI.array() - z.mu * z.s1.cwiseInverse().array()).matrix());
      Vector tmpI2((1.0 - z.lI.array() - z.mu * z.s2.cwiseInverse().array()).matrix());
      z.b(Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1)) << tmpI1, tmpI2;
    }

    // Set penalty-interior-point constraint values
    if (i.nE > 0) {z.b(Eigen::seq(i.nV+2*i.nE+2*i.nI, i.nV+3*i.nE+2*i.nI-1)) = z.cE + z.r1 - z.r2;}
    if (i.nI > 0) {z.b(Eigen::seq(i.nV+3*i.nE+2*i.nI, i.nV+3*i.nE+3*i.nI-1)) = z.cI + z.s1 - z.s2;}
  }

  // Scalings evaluator
  void evalScalings(Iterate & z, Parameter & p, Input & i, Counter & c)
  {
    // Initialize scalings
    z.fs = 1;
    z.cEs.setOnes(i.nE);
    z.cIs.setOnes(i.nI);

    // Evaluate gradients
    evalGradients(z, i, c);

    // Scale down objective if norm of gradient is too large
    z.fs = p.grad_max / std::max(z.g.template lpNorm<Eigen::Infinity>(), p.grad_max);

    // Loop through equality constraints
    for (Integer j{0}; j < i.nE; ++j)
    {
      // Scale down equality constraint j if norm of gradient is too large
      z.cEs(j) = p.grad_max / std::max(z.JE.row(j).template lpNorm<Eigen::Infinity>(), p.grad_max);
    }

    // Loop through inequality constraints
    for (Integer j{0}; j < i.nI; ++j)
    {
      // Scale down inequality constraint j if norm of gradient is too large
      z.cIs(j) = p.grad_max / std::max(z.JI.row(j).template lpNorm<Eigen::Infinity>(), p.grad_max);
    }
  }

  // Slacks evaluator
  void evalSlacks(Iterate & z, Parameter & p, Input & i)
  {
    // Check for equality constraints
    if (i.nE > 0)
    {
      // Set slacks
      z.r1 = 0.5*((z.mu - z.cE.array()) + (z.cE.array().square() + z.mu*z.mu).sqrt()).matrix();
      z.r2 = 0.5*((z.mu + z.cE.array()) + (z.cE.array().square() + z.mu*z.mu).sqrt()).matrix();

      // Adjust for numerical error
      z.r1.cwiseMax(p.slack_min);
      z.r2.cwiseMax(p.slack_min);
    }

    // Check for inequality constraints
    if (i.nI > 0)
    {
      // Set slacks
      z.s1 = 0.5*((2.0*z.mu - z.cI.array()) + (z.cI.array().square() + 4.0*z.mu*z.mu).sqrt()).matrix();
      z.s2 = 0.5*((2.0*z.mu + z.cI.array()) + (z.cI.array().square() + 4.0*z.mu*z.mu).sqrt()).matrix();

      // Adjust for numerical error
      z.s1.cwiseMax(p.slack_min);
      z.s2.cwiseMax(p.slack_min);
    }
  }

  // Evaluator of x in original space
  void evalXOriginal(Iterate & z, Input & i, Vector & x)
  {
    // Initialize x in original space
    x.setZero(i.n0);

    // Evaluate x in original space
    x(i.I1) = z.x(Eigen::seq(0, i.n1-1));
    x(i.I2) = i.b2;
    x(i.I3) = z.x(Eigen::seq(i.n1, i.n1+i.n3-1));
    x(i.I4) = z.x(Eigen::seq(i.n1+i.n3, i.n1+i.n3+i.n4-1));
    x(i.I5) = z.x(Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1));
  }

  // Gets primal-dual point
  void getSolution(Iterate & z, Input & i, Vector & x, Vector & l)
  {
    evalXOriginal(z, i, x);
    evalLambdaOriginal(z, i, l);
  }

  // Initializes Newton matrix
  void initNewtonMatrix(Iterate & z, Input & i)
  {
    // Allocate memory
    z.A.resize(i.nA, i.nA);// FIXME, z.Hnnz+5*i.nE+5*i.nI+z.JEnnz+z.JInnz);

    // Initialize interior-point Hessians
    z.A(Eigen::seq(i.nV, i.nV+2*i.nE-1), Eigen::seq(i.nV, i.nV+2*i.nE-1)) = SparseMatrix::Identity(2*i.nE, 2*i.nE);
    z.A(Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1), Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1)) = SparseMatrix::Identity(2*i.nI, 2*i.nI);

    // Check for equality constraints
    if (i.nE > 0)
    {
      // Initialize constraint Jacobian
      z.A(Eigen::seq(i.nV+2*i.nE+2*i.nI, i.nV+3*i.nE+2*i.nI-1), Eigen::seq(i.nV, i.nV+i.nE-1)) =  SparseMatrix::Identity(i.nE, i.nE);
      z.A(Eigen::seq(i.nV+2*i.nE+2*i.nI, i.nV+3*i.nE+2*i.nI-1), Eigen::seq(i.nV+i.nE, i.nV+2*i.nE-1)) = -SparseMatrix::Identity(i.nE, i.nE);
    }

    // Check for inequality constraints
    if (i.nI > 0)
    {
      // Initialize constraint Jacobian
      z.A(Eigen::seq(i.nV+3*i.nE+2*i.nI, i.nV+3*i.nE+3*i.nI-1), Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+i.nI-1)) =  SparseMatrix::Identity(i.nI, i.nI);
      z.A(Eigen::seq(i.nV+3*i.nE+2*i.nI, i.nV+3*i.nE+3*i.nI-1), Eigen::seq(i.nV+2*i.nE+i.nI, i.nV+2*i.nE+2*i.nI-1)) = -SparseMatrix::Identity(i.nI, i.nI);
    }
  }

  // Set interior-point parameter
  void setMu(Iterate & z, Real const mu) {z.mu = mu;}

  // Set primal variables
  void setPrimals(Iterate & z, Input & i, Vector & x, Vector & r1, Vector & r2,
    Vector & s1, Vector & s2, Vector & lE, Vector & lI, Real const f,
    Vector & cE, Vector & cI, Real const phi)
  {
    // Set primal variables
    z.x = x; z.f = f;
    if (i.nE > 0) {z.cE = cE; z.r1 = r1; z.r2 = r2; z.lE = lE;}
    if (i.nI > 0) {z.cI = cI; z.s1 = s1; z.s2 = s2; z.lI = lI;}
    z.phi = phi;
  }

  // Set penalty parameter
  void setRho(Iterate & z, Real const rho) {z.rho = rho;}

  // Set last penalty parameter
  void setRhoLast(Iterate & z, Real const rho) {z.rho_ = rho;}

  // Iterate updater
  void updateIterate(Iterate & z, Parameter & p, Input & i, Counter & c,
    Direction & d, Acceptance & a)
  {
    // Update last quantities
    z.v_   = z.v;
    z.cut_ = (a.p < a.p0);

    // Update iterate quantities
    updatePoint(z, i, d, a);
    evalInfeasibility(z, i);
    evalGradients(z, i, c);
    evalDependent(z, p, i);

    // Update last KKT errors
    z.kkt_.resize(p.opt_err_mem);
    z.kkt_ << z.kkt(1), z.kkt_(Eigen::seq(0, p.opt_err_mem-2));
  }

  // Parameter updater
  void updateParameters(Iterate & z, Parameter & p, Input & i)
  {
    // Check for interior-point parameter update based on optimality error
    while (z.mu > p.mu_min && z.kkt(2) <= std::max({z.mu, p.opt_err_tol-z.mu}))
    {
      // Restrict interior-point parameter increase
      setMuMaxExpZero(p);

      // Update interior-point parameter
      if (z.mu > p.mu_min)
      {
        // Decrease interior-point
        z.mu = std::max(p.mu_min, std::min(p.mu_factor*z.mu, std::pow(z.mu, p.mu_factor_exp)));

        // Evaluate penalty and interior-point parameter dependent quantities
        evalDependent(z, p, i);
      }
    }

    // Check for penalty parameter update based on optimality error
    if ((z.kkt(1) <= p.opt_err_tol && z.v > p.opt_err_tol) ||
          z.v > std::max({1.0, z.v_, p.infeas_max}))
    {
      // Update penalty parameter
      if (z.rho > p.rho_min)
      {
        // Decrease penalty parameter
        z.rho = std::max(p.rho_min, p.rho_factor*z.rho);

        // Evaluate penalty and interior-point parameter dependent quantities
        evalDependent(z, p, i);
      }
    }
  }

  // Primal point updater
  void updatePoint(Iterate & z, Input & i, Direction & d, Acceptance & a)
  {
    // Update primal and dual variables
    z.x += a.p*d.x ;
    if (i.nE > 0) {z.r1 += a.p*d.r1; z.r2 += a.p*d.r2;}
    if (i.nI > 0) {z.s1 += a.p*d.s1; z.s2 += a.p*d.s2;}
    if (i.nE > 0) {z.lE += a.d*d.lE;}
    if (i.nI > 0) {z.lI += a.d*d.lI;}
  }

  // Feasibility violation evaluator
  Real evalViolation(Input & i, Vector & cE, Vector & cI)
  {
    // Initialize violation vector
    Vector vec;

    // Update vector for constraint values
    if (i.nE > 0) {vec = cE;}
    if (i.nI > 0) {
      Vector cIpos(cI.cwiseMax(0.0));
      if (vec.size() > 0) {
        Vector tmp(vec.size() + cIpos.size());
        tmp << vec, cIpos;
        vec = std::move(tmp);
      } else {
        vec = cIpos;
      }
    }

    // Evaluate vector norm
    return vec.template lpNorm<1>();
  }

} // namespace Pipal

#endif /* INCLUDE_PIPAL_ITERATE_HH */
