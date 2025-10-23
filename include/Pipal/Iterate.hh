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
#include "Pipal/Parameter.hh"
#include "Pipal/Counter.hh"

namespace Pipal
{

  // Constructor
  inline void buildIterate(Iterate & z, Parameter & p, Input & i, Counter & c)
  {
    // Initialize quantities
    // z.rho_ = p.rho_init;
    // z.f = 0.0;
    // z.fu = 0.0;
    z.g.setZero(i.nV);
    // z.r1.setZero(i.nE);
    // z.r2.setZero(i.nE);
    // z.cE.setZero(i.nE);
    z.JE.resize(i.nE, i.nV); // SparseMatrix
    // z.JEnnz = 0;
    // z.lE.setZero(i.nE);
    // z.s1.setZero(i.nI);
    // z.s2.setZero(i.nI);
    // z.cI.setZero(i.nI);
    z.JI.resize(i.nI, i.nV); // SparseMatrix
    // z.JInnz = 0;
    // z.lI.setConstant(i.nI, 0.5);
    z.H.resize(i.nV, i.nV); // SparseMatrix
    // z.Hnnz = 0;
    // z.v = 0.0;
    // z.vu = 0.0;
    // z.v0 = 1.0;
    // z.phi = 0.0;
    // z.Annz = 0;
    // z.shift = 0.0;
    z.b.resize(i.nA);
    z.kkt.setZero(3);
    z.kkt_.setConstant(p.opt_err_mem, INFTY);
    // z.err = 0;
    // z.fs = 1.0;
    // z.cEs.setOnes(i.nE);
    // z.cEu.setZero(i.nE);
    // z.cIs.setOnes(i.nI);
    // z.cIu.setZero(i.nI);
    z.A.resize(i.nA, i.nA);
    // z.shift22 = 0.0;
    // z.v_ = 0.0;
    // z.cut_ = false;

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
    z.Hnnz  = static_cast<Integer>(z.H.nonZeros());
    z.JEnnz = static_cast<Integer>(z.JE.nonZeros());
    z.JInnz = static_cast<Integer>(z.JI.nonZeros());
    initNewtonMatrix(z, i);
    evalNewtonMatrix(z, p, i, c);
  }

  // Termination checker
  inline Integer checkTermination(Iterate const & z, Parameter const & p, Input const & i, Counter const & c)
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
  inline void evalDependent(Iterate & z, Parameter & p, Input & i)
  {
    // Evaluate quantities dependent on penalty and interior-point parameters
    evalSlacks(z, p, i);
    evalMerit(z, i);
    evalKKTErrors(z, i);
  }

  // Function evaluator
  inline void evalFunctions(Iterate & z, Input & i, Counter & c)
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
      z.cI.head(i.n3) = i.l3 - z.x.segment(i.n1, i.n3);
    }
    if (i.n4 > 0) {
      z.cI.segment(i.n3, i.n4) = -i.u4 + z.x.segment(i.n1+i.n3, i.n4);
    }
    if (i.n5 > 0) {
      z.cI.segment(i.n3+i.n4,      i.n5) =  i.l5 - z.x.segment(i.n1+i.n3+i.n4, i.n5);
      z.cI.segment(i.n3+i.n4+i.n5, i.n5) = -i.u5 + z.x.segment(i.n1+i.n3+i.n4, i.n5);
    }
    if (i.n7 > 0) {
      z.cI.segment(i.n3+i.n4+i.n5+i.n5, i.n7) = i.l7 - c_orig(i.I7);
    }
    if (i.n8 > 0) {
      z.cI.segment(i.n3+i.n4+i.n5+i.n5+i.n7, i.n8) = -i.u8 + c_orig(i.I8);
    }
    if (i.n9 > 0) {
      z.cI.segment(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8,      i.n9) =  i.l9 - c_orig(i.I9);
      z.cI.segment(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n9) = -i.u9 + c_orig(i.I9);
    }

    // Store unscaled quantities
    z.fu = z.f;
    if (i.nE > 0) {z.cEu = z.cE;}
    if (i.nI > 0) {z.cIu = z.cI;}

    // Scale quantities
    z.f *= z.fs;
    if (i.nE > 0) {z.cE = (z.cEs*z.cE).matrix();}
    if (i.nI > 0) {z.cI = (z.cIs*z.cI).matrix();}

  }


  // Create the function to inject a block into a sparse matrix with given indices
  inline void insert_block(SparseMatrix & mat_sparse, Matrix const & mat_dense, Integer const row_offset,
    Integer const col_offset)
  {
    const Integer cols{static_cast<Integer>(mat_dense.cols())}, rows{static_cast<Integer>(mat_dense.rows())};
    for (Integer r{0}; r < rows; ++r) {
      for (Integer c{0}; c < cols; ++c) {
        mat_sparse.coeffRef(r + row_offset, c + col_offset) = mat_dense(r, c);
      }
    }
  }

  // Create the function to inject a block into a sparse matrix with given indices
  inline void insert_block(SparseMatrix & mat_sparse, Matrix const & mat_dense, Indices const & idx_row,
    Indices const & idx_col, Integer const row_offset, Integer const col_offset)
  {
    for (Integer r_idx{0}; r_idx < idx_row.size(); ++r_idx) {
      Integer r{idx_row[r_idx] + row_offset};
      for (Integer c_idx{0}; c_idx < idx_col.size(); ++c_idx) {
          
        mat_sparse.coeffRef(r, idx_col[c_idx] + col_offset) = mat_dense.coeff(r_idx, c_idx);
      }
    }
  }

  // Gradient evaluator
  inline void evalGradients(Iterate & z, Input & i, Counter & c)
  {
    // Evaluate x in original space
    Vector x_orig;
    evalXOriginal(z, i, x_orig);

    // Initialize/Reset evaluation flag
    z.err = 0;

    // Increment gradient evaluation counter
    incrementGradientCount(c);

    // Try AMPL gradients evaluation
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
      Integer col_offset{0};
      insert_block(z.JE, J_orig(i.I6, i.I1), 0, col_offset);
      col_offset += i.I1.size();
      insert_block(z.JE, J_orig(i.I6, i.I3), 0, col_offset);
      col_offset += i.I3.size();
      insert_block(z.JE, J_orig(i.I6, i.I4), 0, col_offset);
      col_offset += i.I4.size();
      insert_block(z.JE, J_orig(i.I6, i.I5), 0, col_offset);
    }

    // Initialize inequality constraint Jacobian
    //if (i.nI > 0) {z.JI.resize(i.nI, i.nV);}

    // Set inequality constraint Jacobian
    if (i.n3 > 0) {
      Integer diag_len{std::min(i.n3, i.n1+i.n3)};
      for (Integer k{0}; k < diag_len; ++k)  {z.JI.coeffRef(k, i.n1+k) = -1.0;}
    }
    if (i.n4 > 0) {
      Integer tmp{i.n1+i.n3};
      for (Integer k{0}; k < i.n4; ++k) {z.JI.coeffRef(i.n3+k, tmp+k) = 1.0;}
    }
    if (i.n5 > 0) {
      Integer tmp_row{i.n3+i.n4}, tmp_col{i.n1+i.n3+i.n4};
      for (Integer k{0}; k < i.n5; ++k) {
        z.JI.coeffRef(tmp_row+k,      tmp_col+k) = -1.0;
        z.JI.coeffRef(tmp_row+i.n5+k, tmp_col+k) = 1.0;
      }
    }
    if (i.n7 > 0) {
      Integer row_offset{i.n3+i.n4+i.n5+i.n5}, col_offset{0};
      insert_block(z.JI, -J_orig(i.I7, i.I1), i.I7, i.I1, row_offset, col_offset);
      col_offset += i.I1.size();
      insert_block(z.JI, J_orig(i.I7, i.I3), i.I7, i.I3, row_offset, col_offset);
      col_offset += i.I3.size();
      insert_block(z.JI, J_orig(i.I7, i.I4), i.I7, i.I4, row_offset, col_offset);
      col_offset += i.I4.size();
      insert_block(z.JI, J_orig(i.I7, i.I5), i.I7, i.I5, row_offset, col_offset);
    }
    if (i.n8 > 0) {
      Integer row_offset{i.n3+i.n4+i.n5+i.n5+i.n7}, col_offset{0};
      insert_block(z.JI, J_orig(i.I8, i.I1), i.I8, i.I1, row_offset, col_offset);
      col_offset += i.I1.size();
      insert_block(z.JI, J_orig(i.I8, i.I3), i.I8, i.I3, row_offset, col_offset);
      col_offset += i.I3.size();
      insert_block(z.JI, J_orig(i.I8, i.I4), i.I8, i.I4, row_offset, col_offset);
      col_offset += i.I4.size();
      insert_block(z.JI, J_orig(i.I8, i.I5), i.I8, i.I5, row_offset, col_offset);
    }
    if (i.n9 > 0) {
      Integer row_offset{i.n3+i.n4+i.n5+i.n5+i.n7+i.n8}, col_offset{0};
      insert_block(z.JI, -J_orig(i.I9, i.I1), i.I9, i.I1, row_offset, col_offset);
      col_offset += i.I1.size();
      insert_block(z.JI, -J_orig(i.I9, i.I3), i.I9, i.I3, row_offset, col_offset);
      col_offset += i.I3.size();
      insert_block(z.JI, -J_orig(i.I9, i.I4), i.I9, i.I4, row_offset, col_offset);
      col_offset += i.I4.size();
      insert_block(z.JI, -J_orig(i.I9, i.I5), i.I9, i.I5, row_offset, col_offset);
      row_offset += i.n9; // next row block
      col_offset = 0;
      insert_block(z.JI, J_orig(i.I9, i.I1), i.I9, i.I1, row_offset, col_offset);
      col_offset += i.I1.size();
      insert_block(z.JI, J_orig(i.I9, i.I3), i.I9, i.I3, row_offset, col_offset);
      col_offset += i.I3.size();
      insert_block(z.JI, J_orig(i.I9, i.I4), i.I9, i.I4, row_offset, col_offset);
      col_offset += i.I4.size();
      insert_block(z.JI, J_orig(i.I9, i.I5), i.I9, i.I5, row_offset, col_offset);
    }

    // Scale objective gradient
    z.g *= z.fs;

    // Scale constraint Jacobians
    if (i.nE > 0) {
      for (Integer k{0}; k < z.JE.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(z.JE, k); it; ++it) {it.valueRef() *= z.cEs[it.row()];}
      }
    }
    if (i.nI > 0) {
      for (Integer k{0}; k < z.JI.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(z.JI, k); it; ++it) {it.valueRef() *= z.cIs[it.row()];}
      }
    }
  }

  // Hessian evaluator
  inline void evalHessian(Iterate & z, Input & i, Counter & c)
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
      return;
    }

    // Set Hessian of the Lagrangian
    Integer row_offset{0}, col_offset{0};
    insert_block(z.H, H_orig(i.I1, i.I1), i.I1, i.I1, row_offset, col_offset);
    col_offset += i.I1.size();
    insert_block(z.H, H_orig(i.I1, i.I3), i.I1, i.I3, row_offset, col_offset);
    col_offset += i.I3.size();
    insert_block(z.H, H_orig(i.I1, i.I4), i.I1, i.I4, row_offset, col_offset);
    col_offset += i.I4.size();
    insert_block(z.H, H_orig(i.I1, i.I5), i.I1, i.I5, row_offset, col_offset);
    row_offset += i.I1.size();
    col_offset = 0; // next row block
    insert_block(z.H, H_orig(i.I3, i.I1), i.I3, i.I1, row_offset, col_offset);
    col_offset += i.I1.size();
    insert_block(z.H, H_orig(i.I3, i.I3), i.I3, i.I3, row_offset, col_offset);
    col_offset += i.I3.size();
    insert_block(z.H, H_orig(i.I3, i.I4), i.I3, i.I4, row_offset, col_offset);
    col_offset += i.I4.size();
    insert_block(z.H, H_orig(i.I3, i.I5), i.I3, i.I5, row_offset, col_offset);
    row_offset += i.I3.size();
    col_offset = 0; // next row block
    insert_block(z.H, H_orig(i.I4, i.I1), i.I4, i.I1, row_offset, col_offset);
    col_offset += i.I1.size();
    insert_block(z.H, H_orig(i.I4, i.I3), i.I4, i.I3, row_offset, col_offset);
    col_offset += i.I3.size();
    insert_block(z.H, H_orig(i.I4, i.I4), i.I4, i.I4, row_offset, col_offset);
    col_offset += i.I4.size();
    insert_block(z.H, H_orig(i.I4, i.I5), i.I4, i.I5, row_offset, col_offset);
    row_offset += i.I4.size();
    col_offset = 0; // next row block
    insert_block(z.H, H_orig(i.I5, i.I1), i.I5, i.I1, row_offset, col_offset);
    col_offset += i.I1.size();
    insert_block(z.H, H_orig(i.I5, i.I3), i.I5, i.I3, row_offset, col_offset);
    col_offset += i.I3.size();
    insert_block(z.H, H_orig(i.I5, i.I4), i.I5, i.I4, row_offset, col_offset);
    col_offset += i.I4.size();
    insert_block(z.H, H_orig(i.I5, i.I5), i.I5, i.I5, row_offset, col_offset);

    // Rescale H
    z.H *= z.rho*z.fs;
  }

  // Infeasibility evaluator
  inline void evalInfeasibility(Iterate & z, Input const & i)
  {
    // Evaluate scaled and unscaled feasibility violations
    z.v  = evalViolation(i, z.cE, z.cI) / std::max(1.0, z.v0);
    z.vu = evalViolation(i, z.cEu, z.cIu);
  }

  // KKT error evaluator
  inline Real evalKKTError(Iterate & z, Input const & i, Real const rho, Real const mu)
  {
    // Initialize optimality vector
    Vector kkt(i.nV+2*i.nE+2*i.nI);
    kkt.setZero();

    // Set gradient of penalty objective
    kkt.head(i.nV) = rho*z.g;

    // Set gradient of Lagrangian for constraints
    if (i.nE > 0) {kkt.head(i.nV) += (z.lE.matrix().transpose()*z.JE).transpose();}
    if (i.nI > 0) {kkt.head(i.nV) += (z.lI.matrix().transpose()*z.JI).transpose();}

    // Set complementarity for constraint slacks
    if (i.nE > 0) {
      kkt.segment(i.nV, 2*i.nE) << z.r1*(1.0 + z.lE) - mu, z.r2*(1.0 - z.lE) - mu;
    }
    if (i.nI > 0) {
      kkt.segment(i.nV+2*i.nE, 2*i.nI) << z.s1*z.lI - mu, z.s2*(1.0 - z.lI) - mu;
    }

    // Scale complementarity
    if (rho > 0.0) {
      kkt *= (1.0 / std::max(1.0, (rho*z.g).template lpNorm<Eigen::Infinity>()));
    }

    // Evaluate optimality error
    return kkt.template lpNorm<Eigen::Infinity>();
  }

  // KKT errors evaluator
  inline void evalKKTErrors(Iterate & z, Input const & i)
  {
    // Loop to compute optimality errors
    z.kkt << evalKKTError(z, i, 0, 0), evalKKTError(z, i, z.rho, 0), evalKKTError(z, i, z.rho, z.mu);
  }

  // Evaluator of lambda in original space
  inline void evalLambdaOriginal(Iterate const & z, Input const & i, Vector & l)
  {
    // Initialize multipliers in original space
    l.setZero(i.nE+i.n7+i.n8+i.n9);

    // Scale equality constraint multipliers
    if (i.nE > 0) {l(i.I6) = z.lE*(z.cEs/(z.rho*z.fs));}

    // Scale inequality constraint multipliers
    Array lI;
    if (i.n7+i.n8+i.n9 > 0) {lI = z.lI*(z.cIs/(z.rho*z.fs));}

    // Set inequality constraint multipliers in original space
    if (i.n7 > 0) {
      l(i.I7) = -lI.segment(i.n3+i.n4+i.n5+i.n5, i.n7);
    }
    if (i.n8 > 0) {
      l(i.I8) = lI.segment(i.n3+i.n4+i.n5+i.n5+i.n7, i.n8);
    }
    if (i.n9 > 0) {
      l(i.I9) = lI.segment(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n9)
                -lI.segment(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8, i.n9);
    }
  }

  // Matrices evaluator
  inline void evalMatrices(Iterate & z, Parameter & p, Input & i, Counter & c)
  {
    // Evaluate Hessian and Newton matrices
    evalHessian(z, i, c);
    evalNewtonMatrix(z, p, i, c);
  }

  // Merit evaluator
  inline void evalMerit(Iterate & z, Input const & i)
  {
    // Initialize merit for objective
    z.phi = z.rho*z.f;

    // Update merit for slacks
    if (i.nE > 0) {
      Array r_all(2*i.nE); r_all << z.r1, z.r2;
      z.phi -= z.mu * r_all.log().sum() - r_all.sum();
    }
    if (i.nI > 0) {
      Array s_all(2*i.nI); s_all << z.s1, z.s2;
      z.phi -= z.mu * s_all.log().sum() - z.s2.sum();
    }
  }

  // Newton matrix evaluator
  inline void evalNewtonMatrix(Iterate & z, Parameter & p, Input const & i, Counter & c)
  {
    // Check for equality constraints
    if (i.nE > 0)
    {
      // Set diagonal terms
      for (Integer j{0}; j < i.nE; ++j) {
        z.A.coeffRef(i.nV+j, i.nV+j)           = (1.0 + z.lE(j))/z.r1(j);
        z.A.coeffRef(i.nV+i.nE+j, i.nV+i.nE+j) = (1.0 - z.lE(j))/z.r2(j);
      }

      // Set constraint Jacobian
      insert_block(z.A, z.JE, i.nV+2*i.nE+2*i.nI, 0);
    }

    // Check for inequality constraints
    if (i.nI > 0)
    {
      // Set diagonal terms
      for (Integer j{0}; j < i.nI; ++j) {
        z.A.coeffRef(i.nV+2*i.nE+j,      i.nV+2*i.nE+j)      = z.lI(j)/z.s1(j);
        z.A.coeffRef(i.nV+2*i.nE+i.nI+j, i.nV+2*i.nE+i.nI+j) = (1.0 - z.lI(j))/z.s2(j);
      }

      // Set constraint Jacobian
      insert_block(z.A, z.JI, i.nV+3*i.nE+2*i.nI, 0);
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
      insert_block(z.A, z.H+z.shift*Matrix::Identity(i.nV, i.nV), 0, 0);;

      // Set diagonal terms
      for (Integer j{0}; j < i.nE; ++j) {
        z.A.coeffRef(i.nV+2*i.nE+2*i.nI+j, i.nV+2*i.nE+2*i.nI+j) = -z.shift22;
      }

      // Set diagonal terms
      for (Integer j{0}; j < i.nI; ++j) {
        z.A.coeffRef(i.nV+3*i.nE+2*i.nI+j, i.nV+3*i.nE+2*i.nI+j) = -z.shift22;
      }

      // Set number of nonzeros in (upper triangle of) Newton matrix
      z.Annz = static_cast<Integer>(z.A.nonZeros());

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
    z.H.diagonal().array() += z.shift;
  }

  // Newton right-hand side evaluator
  inline void evalNewtonRhs(Iterate & z, Input const & i)
  {
    // Initialize right-hand side vector
    z.b.setZero();

    // Set gradient of objective
    for (Integer k{0}; k < i.nV; ++k) {
      for (Integer k{0}; k < i.nV; ++k) {z.b.coeffRef(k) = z.rho*z.g(k);}

      // Set gradient of Lagrangian for constraints
      if (i.nE > 0) {
        Vector tmp((z.lE.matrix().transpose()*z.JE).transpose());
        for (Integer k{0}; k < tmp.size(); ++k) {z.b.coeffRef(k) += tmp[k];}
      }
      if (i.nI > 0) {
        Vector tmp((z.lI.matrix().transpose()*z.JI).transpose());
        for (Integer k{0}; k < tmp.size(); ++k) {z.b.coeffRef(k) += tmp[k];}
      }
    }

    // Set complementarity for constraint slacks
    if (i.nE > 0) {
      // Compute element-wise complementarity terms with safe element-wise division
      Vector tmp(2*i.nE);
      tmp << 1.0 + z.lE - z.mu * z.r1.cwiseInverse(), 1.0 - z.lE - z.mu * z.r2.cwiseInverse();
      for (Integer k{0}; k < 2*i.nE; ++k) {z.b.insert(i.nV + k) = tmp[k];}
    }
    if (i.nI > 0) {
      // Compute element-wise complementarity terms with safe element-wise division
      Vector tmp(2*i.nI);
      tmp << z.lI - z.mu * z.s1.cwiseInverse(), 1.0 - z.lI - z.mu * z.s2.cwiseInverse();
      for (Integer k{0}; k < 2*i.nI; ++k) {z.b.insert(i.nV + 2*i.nE + k) = tmp[k];}
    }

    // Set penalty-interior-point constraint values
    if (i.nE > 0) {
      const Vector tmp(z.cE + z.r1 - z.r2);
      const Integer offset{i.nV+2*i.nE+2*i.nI};
      for (Integer k{0}; k < i.nE; ++k) {z.b.coeffRef(offset+k) = tmp[k];}
    }
    if (i.nI > 0) {
      const Vector tmp(z.cI + z.s1 - z.s2);
      const Integer offset{i.nV+3*i.nE+2*i.nI};
      for (Integer k{0}; k < i.nI; ++k) {z.b.coeffRef(offset+k) = tmp[k];}
    }
  }

  // Scalings evaluator
  inline void evalScalings(Iterate & z, Parameter & p, Input & i, Counter & c)
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
      Real row_inf_norm{0.0};
      for (SparseMatrix::InnerIterator it(z.JE, j); it; ++it) {
        row_inf_norm = std::max(row_inf_norm, std::abs(it.value()));
      }
      z.cEs(j) = p.grad_max / std::max(row_inf_norm, p.grad_max);
    }

    // Loop through inequality constraints
    for (Integer j{0}; j < i.nI; ++j)
    {
      // Scale down inequality constraint j if norm of gradient is too large
      Real row_inf_norm{0.0};
      for (SparseMatrix::InnerIterator it(z.JI, j); it; ++it) {
        row_inf_norm = std::max(row_inf_norm, std::abs(it.value()));
      }
      z.cIs(j) = p.grad_max / std::max(row_inf_norm, p.grad_max);
    }
  }

  // Slacks evaluator
  inline void evalSlacks(Iterate & z, Parameter & p, Input const & i)
  {
    // Check for equality constraints
    if (i.nE > 0)
    {
      // Set slacks
      z.r1 = 0.5*((z.mu - z.cE) + (z.cE.square() + z.mu*z.mu).sqrt());
      z.r2 = 0.5*((z.mu + z.cE) + (z.cE.square() + z.mu*z.mu).sqrt());

      // Adjust for numerical error
      z.r1.cwiseMax(p.slack_min);
      z.r2.cwiseMax(p.slack_min);
    }

    // Check for inequality constraints
    if (i.nI > 0)
    {
      // Set slacks
      z.s1 = 0.5*((2.0*z.mu - z.cI) + (z.cI.square() + 4.0*z.mu*z.mu).sqrt());
      z.s2 = 0.5*((2.0*z.mu + z.cI) + (z.cI.square() + 4.0*z.mu*z.mu).sqrt());

      // Adjust for numerical error
      z.s1.cwiseMax(p.slack_min);
      z.s2.cwiseMax(p.slack_min);
    }
  }

  // Evaluator of x in original space
  inline void evalXOriginal(Iterate & z, Input const & i, Vector & x)
  {
    // Initialize x in original space
    x.setZero(i.n0);

    // Evaluate x in original space
    x(i.I1) = z.x.head(i.n1);
    x(i.I2) = i.b2;
    x(i.I3) = z.x.segment(i.n1, i.n3);
    x(i.I4) = z.x.segment(i.n1+i.n3, i.n4);
    x(i.I5) = z.x.segment(i.n1+i.n3+i.n4, i.n5);
  }

  // Gets primal-dual point
  inline void getSolution(Iterate & z, Input & i, Vector & x, Vector & l)
  {
    evalXOriginal(z, i, x);
    evalLambdaOriginal(z, i, l);
  }

  // Initializes Newton matrix
  inline void initNewtonMatrix(Iterate &z, Input const &i)
  {
    // Allocate memory
    z.A.reserve(z.Hnnz + 5*i.nE + 5*i.nI + z.JEnnz + z.JInnz);

    // Initialize interior-point Hessians
    {
      Integer diag;
      for (Integer k{0}; k < 2*i.nE; ++k) {
        diag = i.nV+k;
        z.A.coeffRef(diag, diag) = 1.0;
      }

      for (Integer k{0}; k < 2*i.nI; ++k) {
        diag = i.nV+2*i.nE+k;
        z.A.coeffRef(diag, diag) = 1.0;
      }
    }

    // Check for constraints
    if (i.nE > 0)
    {
      // Initialize constraint Jacobian
      Integer row, col;
      for (Integer k{0}; k < i.nE; ++k) {
        row = i.nV+2*i.nE+2*i.nI+k, col = i.nV+k;
        z.A.coeffRef(row, col) = 1.0;
        z.A.coeffRef(row, col+i.nE) = -1.0;
      }
    }

    // Check for inequality constraints
    if (i.nI > 0)
    {
      // Initialize constraint Jacobian
      Integer row, col;
      for (Integer k{0}; k < i.nI; ++k) {
        row = i.nV+3*i.nE+2*i.nI+k, col = i.nV+2*i.nE+k;
        z.A.coeffRef(row, col) = 1.0;
        z.A.coeffRef(row, col+i.nI) = -1.0;
      }
    }

    // Optional: compress after all insertions for efficient arithmetic
    //z.A.makeCompressed();
  }

  // Set interior-point parameter
  inline void setMu(Iterate & z, Real const mu) {z.mu = mu;}

  // Set primal variables
  inline void setPrimals(Iterate & z, Input const & i, Vector const & x, Array const & r1,
    Array const & r2, Array const  & s1, Array const & s2, Array const  & lE, Array const & lI,
    Real const f, Array const  & cE, Array const & cI, Real const phi)
  {
    // Set primal variables
    z.x = x; z.f = f;
    if (i.nE > 0) {z.cE = cE; z.r1 = r1; z.r2 = r2; z.lE = lE;}
    if (i.nI > 0) {z.cI = cI; z.s1 = s1; z.s2 = s2; z.lI = lI;}
    z.phi = phi;
  }

  // Set penalty parameter
  inline void setRho(Iterate & z, Real const rho) {z.rho = rho;}

  // Set last penalty parameter
  inline void setRhoLast(Iterate & z, Real const rho) {z.rho_ = rho;}

  // Iterate updater
  inline void updateIterate(Iterate & z, Parameter & p, Input & i, Counter & c, Direction const & d, Acceptance const  & a)
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
    //z.kkt_.resize(p.opt_err_mem);
    z.kkt_ << z.kkt(1), z.kkt_.head(p.opt_err_mem-1);
  }

  // Parameter updater
  inline void updateParameters(Iterate & z, Parameter & p, Input & i)
  {
    // Check for interior-point parameter update based on optimality error
    while (z.mu > p.mu_min && z.kkt(2) <= std::max(z.mu, p.opt_err_tol-z.mu))
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
  inline void updatePoint(Iterate & z, Input const & i, Direction const & d, Acceptance const & a)
  {
    // Update primal and dual variables
    z.x += a.p*d.x ;
    if (i.nE > 0) {z.r1 += a.p*d.r1; z.r2 += a.p*d.r2;}
    if (i.nI > 0) {z.s1 += a.p*d.s1; z.s2 += a.p*d.s2;}
    if (i.nE > 0) {z.lE += a.d*d.lE;}
    if (i.nI > 0) {z.lI += a.d*d.lI;}
  }

  // Feasibility violation evaluator
  inline Real evalViolation(Input const & i, Array const & cE, Array const & cI)
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
