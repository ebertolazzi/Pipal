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
#include "Pipal/Defines.hh"
#include "Pipal/Problem.hh"
#include "Pipal/Parameter.hh"
#include "Pipal/Input.hh"
#include "Pipal/Counter.hh"

// Standard libraries
#include <type_traits>
#include <limits>

// Eigen library
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Pipal
{
  /**
   * \brief Class for managing the current iterate of the solver.
   * \tparam Real The floating-point type.
   * \tparam Integer The integer type.
   */
  template<typename Real, typename Integer>
  struct Iterate
  {
    static_assert(std::is_floating_point<Real>::value,
      "Pipal::Iterate: template argument 'Real' must be a floating-point type.");
    static_assert(std::is_integral<Integer>::value,
      "Pipal::Iterate: template argument 'Integer' must be an integer type.");

    static constexpr Real INF{std::numeric_limits<Real>::infinity()}; /*!< Infinity value. */
    static constexpr Real NAN{std::numeric_limits<Real>::quiet_NaN()}; /*!< Not-a-number value. */

    using Vector       = Eigen::Vector<Real, Eigen::Dynamic>;
    using Matrix       = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseVector = Eigen::Vector<Real, Eigen::Dynamic>;
    using SparseMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

    Vector       x;     // Primal point
    Real         rho;   // Penalty parameter value
    Real         rho_;  // Penalty parameter last value
    Real         mu;    // Interior-point parameter value
    Real         f;     // Objective function value (scaled)
    Real         fu;    // Objective function value (unscaled)
    Vector       g;     // Objective gradient value
    Vector       r1;    // Equality constraint slack value
    Vector       r2;    // Equality constraint slack value
    Vector       cE;    // Equality constraint value (scaled)
    SparseVector JE;    // Equality constraint Jacobian value
    Integer      JEnnz; // Equality constraint Jacobian nonzeros
    Vector       lE;    // Equality constraint multipliers
    Vector       s1;    // Inequality constraint slack value
    Vector       s2;    // Inequality constraint slack value
    Vector       cI;    // Inequality constraint value (scaled)
    SparseVector JI;    // Inequality constraint Jacobian value
    Integer      JInnz; // Inequality constraint Jacobian nonzeros
    Vector       lI;    // Inequality constraint multipliers
    SparseMatrix H;     // Hessian of Lagrangian
    Integer      Hnnz;  // Hessian of Lagrangian nonzeros
    Real         v;     // Feasibility violation measure value (scaled)
    Real         vu;    // Feasibility violation measure value (unscaled)
    Real         v0;    // Feasibility violation measure initial value
    Real         phi;   // Merit function value
    SparseMatrix AL;    // Newton matrix L-factor in LDL factorization
    SparseMatrix AD;    // Newton matrix D-factor in LDL factorization
    SparseMatrix AP;    // Newton matrix P-factor in LDL factorization
    SparseMatrix AS;    // Newton matrix S-factor in LDL factorization
    Integer      Annz;  // Newton matrix (upper triangle) nonzeros
    Real         shift; // Hessian shift value
    SparseVector b;     // Newton right-hand side
    Vecor        kkt;   // KKT errors
    Vecor        kkt_;  // KKT errors last value
    Integer      err;   // Function evaluation error flag

    Real         fs;      // Objective scaling factor
    Vector       cEs;     // Equality constraint scaling factors
    Vector       cEu;     // Equality constraint value (unscaled)
    Vector       cIs;     // Inequality constraint scaling factors
    Vector       cIu;     // Inequality constraint value (unscaled)
    SparseMatrix A;       // Newton matrix
    Integer      shift22; // Newton matrix (2,2)-block shift value
    Real         v_;      // Feasibility violation measure last value
    bool         cut_;    // Boolean value for last backtracking line search

    // Constructor
    Iterate(Parameter<Real> const & p, Input<Real> const & i, Counter const & c)
    {
      // Initialize quantities
      this->x     = i.x0;
      this->rho   = p.rho_init;
      this->mu    = p.mu_init;
      this->lE.setZero(i.nE);
      this->lI.setConstant(i.nI, 0.5);
      this->err   = 0;
      this->evalScalings(p, i, c);
      this->evalFunctions(i, c);
      this->evalGradients(i, c);
      this->evalDependent(p, i);
      this->v0    = 1;
      this->evalInfeasibility(i);
      this->v0    = this->v;
      this->evalInfeasibility(i);
      this->v_    = this->v;
      this->shift = 0;
      this->kkt_.setConstant(p.opt_err_mem, INF);
      this->cut_  = 0;
      this->evalHessian(i, c);
      this->Hnnz  = nnz(this->H);
      this->JEnnz = nnz(this->JE);
      this->JInnz = nnz(this->JI);
      this->initNewtonMatrix(i);
      this->evalNewtonMatrix(p, i, c);
    }

    // Termination checker
    Integer checkTermination(Parameter<Real> const & p, Input<Real> const & i, Counter const & c)
    {
      // Update termination based on optimality error of nonlinear optimization problem
      if (this->kkt(1) <= p.opt_err_tol && this->v <= p.opt_err_tol) {return 1;}

      // Update termination based on optimality error of feasibility problem
      if (this->kkt(0) <= p.opt_err_tol && this->v > p.opt_err_tol) {return 2;}

      // Update termination based on iteration count
      if (c.k >= p.iter_max) {return 3;}

      // Update termination based on invalid bounds
      if (i.vi > 0) {return 4;}

      // Update termination based on function evaluation error
      if (this->err > 0) {return 5;}

      return 0;
    }

    // Dependent quantity evaluator
    void evalDependent(Parameter<Real> const & p, Input<Real> const & i)
    {
      // Evaluate quantities dependent on penalty and interior-point parameters
      this->evalSlacks(p, i);
      this->evalMerit(i);
      this->evalKKTErrors(i);
    }

    // Function evaluator
    void evalFunctions(Input<Real> const & i, Counter const & c)
    {
      // Evaluate x in original space
      x_orig = this->evalXOriginal(i);

      // Initialize/Reset evaluation flag
      this->err = 0;

      // Increment function evaluation counter
      c.incrementFunctionCount();

      // Try AMPL functions evaluation
      try
      {
        // Evaluate AMPL functions
        i.f_orig(x_orig, this->f);
        i.c_orig(x_orig, c_orig);
      }
      catch (...)
      {
        // Set evaluation flag, default values, and return
        this->err = 1;
        this->f   = NAN;
        this->cE.setConstant(i.nE, NAN);
        this->cI.setConstant(i.nI, NAN);
        this->fu  = NAN;
        this->cEu.setConstan(i.nE, NAN);
        this->cIu.setConstan(i.nI, NAN);
        return;
      }

      // Set equality constraint values
      if (i.nE > 0) {this->cE = c_orig(i.I6) - i.b6;}

      // Initialize inequality constraint values
      if (i.nI > 0) {this->cI.setZero(i.nI);}

      // Set inequality constraint values
      if (i.n3 > 0) {
        this->cI(Eigen::seq(0, i.n3-1)) = i.l3 - this->x(Eigen::seq(i.n1, i.n1+i.n3-1));
      }
      if (i.n4 > 0) {
        this->cI(Eigen::seq(i.n3, i.n3+i.n4-1)) = -i.u4 + this->x(Eigen::seq(i.n1+i.n3, i.n1+i.n3+i.n4-1));
      }
      if (i.n5 > 0) {
        this->cI(Eigen::seq(i.n3+i.n4, i.n3+i.n4+i.n5-1)) = i.l5 - this->x(Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1));
        this->cI(Eigen::seq(i.n3+i.n4+i.n5, i.n3+i.n4+i.n5+i.n5-1)) = -i.u5 + this->x(Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1));
      }
      if (i.n7 > 0) {
        this->cI(Eigen::seq(i.n3+i.n4+i.n5+i.n5, i.n3+i.n4+i.n5+i.n5+i.n7-1)) = i.l7 - c_orig(i.I7);
      }
      if (i.n8 > 0) {
        this->cI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8-1)) = -i.u8 + c_orig(i.I8);
      }
      if (i.n9 > 0) {
        this->cI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9-1)) =  i.l9 - c_orig(i.I9);
        this->cI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9+i.n9-1)) = -i.u9 + c_orig(i.I9);
      }

      // Store unscaled quantities
      this->fu = this->f;
      if (i.nE > 0) {this->cEu = this->cE;}
      if (i.nI > 0) {this->cIu = this->cI;}

      // Scale quantities
      this->f *= this->fs;
      if (i.nE > 0) {this->cE = (this->cEs.array()*this->cE.array()).matrix();}
      if (i.nI > 0) {this->cI = (this->cIs.array()*this->cI.array()).matrix();}

    }

    // Gradient evaluator
    void evalGradients(Input<Real> const & i, Counter const & c)
    {
      // Evaluate x in original space
      x_orig = this->evalXOriginal(i);

      // Initialize/Reset evaluation flag
      this->err = 0;

      // Increment gradient evaluation counter
      c.incrementGradientCount;

      // Try AMPL gradients evaluation
      this->g.resize(i.nV);
      this->JE.resize(i.nE, i.nV);
      this->JI.resize(i.nI, i.nV);
      try
      {
        // Evaluate AMPL gradients
        i.g_orig(x_orig, g_orig);
        i.J_orig(x_orig, J_orig);
      }
      catch (...)
      {
        // Set evaluation flag, default values, and return
        this->err = 1;
        return;
      }

      // Set objective gradient
      this->g << g_orig(i.I1); g_orig(i.I3); g_orig(i.I4); g_orig(i.I5);

      // Set equality constraint Jacobian
      if (i.nE > 0) {this->JE = << J_orig(i.I6, i.I1) J_orig(i.I6, i.I3) J_orig(i.I6, i.I4) J_orig(i.I6, i.I5);}

      // Initialize inequality constraint Jacobian
      if (i.nI > 0) {this->JI.resize(i.nI, i.nV);}

      // Set inequality constraint Jacobian
      if (i.n3 > 0) {
        this->JI(Eigen::seq(0, i.n3-1), Eigen::seq(i.n1, i.n1+i.n3-1)) = -speye(i.n3);
      }
      if (i.n4 > 0) {
        this->JI(Eigen::seq(i.n3, i.n3+i.n4-1), Eigen::seq(i.n1+i.n3, i.n1+i.n3+i.n4-1)) = speye(i.n4);
      }
      if (i.n5 > 0) {
        this->JI(Eigen::seq(i.n3+i.n4, i.n3+i.n4+i.n5-1), Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1)) = -speye(i.n5);
        this->JI(Eigen::seq(i.n3+i.n4+i.n5, i.n3+i.n4+i.n5+i.n5-1), Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1)) = speye(i.n5);
      }
      if (i.n7 > 0) {
        this->JI(Eigen::seq(i.n3+i.n4+i.n5+i.n5, i.n3+i.n4+i.n5+i.n5+i.n7-1), Eigen::seq(0, i.n1+i.n3+i.n4+i.n5-1)) <<
          -J_orig(i.I7, i.I1) J_orig(i.I7,i.I3) J_orig(i.I7,i.I4) J_orig(i.I7,i.I5);
      }
      if (i.n8 > 0) {
        this->JI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8-1), Eigen::seq(0, i.n1+i.n3+i.n4+i.n5-1)) <<
          J_orig(i.I8,i.I1) J_orig(i.I8,i.I3) J_orig(i.I8,i.I4) J_orig(i.I8,i.I5);
      }
      if (i.n9 > 0) {
        this->JI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9-1), Eigen::seq(0, i.n1+i.n3+i.n4+i.n5-1)) <<
          -J_orig(i.I9,i.I1) J_orig(i.I9,i.I3) J_orig(i.I9,i.I4) J_orig(i.I9,i.I5);
        this->JI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9+i.n9-1), Eigen::seq(0, i.n1+i.n3+i.n4+i.n5-1)) <<
          J_orig(i.I9,i.I1) J_orig(i.I9,i.I3) J_orig(i.I9,i.I4) J_orig(i.I9,i.I5);
      }

      // Scale objective gradient
      this->g *= this->fs;

      // Scale constraint Jacobians
      if (i.nE > 0) {
        for (Integer k{0}; k < this->JE.outerSize(); ++k) {
          for (SparseMatrix::InnerIterator it(this->JE, k); it; ++it) {
            it.valueRef() *= this->cEs(it.row());
          }
        }
      }
      if (i.nI > 0) {
        for (Integer k{0}; k < this->JI.outerSize(); ++k) {
          for (SparseMatrix::InnerIterator it(this->JI, k); it; ++it) {
            it.valueRef() *= this->cEs(it.row());
          }
        }
      }
    }

    // Hessian evaluator
    void evalHessian(Input<Real> const & i, Counter const & c)
    {
      // Evaluate lambda in original space
      Vector l_orig, x_orig;
      this->evalLambdaOriginal(i, l_orig);
      this->evalXOriginal(i, x_orig);

      // Initialize/Reset evaluation flag
      this->err = 0;

      // Increment Hessian evaluation counter
      c.incrementHessianCount();

      // Try AMPL Hessian evaluation
      try
      {
        // Evaluate H_orig
        if (i.nE+i.n7+i.n8+i.n9 == 0) {i.H_orig(x_orig, H_orig);}
        else {i.H_orig(x_orig, l_orig, H_orig);}
      }
      catch (...)
      {
        // Set evaluation flag, default values, and return
        this->err = 1;
        this->H.resize(i.nV, i.nV);
        return;
      }

      // Set Hessian of the Lagrangian
      this->H <<
        H_orig(i.I1, i.I1), H_orig(i.I1, i.I3), H_orig(i.I1, i.I4), H_orig(i.I1, i.I5),
        H_orig(i.I3, i.I1), H_orig(i.I3, i.I3), H_orig(i.I3, i.I4), H_orig(i.I3, i.I5),
        H_orig(i.I4, i.I1), H_orig(i.I4, i.I3), H_orig(i.I4, i.I4), H_orig(i.I4, i.I5),
        H_orig(i.I5, i.I1), H_orig(i.I5, i.I3), H_orig(i.I5, i.I4), H_orig(i.I5, i.I5);

      // Rescale H
      this->H *= this->rho*this->fs;
    }

    // Infeasibility evaluator
    void evalInfeasibility(Input<Real> const & i)
    {
      // Evaluate scaled and unscaled feasibility violations
      this->v  = this->evalViolation(i,this->cE ,this->cI) / std::max(1,this->v0);
      this->vu = this->evalViolation(i,this->cEu,this->cIu);
    }

    // KKT error evaluator
    Real evalKKTError(Input<Real> const & i, Real const rho, Real const mu)
    {
      // Initialize optimality vector
      Vector kkt(i.nV+2*i.nE+2*i.nI);
      kkt.setZero();

      // Set gradient of penalty objective
      kkt(Eigen::seq(0, i.nV-1)) = rho*this->g;

      // Set gradient of Lagrangian for constraints
      if (i.nE > 0) {kkt(Eigen::seq(0, i.nV-1)) = kkt(Eigen::seq(0, i.nV-1)) + (this->lE.transpose()*this->JE).transpose();} // OPTMIZE
      if (i.nI > 0) {kkt(Eigen::seq(0, i.nV-1)) = kkt(Eigen::seq(0, i.nV-1)) + (this->lI.transpose()*this->JI).transpose();} // OPTMIZE

      // Set complementarity for constraint slacks
      if (i.nE > 0) {
        kkt(Eigen::seq(i.nV, i.nV+2*i.nE-1)) <<
          (this->r1.array() * (1 + this->lE).array()).matrix() - mu,
          (this->r2.array() * (1 - this->lE).array()).matrix() - mu;
      }
      if (i.nI > 0) {
        kkt(Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1)) <<
          (this->s1.array() * (0 + this->lI).array()).matrix() - mu,
          (this->s2.array() * (1 - this->lI).array()).matrix() - mu;
      }

      // Scale complementarity
      if (rho > 0) {kkt = (1.0 / std::max(1, (rho*this->g).array().abs().maxCoeff()))*kkt.maxCoeff();}

      // Evaluate optimality error
      return kkt.array().abs().maxCoeff();
    }

    // KKT errors evaluator
    void evalKKTErrors(Input<Real> const & i)
    {
      // Loop to compute optimality errors
      this->kkt(1) = this->evalKKTError(i, 0, 0);
      this->kkt(2) = this->evalKKTError(i, this->rho, 0);
      this->kkt(3) = this->evalKKTError(i, this->rho, this->mu);
    }

    // Evaluator of lambda in original space
    void evalLambdaOriginal(Input<Real> const & i, Vector & l)
    {
      // Initialize multipliers in original space
      l.setZero(i.nE+i.n7+i.n8+i.n9);

      // Scale equality constraint multipliers
      if (i.nE > 0) {lE = (this->lE.array()*(this->cEs/(this->rho*this->fs)).array()).matrix();}

      // Set equality constraint multipliers in original space
      if (i.nE > 0) {l(i.I6) = lE;}

      // Scale inequality constraint multipliers
      if (i.n7+i.n8+i.n9 > 0) {lI = (this->lI.array()*(this->cIs/(this->rho*this->fs)).array()).matrix();}

      // Set inequality constraint multipliers in original space
      if (i.n7 > 0) {
        l(i.I7) = -lI(Eigen::seq(i.n3+i.n4+i.n5+i.n5, i.n3+i.n4+i.n5+i.n5+i.n7-1));
      }
      if (i.n8 > 0) {
        l(i.I8) = +lI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8-1));
      }
      if (i.n9 > 0) {
        l(i.I9) = -lI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9-1))
                  +lI(Eigen::seq(i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9, i.n3+i.n4+i.n5+i.n5+i.n7+i.n8+i.n9+i.n9-1));
      }
    }

    // Matrices evaluator
    void evalMatrices(Parameter<Real> const & p, Input<Real> const & i, Counter const & c)
    {
      // Evaluate Hessian and Newton matrices
      this->evalHessian(i, c);
      this->evalNewtonMatrix(p, i, c);
    }

    // Merit evaluator
    void evalMerit(Input<Real> const & i)
    {
      // Initialize merit for objective
      this->phi = this->rho*this->f;

      // Update merit for slacks
      if (i.nE > 0) {
        Vector r_all(this->r1.size() + this->r2.size()); r_all << this->r1, this->r2;
        this->phi -= this->mu * r_all.array().log().sum() + r_all.sum();
      }
      if (i.nI > 0) {
        Eigen::VectorXd s_all(this->s1.size() + this->s2.size()); s_all << this->s1, this->s2;
        this->phi -= this->mu * s_all.array().log().sum() + this->s2.sum();
      }
    }

    // Newton matrix evaluator
    void evalNewtonMatrix(Parameter<Real> const & p, Input<Real> const & i, Counter const & c)
    {
      // Check for equality constraints
      if (i.nE > 0)
      {
        // Set diagonal terms
        for (Integer j{0}; j < i.nE, ++j) {
          this->A(i.nV+     j,i.nV+     j) = (1+this->lE(j))/this->r1(j);
          this->A(i.nV+i.nE+j,i.nV+i.nE+j) = (1-this->lE(j))/this->r2(j);
        }

        // Set constraint Jacobian
        this->A(1+i.nV+2*i.nE+2*i.nI:i.nV+3*i.nE+2*i.nI, Eigen::seq(0, i.nV-1)) = this->JE;
      }

      // Check for inequality constraints
      if (i.nI > 0)
      {
        // Set diagonal terms
        for (Integer j{0}; j < i.nI, ++j) {
          this->A(i.nV+2*i.nE+j, i.nV+2*i.nE+j)           = (0+this->lI(j))/this->s1(j);
          this->A(i.nV+2*i.nE+i.nI+j, i.nV+2*i.nE+i.nI+j) = (1-this->lI(j))/this->s2(j);
        }

        // Set constraint Jacobian
        this->A(1+i.nV+3*i.nE+2*i.nI:i.nV+3*i.nE+3*i.nI, Eigen::seq(0, i.nV-1)) = this->JI;
      }

      // Set minimum potential shift
      min_shift = std::max(p.shift_min, p.shift_factor1*this->shift);

      // Initialize Hessian modification
      if (this->cut_ == 1) {this->shift = std::min(p.shift_max,min_shift/p.shift_factor2);
      else {this->shift = 0;}

      // Initialize inertia correction loop
      bool done{false};
      this->shift22 = 0;

      // Loop until inertia is correct
      while (!done && this->shift < p.shift_max)
      {
        // Set Hessian of Lagrangian
        this->A(Eigen::seq(0, i.nV-1), Eigen::seq(0, i.nV-1)) = this->H+this->shift*speye(i.nV);

        // Set diagonal terms
        for (Integer j{0}; j < i.nE; ++j) {
          this->A(i.nV+2*i.nE+2*i.nI+j, i.nV+2*i.nE+2*i.nI+j) = -this->shift22;
        }

        // Set diagonal terms
        for (Integer j{0}; j < i.nI; ++j) {
          this->A(i.nV+3*i.nE+2*i.nI+j, i.nV+3*i.nE+2*i.nI+j) = -this->shift22;
        }

        // Set number of nonzeros in (upper triangle of) Newton matrix
        this->Annz = nnz<MatrixView::TRIL>(this->A);

        // Factor primal-dual matrix
        [this->AL, this->AD, this->AP, this->AS, neig] = ldl(tril(this->A), p.pivot_thresh, 'vector');

        // Increment factorization counter
        c.incrementFactorizationCount();

        // Set number of nonnegative eigenvalues
        peig = i.nA - neig;

        // Check inertia
        if (peig < i.nV+2*i.nE+2*i.nI) {this->shift = std::max(min_shift,this->shift/p.shift_factor2);}
        else if (neig < i.nE+i.nI & this->shift22 == 0) {this->shift22 = p.shift_min;}
        else {done = true;}
      }

      // Update Hessian
      this->H = this->H+this->shift*speye(i.nV);
    }

    // Newton right-hand side evaluator
    void evalNewtonRhs(Input<Real> const & i)
    {
      // Initialize right-hand side vector
      this->b.setZero(i.nA);

      // Set gradient of objective
      this->b(Eigen::seq(0, i.nV-1)) = this->rho*this->g;

      // Set gradient of Lagrangian for constraints
      if (i.nE > 0) {this->b(Eigen::seq(0, i.nV-1)) = this->b(Eigen::seq(0, i.nV-1)) + (this->lE.transpose()*this->JE).transpose();}
      if (i.nI > 0) {this->b(Eigen::seq(0, i.nV-1)) = this->b(Eigen::seq(0, i.nV-1)) + (this->lI.transpose()*this->JI).transpose();}

      // Set complementarity for constraint slacks
      if (i.nE > 0) {this->b(Eigen::seq(i.nV, i.nV+2*i.nE-1)) <<
        1 + this->lE - this->mu/this->r1, 1 - this->lE - this->mu/this->r2;
      }
      if (i.nI > 0) {this->b(Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1)) <<
        this->lI - this->mu/this->s1, 1 - this->lI - this->mu/this->s2;
      }

      // Set penalty-interior-point constraint values
      if (i.nE > 0) {this->b(Eigen::seq(i.nV+2*i.nE+2*i.nI, i.nV+3*i.nE+2*i.nI-1)) = this->cE + this->r1 - this->r2;}
      if (i.nI > 0) {this->b(Eigen::seq(i.nV+3*i.nE+2*i.nI, i.nV+3*i.nE+3*i.nI-1)) = this->cI + this->s1 - this->s2;}
    }

    // Scalings evaluator
    void evalScalings(Parameter<Real> const & p, Input<Real> const & i, Counter const & c)
    {
      // Initialize scalings
      this->fs = 1;
      this->cEs.setOnes(i.nE);
      this->cIs.setOnes(i.nI);

      // Evaluate gradients
      this->evalGradients(i,c);

      // Scale down objective if norm of gradient is too large
      this->fs = p.grad_max / std::max(this->g.array().abs().maxCoeff(), p.grad_max);

      // Loop through equality constraints
      for (Integer j{0}; j < i.nE; ++j)
      {
        // Scale down equality constraint j if norm of gradient is too large
        this->cEs(j) = p.grad_max / std::max(this->JE.row(j).array().abs().maxCoeff(), p.grad_max);
      }

      // Loop through inequality constraints
      for (Integer j{0}; j < i.nI; ++j)
      {
        // Scale down inequality constraint j if norm of gradient is too large
        this->cIs(j) = p.grad_max / std::max(this->JI.row(j).array().abs().maxCoeff(), p.grad_max);
      }
    }

    // Slacks evaluator
    void evalSlacks(Parameter<Real> const & p, Input<Real> const & i)
    {
      // Check for equality constraints
      if (i.nE > 0)
      {
        // Set slacks
        this->r1 = 0.5*(this->mu - this->cE + std::sqrt(this->cE.array().square().matrix() + this->mu*this->mu));
        this->r2 = 0.5*(this->mu + this->cE + std::sqrt(this->cE.array().square().matrix() + this->mu*this->mu));

        // Adjust for numerical error
        this->r1 = std::max(this->r1, p.slack_min);
        this->r2 = std::max(this->r2, p.slack_min);
      }

      // Check for inequality constraints
      if (i.nI > 0)
      {
        // Set slacks
        this->s1 = 0.5*(2.0*this->mu - this->cI + std::sqrt(this->cI.array().square().matrix() + 4.0*this->mu*this->mu));
        this->s2 = 0.5*(2.0*this->mu + this->cI + std::sqrt(this->cI.array().square().matrix() + 4.0*this->mu*this->mu));

        // Adjust for numerical error
        this->s1 = std::max(this->s1, p.slack_min);
        this->s2 = std::max(this->s2, p.slack_min);
      }
    }

    // Evaluator of x in original space
    void evalXOriginal(Input<Real> const & i, Vector & x)
    {
      // Initialize x in original space
      x.setZero(i.n0);

      // Evaluate x in original space
      x(i.I1) = this->x(Eigen::seq(0, i.n1-1));
      x(i.I2) = i.b2;
      x(i.I3) = this->x(Eigen::seq(i.n1, i.n1+i.n3-1));
      x(i.I4) = this->x(Eigen::seq(i.n1+i.n3, i.n1+i.n3+i.n4-1));
      x(i.I5) = this->x(Eigen::seq(i.n1+i.n3+i.n4, i.n1+i.n3+i.n4+i.n5-1));
    }

    // Gets primal-dual point
    void getSolution(Input<Real> const & i, Vector & x, Vector & l) const
    {
      this->evalXOriginal(i, x);
      this->evalLambdaOriginal(i, l);
    }

    // Initializes Newton matrix
    void initNewtonMatrix(Input<Real> const & i)
    {
      // Allocate memory
      this->A = spalloc(i.nA,i.nA,this->Hnnz+5*i.nE+5*i.nI+this->JEnnz+this->JInnz);

      // Initialize interior-point Hessians
      this->A(Eigen::seq(i.nV, i.nV+2*i.nE-1), Eigen::seq(i.nV, i.nV+2*i.nE-1)) = speye(2*i.nE);
      this->A(Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1), Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1)) = speye(2*i.nI);

      // Check for equality constraints
      if (i.nE > 0)
      {
        // Initialize constraint Jacobian
        this->A(Eigen::seq(i.nV+2*i.nE+2*i.nI, i.nV+3*i.nE+2*i.nI-1), Eigen::seq(i.nV, i.nV+i.nE-1)) =  speye(i.nE);
        this->A(Eigen::seq(i.nV+2*i.nE+2*i.nI, i.nV+3*i.nE+2*i.nI-1), Eigen::seq(i.nV+i.nE, i.nV+2*i.nE-1)) = -speye(i.nE);
      }

      // Check for inequality constraints
      if (i.nI > 0)
      {
        // Initialize constraint Jacobian
        this->A(Eigen::seq(i.nV+3*i.nE+2*i.nI, i.nV+3*i.nE+3*i.nI-1), Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+i.nI-1)) =  speye(i.nI);
        this->A(Eigen::seq(i.nV+3*i.nE+2*i.nI, i.nV+3*i.nE+3*i.nI-1), Eigen::seq(i.nV+2*i.nE+i.nI, i.nV+2*i.nE+2*i.nI-1)) = -speye(i.nI);
      }
    }

    // Set interior-point parameter
    void setMu(Real const mu) {this->mu = mu;}

    // Set primal variables
    void setPrimals(Input<Real> const & i, Vector const & x, Vector const & r1, Vector const & r2,
      Vector const & s1, Vector const & s2, Vector const & lE, Vector const & lI, Real const f,
      Vector const & cE, Vector const & cI, Real const phi)
    {
      // Set primal variables
      this->x = x; this->f = f;
      if (i.nE > 0) {this->cE = cE; this->r1 = r1; this->r2 = r2; this->lE = lE;}
      if (i.nI > 0) {this->cI = cI; this->s1 = s1; this->s2 = s2; this->lI = lI;}
      this->phi = phi;
    }

    // Set penalty parameter
    void setRho(Real const rho) {this->rho = rho;}

    // Set last penalty parameter
    void setRhoLast(Real const rho) {this->rho_ = rho;}

    // Iterate updater
    void updateIterate(Parameter<Real> const & p, Input<Real> const & i, Counter const & c,
      Direction<Real> const & d, Acceptance<Real> const & a)
    {
      // Update last quantities
      this->v_   = this->v;
      this->cut_ = (a.p < a.p0);

      // Update iterate quantities
      this->updatePoint(i, d, a);
      this->evalInfeasibility(i);
      this->evalGradients(i, c);
      this->evalDependent(p, i);

      // Update last KKT errors
      this->kkt_.resize(p.opt_err_mem);
      this->kkt_ << this->kkt(1), this->kkt_(Eigen::seq(0, p.opt_err_mem-2));
    }

    // Parameter updater
    void updateParameters(Parameter<Real> const & p, Input<Real> const & i)
    {
      // Check for interior-point parameter update based on optimality error
      while (this->mu > p.mu_min && this->kkt(2) <= std::max({this->mu, p.opt_err_tol-this->mu}))
      {
        // Restrict interior-point parameter increase
        p.setMuMaxExpZero();

        // Update interior-point parameter
        if this->mu > p.mu_min
        {
          // Decrease interior-point
          this->mu = std::max(p.mu_min, std::min(p.mu_factor*this->mu, std::pow(this->mu, p.mu_factor_exp)));

          // Evaluate penalty and interior-point parameter dependent quantities
          this->evalDependent(p, i);
        }
      }

      // Check for penalty parameter update based on optimality error
      if (this->kkt(1) <= p.opt_err_tol && this->v > p.opt_err_tol) || this->v > std::max({1.0, this->v_, p.infeas_max})
      {
        // Update penalty parameter
        if this->rho > p.rho_min
        {
          // Decrease penalty parameter
          this->rho = std::max(p.rho_min, p.rho_factor*this->rho);

          // Evaluate penalty and interior-point parameter dependent quantities
          this->evalDependent(p, i);
        }
      }
    }

    // Primal point updater
    void updatePoint(Input<Real> const & i, Direction<Real> const & d, Acceptance<Real> const & a)
    {
      // Update primal and dual variables
      this->x += a.p*d.x ;
      if (i.nE > 0) {this->r1 += a.p*d.r1; this->r2 += a.p*d.r2;}
      if (i.nI > 0) {this->s1 += a.p*d.s1; this->s2 += a.p*d.s2;}
      if (i.nE > 0) {this->lE += a.d*d.lE;}
      if (i.nI > 0) {this->lI += a.d*d.lI;}
    }

    // Feasibility violation evaluator
    static Real evalViolation(Input<Real> const & i, Vector const & cE, Vector const & cI)
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
      return vec.lpNorm<1>();
    }

  }; // struct Iterate

} // namespace Pipal

#endif /* INCLUDE_PIPAL_ITERATE_HH */
