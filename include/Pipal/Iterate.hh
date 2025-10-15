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
#include "Pipal.hh"
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

    using Vector = Eigen::Vector<Real, Eigen::Dynamic>;
    using Mask   = Eigen::Array<bool, Eigen::Dynamic, 1>;

    x     /*!< Primal point. */
    rho   /*!< Penalty parameter value. */
    rho_  /*!< Penalty parameter last value. */
    mu    /*!< Interior-point parameter value. */
    f     /*!< Objective function value (scaled). */
    fu    /*!< Objective function value (unscaled). */
    g     /*!< Objective gradient value. */
    r1    /*!< Equality constraint slack value. */
    r2    /*!< Equality constraint slack value. */
    cE    /*!< Equality constraint value (scaled). */
    JE    /*!< Equality constraint Jacobian value. */
    JEnnz /*!< Equality constraint Jacobian nonzeros. */
    lE    /*!< Equality constraint multipliers. */
    s1    /*!< Inequality constraint slack value. */
    s2    /*!< Inequality constraint slack value. */
    cI    /*!< Inequality constraint value (scaled). */
    JI    /*!< Inequality constraint Jacobian value. */
    JInnz /*!< Inequality constraint Jacobian nonzeros. */
    lI    /*!< Inequality constraint multipliers. */
    H     /*!< Hessian of Lagrangian. */
    Hnnz  /*!< Hessian of Lagrangian nonzeros. */
    v     /*!< Feasibility violation measure value (scaled). */
    vu    /*!< Feasibility violation measure value (unscaled). */
    v0    /*!< Feasibility violation measure initial value. */
    phi   /*!< Merit function value. */
    AL    /*!< Newton matrix L-factor in LDL factorization. */
    AD    /*!< Newton matrix D-factor in LDL factorization. */
    AP    /*!< Newton matrix P-factor in LDL factorization. */
    AS    /*!< Newton matrix S-factor in LDL factorization. */
    Annz  /*!< Newton matrix (upper triangle) nonzeros. */
    shift /*!< Hessian shift value. */
    b     /*!< Newton right-hand side. */
    kkt   /*!< KKT errors. */
    kkt_  /*!< KKT errors last value. */
    err   /*!< Function evaluation error flag. */

    fs      /*!< Objective scaling factor. */
    cEs     /*!< Equality constraint scaling factors. */
    cEu     /*!< Equality constraint value (unscaled). */
    cIs     /*!< Inequality constraint scaling factors. */
    cIu     /*!< Inequality constraint value (unscaled). */
    A       /*!< Newton matrix. */
    shift22 /*!< Newton matrix (2,2)-block shift value. */
    v_      /*!< Feasibility violation measure last value. */
    cut_    /*!< Boolean value for last backtracking line search. */

    /**
     * \brief Fill the iterate structure with problem data.
     * \param[in] problem Problem object.
     * \param[in] parameter Parameter structure.
     * \param[in] input Input structure.
     * \param[in] counter Counter structure.
     */
    void fill_iterate(Problem<Real> const & problem, Parameter<Real> const & parameter,
      Input<Real, Integer> const & input, Counter<Integer> const & counter)
    {
      #define CMD "Pipal::Solver::fill_iterate(...): "

      this->x  = input.x0;
      this->rho = parameter.rho_init;
      this->mu  = parameter.mu_init;
      this->lE  = zeros(input.nE, 1);
      this->lI  = (1/2)*ones(input.nI, 1);
      this->err = 0;
      this->evalScalings(problem, input, counter);
      this->evalFunctions(input, counter);
      this->eval_gradients(input, counter);
      this->evalDependent(problem, input);
      this->v0 = 1;
      this->evalInfeasibility(input);
      this->v0 = this->v;
      this->evalInfeasibility(input);
      this->v_ = this->v;
      this->shift = 0;
      this->kkt_ = inf*ones(parameter.opt_err_mem, 1);
      this->cut_ = 0;
      this->evalHessian(input, counter);
      this->Hnnz = nnz(this->H);
      this->JEnnz = nnz(this->JE);
      this->JInnz = nnz(this->JI);
      this->initNewtonMatrix(input);
      this->evalNewtonMatrix(problem, input, counter);

      #undef CMD
    }

    // Termination checker
    Integer check_termination(Parameter<Real> const & parameter, Input<Real, Integer> const & input,
      Counter<Integer> const & counter)
    {
      // Update termination based on optimality error of nonlinear optimization problem
      if (this->kkt(1) <= parameter.opt_err_tol && this->v <= parameter.opt_err_tol) {return 1;}

      // Update termination based on optimality error of feasibility problem
      if (this->kkt(0) <= parameter.opt_err_tol && this->v > parameter.opt_err_tol) {return 2;}

      // Update termination based on iteration count
      if (counter.k >= parameter.iter_max) {return 3;}

      // Update termination based on invalid bounds
      if (input.vi > 0) {return 4;}

      // Update termination based on function evaluation error
      if (this->err > 0) {return 5;}

      return 0;
    }

    // Dependent quantity evaluator
    void eval_dependent(Parameter<Real> const & parameter, Input<Real, Integer> const & input)
    {
      // Evaluate quantities dependent on penalty and interior-point parameters
      this->eval_slacks(parameter, input);
      this->eval_merit(parameter, input);
      this->eval_kkt_errors(parameter, input);
    }

    // Function evaluator
    bool eval_functions(Problem<Real> const & problem, Input<Real, Integer> const & input, Counter<Integer> & counter)
    {
      #define CMD "Pipal::Iterate::eval_functions(...): "

      // Evaluate x in original space
      Vector x_orig(this->eval_x_original(input));

      // Initialize/Reset evaluation flag
      this->err = 0;

      // Try AMPL functions evaluation
      ++counter.function_evaluations;
      Vector c_orig;
      try
      {
        // Evaluate AMPL functions
        PIPAL_ASSERT(problem.objective(x_orig, this->f),
          CMD "error in evaluating objective function");
        PIPAL_ASSERT(problem.constraints(x_orig, c_orig),
          CMD "error in evaluating constraint functions");
      }
      catch (...)
      {
        // Set evaluation flag, default values, and return
        this->err = 1;
        // FIXME: this->f   = NAN;
        // FIXME: this->fu  = NAN;
        // FIXME: this->cE.setConstant(input.nE, NAN);
        // FIXME: this->cInput.setConstant(input.nI, NAN);
        // FIXME: this->cEu.setConstant(input.nE, NAN);
        // FIXME: this->cIu.setConstant(input.nI, NAN);
        return false;
      }

      // Set equality constraint values
      if (input.nE > 0) {this->cE = c_orig(input.i6) - input.b6;}

      // Initialize inequality constraint values
      if (input.nI > 0) {this->cInput.setZero(input.nI);}

      // Set inequality constraint values
      Integer idx_c_ini{0}, idx_c_end{input.n3-1}, idx_x_ini{input.n1}, idx_x_end{input.n1+input.n3-1};
      if (input.n3 > 0) {
        this->cI(Eigen::seq(idx_c_ini, idx_c_end)) = input.l3 - this->x(Eigen::seq(idx_x_ini, idx_x_end));
      }
      idx_c_ini += input.n3; idx_c_end += input.n4; idx_x_ini += input.n3; idx_x_end += input.n3;
      if (input.n4 > 0) {
        this->cI(Eigen::seq(idx_c_ini, idx_c_end)) = -input.u4 + this->x(Eigen::seq(idx_x_ini, idx_x_end));
      }
      idx_c_ini += input.n4; idx_c_end += input.n5; idx_x_ini += input.n4; idx_x_end += input.n4;
      if (input.n5 > 0) {
        Vector x_seq(this->x(Eigen::seq(idx_x_ini, idx_x_end)));
        this->cI(Eigen::seq(idx_c_ini, idx_c_end)) =  input.l5 - x_seq;
        this->cI(Eigen::seq(idx_c_ini+input.n5, idx_c_end+input.n5)) = -input.u5 + x_seq;
      }
      idx_c_ini += 2*input.n5; idx_c_end += (input.n5 + input.n7);
      if (input.n7 > 0) {
        this->cI(Eigen::seq(idx_c_ini, idx_c_end)) = input.l7 - c_orig(input.i7);
      }
      idx_c_ini += input.n7; idx_c_end += input.n8;
      if (input.n8 > 0) {
        this->cI(Eigen::seq(idx_c_ini, idx_c_end)) = -input.u8 + c_orig(input.i8);
      }
      idx_c_ini += input.n8; idx_c_end += (input.n8 + input.n9);
      if (input.n9 > 0) {
        this->cI(Eigen::seq(idx_c_ini, idx_c_end))                   =  input.l9 - c_orig(input.i9);
        this->cI(Eigen::seq(idx_c_ini+input.n9, idx_c_end+input.n9)) = -input.u9 + c_orig(input.i9);
      }

      // Store unscaled quantities
      this->fu = this->f;
      if (input.nE > 0) {this->cEu = this->cE;}
      if (input.nI > 0) {this->cIu = this->cI;}

      // Scale quantities
      this->f *= this->fs;
      if (input.nE > 0) {this->cE = (this->cEs.array() * this->cE.array()).matrix();}
      if (input.nI > 0) {this->cI = (this->cIs.array() * this->cInput.array()).matrix();}

      return true;

      #undef CMD
    }

    // Gradient evaluator
    bool eval_gradients(Problem<Real> const & problem, Input<Real, Integer> const & input, Counter<Integer> & counter)
    {
      #define CMD "Pipal::Iterate::eval_gradients(...): "

      // Evaluate x in original space
      Vector x_orig(this->eval_x_original(input));

      // Initialize/Reset evaluation flag
      this->err = 0;

      // Try AMPL gradients evaluation
      ++counter.gradient_evaluations;
      Vector g_orig;
      Matrix J_orig;
      try
      {
        // Evaluate AMPL gradients
        problem.g_orig(x_orig, g_orig);
        problem.J_orig(x_orig, J_orig);
      }
      catch (...)
      {
        // Set evaluation flag, default values, and return
        this->err = 1;
        // FIXME: this->g  = sparse(input.nV,1);
        // FIXME: this->JE = sparse(input.nE,input.nV);
        // FIXME: this->JI = sparse(input.nI,input.nV);
        return false;
      }

      // Set objective gradient
      this->g.resize(input.nV);
      this->g << g_orig(input.i1), g_orig(input.i3), g_orig(input.i4), g_orig(input.i5);

      // Set equality constraint Jacobian
      if (input.nE > 0) {this->JE = [J_orig(input.i6, input.i1) J_orig(input.i6, input.i3) J_orig(input.i6, input.i4) J_orig(input.i6, input.i5)];}

      // Initialize inequality constraint Jacobian
      if (input.nI > 0) {this->JI = sparse(input.nI, input.nV)};

      // Set inequality constraint Jacobian
      if (input.n3 > 0)
        this->JI(1:input.n3,1+input.n1:input.n1+input.n3) = -speye(input.n3); end;
      if (input.n4 > 0)
        this->JI(1+input.n3:input.n3+input.n4,1+input.n1+input.n3:input.n1+input.n3+input.n4) =  speye(input.n4); end;
      if (input.n5 > 0)
        this->JI(1+input.n3+input.n4:input.n3+input.n4+input.n5,1+input.n1+input.n3+input.n4:input.n1+input.n3+input.n4+input.n5) = -speye(input.n5);
        this->JI(1+input.n3+input.n4+input.n5:input.n3+input.n4+input.n5+input.n5,1+input.n1+input.n3+input.n4:input.n1+input.n3+input.n4+input.n5) =  speye(input.n5); end;
      if (input.n7 > 0)
        this->JI(1+input.n3+input.n4+input.n5+input.n5:input.n3+input.n4+input.n5+input.n5+input.n7,1:input.n1+input.n3+input.n4+input.n5) = -[J_orig(input.i7,input.i1) J_orig(input.i7,input.i3) J_orig(input.i7,input.i4) J_orig(input.i7,input.i5)]; end;
      if (input.n8 > 0)
        this->JI(1+input.n3+input.n4+input.n5+input.n5+input.n7:input.n3+input.n4+input.n5+input.n5+input.n7+input.n8,1:input.n1+input.n3+input.n4+input.n5) =  [J_orig(input.i8,input.i1) J_orig(input.i8,input.i3) J_orig(input.i8,input.i4) J_orig(input.i8,input.i5)]; end;
      if (input.n9 > 0)
        this->JI(1+input.n3+input.n4+input.n5+input.n5+input.n7+input.n8:input.n3+input.n4+input.n5+input.n5+input.n7+input.n8+input.n9,1:input.n1+input.n3+input.n4+input.n5) = -[J_orig(input.i9,input.i1) J_orig(input.i9,input.i3) J_orig(input.i9,input.i4) J_orig(input.i9,input.i5)];
        this->JI(1+input.n3+input.n4+input.n5+input.n5+input.n7+input.n8+input.n9:input.n3+input.n4+input.n5+input.n5+input.n7+input.n8+input.n9+input.n9,1:input.n1+input.n3+input.n4+input.n5) =  [J_orig(input.i9,input.i1) J_orig(input.i9,input.i3) J_orig(input.i9,input.i4) J_orig(input.i9,input.i5)]; end;

      // Scale objective gradient
      this->g *= this->fs;

      // Scale constraint Jacobians
      if (input.nE > 0) {this->JE = spdiags(this->cEs,0:0,input.nE,input.nE)*this->JE;}
      if (input.nI > 0) {this->JI = spdiags(this->cIs,0:0,input.nI,input.nI)*this->JI;}

      return true;

      #undef CMD
    }

//            // Hessian evaluator
//            function evalHessian(z,i,c)
//
//              // Evaluate lambda in original space
//              l_orig = this->evalLambdaOriginal(i);
//              x_orig = this->eval_x_original(i);
//
//              // Initialize/Reset evaluation flag
//              this->err = 0;
//
//              // Increment Hessian evaluation counter
//              counter.incrementHessianCount;
//
//              // Try AMPL Hessian evaluation
//              try
//
//                // Evaluate H_orig
//                if (input.nE+input.n7+input.n8+input.n9 == 0), H_orig = input.H_orig(x_orig, []);
//                else                           H_orig = input.H_orig(x_orig, l_orig); end;
//
//              catch
//
//                // Set evaluation flag, default values, and return
//                this->err = 1; this->H = sparse(input.nV,input.nV); return;
//
//              end
//
//              // Set Hessian of the Lagrangian
//              this->H = [H_orig(input.i1,input.i1) H_orig(input.i1,input.i3) H_orig(input.i1,input.i4) H_orig(input.i1,input.i5);
//                     H_orig(input.i3,input.i1) H_orig(input.i3,input.i3) H_orig(input.i3,input.i4) H_orig(input.i3,input.i5);
//                     H_orig(input.i4,input.i1) H_orig(input.i4,input.i3) H_orig(input.i4,input.i4) H_orig(input.i4,input.i5);
//                     H_orig(input.i5,input.i1) H_orig(input.i5,input.i3) H_orig(input.i5,input.i4) H_orig(input.i5,input.i5);];
//
//              // Rescale H
//              this->H = this->rho*this->fs*this->H;
//
//            end
//
//            // Infeasibility evaluator
//            function evalInfeasibility(z,i)
//
//              // Evaluate scaled and unscaled feasibility violations
//              this->v  = this->evalViolation(i,this->cE ,this->cI )/max(1,this->v0);
//              this->vu = this->evalViolation(i,this->cEu,this->cIu);
//
//            end
//
//            // KKT error evaluator
//            function v = evalKKTError(z,i,rho,mu)
//
//              // Initialize optimality vector
//              kkt = zeros(input.nV+2*input.nE+2*input.nI,1);
//
//              // Set gradient of penalty objective
//              kkt(1:input.nV) = rho*this->g;
//
//              // Set gradient of Lagrangian for constraints
//              if input.nE > 0, kkt(1:input.nV) = kkt(1:input.nV) + (this->lE'*this->JE)'; end;
//              if input.nI > 0, kkt(1:input.nV) = kkt(1:input.nV) + (this->lI'*this->JI)'; end;
//
//              // Set complementarity for constraint slacks
//              if input.nE > 0, kkt(1+input.nV       :input.nV+2*input.nE       ) = [this->r1.*(1 + this->lE) - mu; this->r2.*(1 - this->lE) - mu]; end;
//              if input.nI > 0, kkt(1+input.nV+2*input.nE:input.nV+2*input.nE+2*input.nI) = [this->s1.*(0 + this->lI) - mu; this->s2.*(1 - this->lI) - mu]; end;
//
//              // Scale complementarity
//              if rho > 0, kkt = (1/max(1,norm(rho*this->g,inf)))*kkt; end;
//
//              // Evaluate optimality error
//              v = norm(kkt,inf);
//
//            end
//
//            // KKT errors evaluator
//            function evalKKTErrors(z,i)
//
//              // Loop to compute optimality errors
//              this->kkt(1) = this->evalKKTError(i,0    ,0   );
//              this->kkt(2) = this->evalKKTError(i,this->rho,0   );
//              this->kkt(3) = this->evalKKTError(i,this->rho,this->mu);
//
//            end
//
//            // Evaluator of lambda in original space
//            function l = evalLambdaOriginal(z,i)
//
//              // Initialize multipliers in original space
//              l = zeros(input.nE+input.n7+input.n8+input.n9,1);
//
//              // Scale equality constraint multipliers
//              if input.nE > 0, lE = this->lE.*(this->cEs/(this->rho*this->fs)); end;
//
//              // Set equality constraint multipliers in original space
//              if input.nE > 0, l(input.i6) = lE; end;
//
//              // Scale inequality constraint multipliers
//              if input.n7+input.n8+input.n9 > 0, lI = this->lInput.*(this->cIs/(this->rho*this->fs)); end;
//
//              // Set inequality constraint multipliers in original space
//              if input.n7 > 0, l(input.i7) = -lI(1+input.n3+input.n4+input.n5+input.n5               :input.n3+input.n4+input.n5+input.n5+input.n7               ); end;
//              if input.n8 > 0, l(input.i8) = +lI(1+input.n3+input.n4+input.n5+input.n5+input.n7          :input.n3+input.n4+input.n5+input.n5+input.n7+input.n8          ); end;
//              if input.n9 > 0, l(input.i9) = -lI(1+input.n3+input.n4+input.n5+input.n5+input.n7+input.n8     :input.n3+input.n4+input.n5+input.n5+input.n7+input.n8+input.n9     )...
//                                     +lI(1+input.n3+input.n4+input.n5+input.n5+input.n7+input.n8+input.n9:input.n3+input.n4+input.n5+input.n5+input.n7+input.n8+input.n9+input.n9); end;
//
//            end
//
//            // Matrices evaluator
//            function evalMatrices(z,p,i,c)
//
//              // Evaluate Hessian and Newton matrices
//              this->evalHessian     (  i,c);
//              this->evalNewtonMatrix(p,i,c);
//
//            end
//
//            // Merit evaluator
//            function evalMerit(z,i)
//
//              // Initialize merit for objective
//              this->phi = this->rho*this->f;
//
//              // Update merit for slacks
//              if input.nE > 0, this->phi = this->phi - this->mu*sum(log([this->r1;this->r2])) + sum([this->r1;this->r2]); end;
//              if input.nI > 0, this->phi = this->phi - this->mu*sum(log([this->s1;this->s2])) + sum(      this->s2 ); end;
//
//            end
//
//            // Newton matrix evaluator
//            function evalNewtonMatrix(z,p,i,c)
//
//              // Check for equality constraints
//              if input.nE > 0
//
//                // Set diagonal terms
//                for j = 1:input.nE, this->A(input.nV+     j,input.nV+     j) = (1+this->lE(j))/this->r1(j);
//                                this->A(input.nV+input.nE+j,input.nV+input.nE+j) = (1-this->lE(j))/this->r2(j); end;
//
//                // Set constraint Jacobian
//                this->A(1+input.nV+2*input.nE+2*input.nI:input.nV+3*input.nE+2*input.nI,1:input.nV) = this->JE;
//
//              end
//
//              // Check for inequality constraints
//              if input.nI > 0
//
//                // Set diagonal terms
//                for j = 1:input.nI, this->A(input.nV+2*input.nE+     j,input.nV+2*input.nE+     j) = (0+this->lI(j))/this->s1(j);
//                                this->A(input.nV+2*input.nE+input.nI+j,input.nV+2*input.nE+input.nI+j) = (1-this->lI(j))/this->s2(j); end;
//
//                // Set constraint Jacobian
//                this->A(1+input.nV+3*input.nE+2*input.nI:input.nV+3*input.nE+3*input.nI,1:input.nV) = this->JI;
//
//              end
//
//              // Set minimum potential shift
//              min_shift = max(p.shift_min,p.shift_factor1*this->shift);
//
//              // Initialize Hessian modification
//              if this->cut_ == 1, this->shift = min(p.shift_max,min_shift/p.shift_factor2); else this->shift = 0; end;
//
//              // Initialize inertia correction loop
//              done = 0; this->shift22 = 0;
//
//              // Loop until inertia is correct
//              while ~done && this->shift < p.shift_max
//
//                // Set Hessian of Lagrangian
//                this->A(1:input.nV,1:input.nV) = this->H+this->shift*speye(input.nV);
//
//                // Set diagonal terms
//                for j = 1:input.nE, this->A(input.nV+2*input.nE+2*input.nI+j,input.nV+2*input.nE+2*input.nI+j) = -this->shift22; end;
//
//                // Set diagonal terms
//                for j = 1:input.nI, this->A(input.nV+3*input.nE+2*input.nI+j,input.nV+3*input.nE+2*input.nI+j) = -this->shift22; end;
//
//                // Set number of nonzeros in (upper triangle of) Newton matrix
//                this->Annz = nnz(tril(this->A));
//
//                // Factor primal-dual matrix
//                [this->AL,this->AD,this->AP,this->AS,neig] = ldl(tril(this->A),p.pivot_thresh,'vector');
//
//                // Increment factorization counter
//                counter.incrementFactorizationCount;
//
//                // Set number of nonnegative eigenvalues
//                peig = input.nA - neig;
//
//                // Check inertia
//                if     peig < input.nV+2*input.nE+2*input.nI        , this->shift   = max(min_shift,this->shift/p.shift_factor2);
//                elseif neig < input.nE+input.nI & this->shift22 == 0, this->shift22 = p.shift_min;
//                else                                      done      = 1; end;
//
//              end
//
//              // Update Hessian
//              this->H = this->H+this->shift*speye(input.nV);
//
//            end
//
//            // Newton right-hand side evaluator
//            function evalNewtonRhs(z,i)
//
//              // Initialize right-hand side vector
//              this->b = zeros(input.nA,1);
//
//              // Set gradient of objective
//              this->b(1:input.nV) = this->rho*this->g;
//
//              // Set gradient of Lagrangian for constraints
//              if input.nE > 0, this->b(1:input.nV) = this->b(1:input.nV) + (this->lE'*this->JE)'; end;
//              if input.nI > 0, this->b(1:input.nV) = this->b(1:input.nV) + (this->lI'*this->JI)'; end;
//
//              // Set complementarity for constraint slacks
//              if input.nE > 0, this->b(1+input.nV       :input.nV+2*input.nE       ) = [1 + this->lE - this->mu./this->r1; 1 - this->lE - this->mu./this->r2]; end;
//              if input.nI > 0, this->b(1+input.nV+2*input.nE:input.nV+2*input.nE+2*input.nI) = [0 + this->lI - this->mu./this->s1; 1 - this->lI - this->mu./this->s2]; end;
//
//              // Set penalty-interior-point constraint values
//              if input.nE > 0, this->b(1+input.nV+2*input.nE+2*input.nI:input.nV+3*input.nE+2*input.nI) = this->cE + this->r1 - this->r2; end;
//              if input.nI > 0, this->b(1+input.nV+3*input.nE+2*input.nI:input.nV+3*input.nE+3*input.nI) = this->cI + this->s1 - this->s2; end;
//
//            end
//
//            // Scalings evaluator
//            function evalScalings(z,p,i,c)
//
//              // Initialize scalings
//              this->fs  = 1;
//              this->cEs = ones(input.nE,1);
//              this->cIs = ones(input.nI,1);
//
//              // Evaluate gradients
//              this->eval_gradients(i,c);
//
//              // Evaluate norm of objective gradient
//              g_norm_inf = norm(this->g,inf);
//
//              // Scale down objective if norm of gradient is too large
//              this->fs = p.grad_max/max(g_norm_inf,p.grad_max);
//
//              // Loop through equality constraints
//              for j = 1:input.nE
//
//                // Evaluate norm of gradient of jth equality constraint
//                JE_j_norm_inf = norm(this->JE(j,:),inf);
//
//                // Scale down equality constraint j if norm of gradient is too large
//                this->cEs(j) = p.grad_max/max(JE_j_norm_inf,p.grad_max);
//
//              end
//
//              // Loop through inequality constraints
//              for j = 1:input.nI
//
//                // Evaluate norm of gradient of jth inequality constraint
//                JI_j_norm_inf = norm(this->JI(j,:),inf);
//
//                // Scale down inequality constraint j if norm of gradient is too large
//                this->cIs(j) = p.grad_max/max(JI_j_norm_inf,p.grad_max);
//
//              end
//
//            end
//
//            // Slacks evaluator
//            function evalSlacks(z,p,i)
//
//              // Check for equality constraints
//              if input.nE > 0
//
//                // Set slacks
//                this->r1 = (1/2)*(this->mu - this->cE + sqrt(this->cE.^2 + this->mu^2));
//                this->r2 = (1/2)*(this->mu + this->cE + sqrt(this->cE.^2 + this->mu^2));
//
//                // Adjust for numerical error
//                this->r1 = max(this->r1,p.slack_min);
//                this->r2 = max(this->r2,p.slack_min);
//
//              end
//
//              // Check for inequality constraints
//              if input.nI > 0
//
//                // Set slacks
//                this->s1 = (1/2)*(2*this->mu - this->cI + sqrt(this->cInput.^2 + 4*this->mu^2));
//                this->s2 = (1/2)*(2*this->mu + this->cI + sqrt(this->cInput.^2 + 4*this->mu^2));
//
//                // Adjust for numerical error
//                this->s1 = max(this->s1,p.slack_min);
//                this->s2 = max(this->s2,p.slack_min);
//
//              end
//
//            end
//
//            // Evaluator of x in original space
//            function x = eval_x_original(z,i)
//
//              // Initialize x in original space
//              x = zeros(input.n0,1);
//
//              // Evaluate x in original space
//              x(input.i1) = this->x(1               :input.n1               );
//              x(input.i2) = input.b2;
//              x(input.i3) = this->x(1+input.n1          :input.n1+input.n3          );
//              x(input.i4) = this->x(1+input.n1+input.n3     :input.n1+input.n3+input.n4     );
//              x(input.i5) = this->x(1+input.n1+input.n3+input.n4:input.n1+input.n3+input.n4+input.n5);
//
//            end
//
//            // Gets primal-dual point
//            function sol = getSolution(z,i)
//              sol.x = this->eval_x_original(i);
//              sol.l = this->evalLambdaOriginal(i);
//            end
//
//            // Initializes Newton matrix
//            function initNewtonMatrix(z,i)
//
//              // Allocate memory
//              this->A = spalloc(input.nA,input.nA,this->Hnnz+5*input.nE+5*input.nI+this->JEnnz+this->JInnz);
//
//              // Initialize interior-point Hessians
//              this->A(1+input.nV       :input.nV+2*input.nE       ,1+input.nV       :input.nV+2*input.nE       ) = speye(2*input.nE);
//              this->A(1+input.nV+2*input.nE:input.nV+2*input.nE+2*input.nI,1+input.nV+2*input.nE:input.nV+2*input.nE+2*input.nI) = speye(2*input.nI);
//
//              // Check for equality constraints
//              if input.nE > 0
//
//                // Initialize constraint Jacobian
//                this->A(1+input.nV+2*input.nE+2*input.nI:input.nV+3*input.nE+2*input.nI,1+input.nV     :input.nV+  input.nE) =  speye(input.nE);
//                this->A(1+input.nV+2*input.nE+2*input.nI:input.nV+3*input.nE+2*input.nI,1+input.nV+input.nE:input.nV+2*input.nE) = -speye(input.nE);
//
//              end
//
//              // Check for inequality constraints
//              if input.nI > 0
//
//                // Initialize constraint Jacobian
//                this->A(1+input.nV+3*input.nE+2*input.nI:input.nV+3*input.nE+3*input.nI,1+input.nV+2*input.nE     :input.nV+2*input.nE+  input.nI) =  speye(input.nI);
//                this->A(1+input.nV+3*input.nE+2*input.nI:input.nV+3*input.nE+3*input.nI,1+input.nV+2*input.nE+input.nI:input.nV+2*input.nE+2*input.nI) = -speye(input.nI);
//
//              end
//
//            end
//
//            // Set interior-point parameter
//            function setMu(z,mu)
//              this->mu = mu;
//            end
//
//            // Set primal variables
//            function setPrimals(z,i,x,r1,r2,s1,s2,lE,lI,f,cE,cI,phi)
//
//              // Set primal variables
//              this->x = x; this->f = f;
//              if input.nE > 0, this->cE = cE; this->r1 = r1; this->r2 = r2; this->lE = lE; end;
//              if input.nI > 0, this->cI = cI; this->s1 = s1; this->s2 = s2; this->lI = lI; end;
//              this->phi = phi;
//
//            end
//
//            // Set penalty parameter
//            function setRho(z,rho)
//              this->rho = rho;
//            end
//
//            // Set last penalty parameter
//            function setRhoLast(z,rho)
//              this->rho_ = rho;
//            end
//
//            // Iterate updater
//            function updateIterate(z,p,i,c,d,a)
//
//              // Update last quantities
//              this->v_   = this->v;
//              this->cut_ = (a.p < a.p0);
//
//              // Update iterate quantities
//              this->updatePoint      (  i,  d,a);
//              this->evalInfeasibility(  i      );
//              this->eval_gradients    (  i,c    );
//              this->evalDependent    (p,i      );
//
//              // Update last KKT errors
//              this->kkt_ = [this->kkt(2); this->kkt_(1:p.opt_err_mem-1)];
//
//            end
//
//            // Parameter updater
//            function updateParameters(z,p,i)
//
//              // Check for interior-point parameter update based on optimality error
//              while this->mu > p.mu_min && this->kkt(3) <= max([this->mu;p.opt_err_tol-this->mu])
//
//                // Restrict interior-point parameter increase
//                p.setMuMaxExpZero;
//
//                // Update interior-point parameter
//                if this->mu > p.mu_min
//
//                  // Decrease interior-point
//                  this->mu = max(p.mu_min,min(p.mu_factor*this->mu,this->mu^p.mu_factor_exp));
//
//                  // Evaluate penalty and interior-point parameter dependent quantities
//                  this->evalDependent(p,i);
//
//                end
//
//              end
//
//              // Check for penalty parameter update based on optimality error
//              if (this->kkt(2) <= p.opt_err_tol && this->v > p.opt_err_tol) || this->v > max([1;this->v_;p.infeas_max])
//
//                // Update penalty parameter
//                if this->rho > p.rho_min
//
//                  // Decrease penalty parameter
//                  this->rho = max(p.rho_min,p.rho_factor*this->rho);
//
//                  // Evaluate penalty and interior-point parameter dependent quantities
//                  this->evalDependent(p,i);
//
//                end
//
//              end
//
//            end
//
//            // Primal point updater
//            function updatePoint(z,i,d,a)
//
//              // Update primal and dual variables
//                           this->x  = this->x  + a.p*d.x ;
//              if input.nE > 0, this->r1 = this->r1 + a.p*d.r1;
//                           this->r2 = this->r2 + a.p*d.r2; end;
//              if input.nI > 0, this->s1 = this->s1 + a.p*d.s1;
//                           this->s2 = this->s2 + a.p*d.s2; end;
//              if input.nE > 0, this->lE = this->lE + a.d*d.lE; end;
//              if input.nI > 0, this->lI = this->lI + a.d*d.lI; end;
//
//            end
//
//          end
//
//          // Class methods (static)
//          methods (Static)
//
//            // Feasibility violation evaluator
//            function v = evalViolation(i,cE,cI)
//
//              // Initialize violation vector
//              vec = [];
//
//              // Update vector for constraint values
//              if input.nE > 0, vec = cE;               end;
//              if input.nI > 0, vec = [vec; max(cI,0)]; end;
//
//              // Evaluate vector norm
//              v = norm(vec,1);
//
//            end
//
//          end
//
//        end

  }; // struct Iterate

} // namespace Pipal

#endif /* INCLUDE_PIPAL_ITERATE_HH */
