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

#ifndef INCLUDE_PIPAL_SOLVER_HH
#define INCLUDE_PIPAL_SOLVER_HH

// Pipal includes
#include "Pipal/Problem.hh"

// Standard libraries
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <type_traits>

// Eigen library
#include <Eigen/Dense>

namespace Pipal {

  /**
   * \brief Solver class for the Pipal library.
   *
   * The Solver class provides an interface for solving the optimization problems using Frank E.
   * Curtis Pipal algorithm. It utilizes the Problem class to define the optimization problem and
   * implements various methods for solving it.
   * \tparam Real The floating-point type.
   */
  template<typename Real>
  class Solver
  {
    static_assert(std::is_floating_point<Real>::value,
      "Pipal::Solver: template argument 'Real' must be a floating-point type.");

  public:
    using Vector = typename Problem<Real>::Vector;
    using Matrix = typename Problem<Real>::Matrix;
    using Mask   = Eigen::Array<bool, Eigen::Dynamic, 1>;

    using ObjectiveFunc           = typename ProblemWrapper<Real>::ObjectiveFunc;
    using ObjectiveGradientFunc   = typename ProblemWrapper<Real>::ObjectiveGradientFunc;
    using ObjectiveHessianFunc    = typename ProblemWrapper<Real>::ObjectiveHessianFunc;
    using ConstraintsFunc         = typename ProblemWrapper<Real>::ConstraintsFunc;
    using ConstraintsJacobianFunc = typename ProblemWrapper<Real>::ConstraintsJacobianFunc;
    using LagrangianHessianFunc   = typename ProblemWrapper<Real>::LagrangianHessianFunc;
    using BoundsFunc              = typename ProblemWrapper<Real>::BoundsFunc;

    using Algorithm = enum class Algorithm : Integer {CONSERVATIVE = 0, ADAPTIVE = 1}; /*!< Algorithm choice. */

  private:

    /**
     * \brief Internal counters for solver statistics.
     */
    struct Counter
    {
      Integer iters{0};          /*!< Number of iterations. */
      Integer function_evals{0}; /*!< Number of function evaluations. */
      Integer gradient_evals{0}; /*!< Number of gradient evaluations. */
      Integer hessian_evals{0};  /*!< Number of Hessian evaluations. */
      Integer jacobian_evals{0}; /*!< Number of Jacobian evaluations. */
      Integer linear_solves{0};  /*!< Number of linear system solves. */
    } m_counters;

    /**
     * \brief Internal parameters for the solver algorithm.
     */
    struct Parameter
    {
      const Real    rhs_bnd{1.0e+18};       /*!< Maximum absolute value allowed for constraint right-hand side. */
      const Real    grad_max{1.0e+02};      /*!< Gradient norm limit for scaling. */
      const Real    infeas_max{1.0e+02};    /*!< Infeasibility limit for penalty parameter update. */
      const Real    nnz_max{2.0e+04};       /*!< Maximum non-zeros in (upper triangle of) Newton matrix. */
      const Integer opt_err_mem{6};         /*!< Optimality error history length. */
      const Real    ls_factor{5.0e-01};     /*!< Line search reduction factor. */
      const Real    ls_thresh{1.0e-08};     /*!< Line search threshold value. */
      const Real    ls_frac{1.0e-02};       /*!< Line search fraction-to-boundary constant. */
      const Real    pivot_thresh{5.0e-01};  /*!< Pivot threshold for LDL factorization. */
      const Real    slack_min{1.0e-20};     /*!< Slack variable bound. */
      const Real    shift_min{1.0e-12};     /*!< Hessian shift (non-zero) minimum value. */
      const Real    shift_factor1{5.0e-01}; /*!< Hessian shift update value (for decreases). */
      const Real    shift_factor2{6.0e-01}; /*!< Hessian shift update value (for increases). */
      const Real    shift_max{1.0e+08};     /*!< Hessian shift maximum value. */
      const Real    rho_init{1.0e-01};      /*!< Penalty parameter initial value. */
      const Real    rho_min{1.0e-12};       /*!< Penalty parameter minimum value. */
      const Real    rho_factor{5.0e-01};    /*!< Penalty parameter reduction factor. */
      const Integer rho_trials{8};          /*!< Penalty parameter number of trial values per iteration. */
      const Real    mu_init{1.0e-01};       /*!< Interior-point parameter initial value. */
      const Real    mu_min{1.0e-12};        /*!< Interior-point parameter minimum value. */
      const Real    mu_factor{1.0e-01};     /*!< Interior-point parameter reduction factor. */
      const Real    mu_factor_exp{1.5};     /*!< Interior-point parameter reduction exponent. */
      const Integer mu_trials{4};           /*!< Interior-point parameter number of trial values per iteration. */
      const Real    mu_max{1.0e-01};        /*!< Interior-point parameter maximum value. */
      const Real    mu_max_exp_0{0.0};      /*!< Interior-point parameter maximum exponent in increases (default). */
      const Real    update_con_1{1.0e-02};  /*!< Steering rule constant 1. */
      const Real    update_con_2{1.0e-02};  /*!< Steering rule constant 2. */
      const Real    update_con_3{1.01};     /*!< Adaptive interior-point rule constant. */
    } m_params;

    struct Input
    {
      Mask    i1; /*!< Mask of free variables. */
      Mask    i2; /*!< Mask of fixed variables. */
      Mask    i3; /*!< Mask of lower bounded variables. */
      Mask    i4; /*!< Mask of upper bounded variables. */
      Mask    i5; /*!< Mask of lower and upper bounded variables. */
      Mask    i6; /*!< Mask of equality constraints. */
      Mask    i7; /*!< Mask of lower bounded constraints. */
      Mask    i8; /*!< Mask of upper bounded constraints. */
      Mask    i9; /*!< Mask of lower and upper bounded constraints. */
      Vector  x0; /*!< Initial guess for the primal variables. */
      Vector  b2; /*!< Right-hand side of fixed variables. */
      Vector  l3; /*!< Right-hand side of lower bounded variables. */
      Vector  u4; /*!< Right-hand side of upper bounded variables. */
      Vector  l5; /*!< Right-hand side of lower half of lower and upper bounded variables. */
      Vector  u5; /*!< Right-hand side of upper half of lower and upper bounded variables. */
      Vector  b6; /*!< Right-hand side of equality constraints. */
      Vector  l7; /*!< Right-hand side of lower bounded constraints. */
      Vector  u8; /*!< Right-hand side of upper bounded constraints. */
      Vector  l9; /*!< Right-hand side of lower half of lower and upper bounded constraints. */
      Vector  u9; /*!< Right-hand side of upper half of lower and upper bounded constraints. */
      Integer n0; /*!< Number of original formulation variables. */
      Integer n1; /*!< Number of free variables. */
      Integer n2; /*!< Number of fixed variables. */
      Integer n3; /*!< Number of lower bounded variables. */
      Integer n4; /*!< Number of upper bounded variables. */
      Integer n5; /*!< Number of lower and upper bounded variables. */
      Integer n6; /*!< Number of equality constraints. */
      Integer n7; /*!< Number of lower bounded constraints. */
      Integer n8; /*!< Number of upper bounded constraints. */
      Integer n9; /*!< Number of lower and upper bounded constraints. */
      Integer nV; /*!< Number of variables. */
      Integer nI; /*!< Number of inequality constraints. */
      Integer nE; /*!< Number of equality constraints. */
      Integer nA; /*!< Size of primal-dual matrix. */
      Integer vi; /*!< Counter for invalid bounds. */
    } m_input;

    std::unique_ptr<Problem<Real>> m_problem; /**< Problem object */

    // Some options for the solver
    Algorithm m_algorithm{Algorithm::ADAPTIVE}; /*!< Algorithm choice. */
    bool      m_verbose{false};      /**< Verbosity flag. */
    Integer   m_max_iterations{100}; /*!< Maximum number of iterations. */
    Real      m_tolerance{1.0e-10};  /*!< Optimality tolerance. */
    Real      m_mu_max_exp{0.0};     /*!< Interior-point parameter maximum exponent in increases. */

  public:
    /**
     * \brief Default constructor for the Pipal class.
     *
     * Initializes the solver with default values for the objective, gradient, constraints, and Jacobian functions.
     */
    Solver() {};

    /**
     * \brief Constructor for the Pipal class.
     *
     * Initializes the solver with the provided objective, gradient, constraints, and Jacobian functions.
     * \param[in] objective Objective function handle.
     * \param[in] objective_gradient Gradient of the objective function handle.
     * \param[in] constraints Constraints function handle.
     * \param[in] constraints_jacobian Jacobian of the constraints function handle.
     * \param[in] lagrangian_hessian Hessian of the Lagrangian function handle.
     * \param[in] primal_lower_bounds Lower bounds on the primal variables handle.
     * \param[in] primal_upper_bounds Upper bounds on the primal variables handle.
     * \param[in] constraints_lower_bounds Lower bounds on the constraints handle.
     * \param[in] constraints_upper_bounds Upper bounds on the constraints handle.
     */
    Solver(ObjectiveFunc const & objective, ObjectiveGradientFunc const & objective_gradient,
      ConstraintsFunc const & constraints, ConstraintsJacobianFunc const & constraints_jacobian,
      LagrangianHessianFunc const & lagrangian_hessian, BoundsFunc const & primal_lower_bounds,
      BoundsFunc const & primal_upper_bounds, BoundsFunc const & constraints_lower_bounds,
      BoundsFunc const & constraints_upper_bounds)
      : m_problem(std::make_unique<ProblemWrapper<Real>>(objective, objective_gradient,
        constraints, constraints_jacobian, lagrangian_hessian, primal_lower_bounds, primal_upper_bounds,
        constraints_lower_bounds, constraints_upper_bounds)) {}

    /**
     * \brief Constructor for the Pipal class (with Hessian).
     *
     * Initializes the solver with the provided objective, gradient, constraints, Jacobian, and Hessian functions.
     * \param[in] objective Objective function handle.
     * \param[in] objective_gradient Gradient of the objective function handle.
     * \param[in] objective_hessian Hessian of the objective function handle.
     * \param[in] constraints Constraints function handle.
     * \param[in] constraints_jacobian Jacobian of the constraints function handle.
     * \param[in] lagrangian_hessian Hessian of the Lagrangian function handle.
     * \param[in] primal_lower_bounds Lower bounds on the primal variables handle.
     * \param[in] primal_upper_bounds Upper bounds on the primal variables handle.
     * \param[in] constraints_lower_bounds Lower bounds on the constraints handle.
     * \param[in] constraints_upper_bounds Upper bounds on the constraints handle.
     */
    Solver(ObjectiveFunc const & objective, ObjectiveGradientFunc const & objective_gradient,
      ObjectiveHessianFunc const & objective_hessian, ConstraintsFunc const & constraints,
      ConstraintsJacobianFunc const & constraints_jacobian, LagrangianHessianFunc const & lagrangian_hessian,
      BoundsFunc const & primal_lower_bounds, BoundsFunc const & primal_upper_bounds,
      BoundsFunc const & constraints_lower_bounds, BoundsFunc const & constraints_upper_bounds)
      : m_problem(std::make_unique<ProblemWrapper<Real>>(objective, objective_gradient, objective_hessian,
        constraints, constraints_jacobian, lagrangian_hessian, primal_lower_bounds, primal_upper_bounds,
        constraints_lower_bounds, constraints_upper_bounds)) {}

    /**
     * \brief Constructor for the Pipal class (with a unique pointer to a Problem object).
     *
     * Initializes the solver with the provided unique pointer to a Problem object.
     * \param[in] problem The unique pointer to the Problem object to use.
     * \warning The pointer is moved into the solver, so it should not be used after this call.
     */
    Solver(std::unique_ptr<Problem<Real>> && problem)
      : m_problem(std::move(problem)) {}

    /**
     * \brief Deleted copy constructor.
     * \note This class is not copyable.
     */
    Solver(Solver const &) = delete;

    /**
     * \brief Deleted assignment operator.
     * \note This class is not assignable.
     */
    Solver & operator=(Solver const &) = delete;

    /**
     * \brief Deleted move constructor.
     * \note This class is not movable.
     */
    Solver(Solver &&) = delete;

    /**
     * \brief Deleted move assignment operator.
     * \note This class is not movable.
     */
    Solver & operator=(Solver &&) = delete;

    /**
     * \brief Destructor for the Pipal class.
     *
     * Cleans up resources used by the solver.
     */
    ~Solver() = default;

    /**
     * \brief Get the counters for the solver.
     * \return A constant reference to the counters structure.
     */
    Counter const & counters() const {return this->m_counters;}

    /**
     * \brief Set the problem to be solved.
     *
     * This method allows the user to specify the problem to be solved.
     * \param[in] problem The problem to set.
     */
    void problem(Problem<Real> const & problem) {
      this->m_problem = std::make_unique<ProblemWrapper<Real>>(problem);
    }

    /**
     * \brief Set the problem to be solved using a unique pointer.
     *
     * This method allows the user to specify the problem to be solved using a unique pointer.
     * \param[in] problem The unique pointer to the problem to set.
     */
    void problem(std::unique_ptr<Problem<Real>> && problem) {this->m_problem = std::move(problem);}

    /**
     * \brief Get the problem being solved.
     * \return A reference to the problem.
     */
    Problem<Real> const & problem() const {return *this->m_problem;}

    /**
     * \brief Get the verbose mode.
     * \return The verbose mode.
     */
    bool verbose_mode() {return this->m_verbose;}

    /**
     * \brief Set the verbose mode.
     * \param[in] t_verbose The verbose mode.
     */
    void verbose_mode(bool const t_verbose) {this->m_verbose = t_verbose;}

    /**
     * \brief Get the aggressive mode.
     * \return The aggressive mode.
     */
    bool aggressive_mode() {return this->m_aggressive;}

    /**
     * \brief Set the adaptive mode.
     */
    void algorithm(Algorithm const t_algorithm) {this->m_algorithm = t_algorithm;}

    /**
     * \brief Set the convergence tolerance for the solver.
     *
     * This method allows the user to specify the tolerance for convergence.
     * The solver will stop when the residuals are below this tolerance.
     *
     * \param[in] t_tolerance The convergence tolerance.
     */
    void tolerance(Real const t_tolerance)
    {
      PIPAL_ASSERT(t_tolerance > static_cast<Real>(0.0),
        "Pipal::Solver::tolerance(...): input value must be positive");
      this->m_tolerance = t_tolerance;
    }

    /**
     * \brief Get the tolerance for convergence.
     * \return The tolerance for convergence.
     */
    Real tolerance() const {return this->m_tolerance;}

    /**
     * \brief Set the maximum number of iterations for the solver.
     *
     * This method allows the user to specify the maximum number of iterations
     * the solver will perform before stopping.
     *
     * \param[in] t_max_iterations The maximum number of iterations.
     */
    void max_iterations(Integer const t_max_iterations)
    {
      PIPAL_ASSERT(t_max_iterations > 0,
        "Pipal::Solver::max_iterations(...): input value must be positive");
      this->m_max_iterations = t_max_iterations;
    }

    /**
     * \brief Get the maximum number of iterations.
     * \return The maximum number of iterations.
     */
    Integer max_iterations() const {return this->m_max_iterations;}

    /**
     * \brief Solves the optimization problem using the interior-point method.
     *
     * This method implements the interior-point algorithm to solve the optimization problem defined
     * by the objective function, constraints, and their respective gradients and Jacobians.
     * \param[in] x_guess Initial guess for the optimization variables.
     * \param[out] x_sol Solution vector.
     * \return True if the optimization was successful, false otherwise.
     */
    bool optimize(Vector const & x_guess, Vector & x_sol)
    {
      #define CMD "Pipal::Solver::optimize(...): "

      // Reset counters
      this->reset_counters();

      // Check that the problem is set
      PIPAL_ASSERT(this->m_problem != nullptr,
        CMD "problem not set, use 'problem(...)' method to set it");

      // Fill input structure
      this->fill_input(x_guess);

      x_sol = x_guess;
      return true;

      #undef CMD
    }

  private:
    /**
     * \brief Reset all internal counters to zero.
     */
    void reset_counters()
    {
      this->m_counters.iters          = 0;
      this->m_counters.function_evals = 0;
      this->m_counters.gradient_evals = 0;
      this->m_counters.hessian_evals  = 0;
      this->m_counters.jacobian_evals = 0;
      this->m_counters.linear_solves  = 0;
    }

    /**
     * \brief Select elements from a vector based on a boolean mask.
     * \param[in] vector The input vector.
     * \param[in] mask The boolean mask.
     * \return The selected elements from the input vector.
     */
    static Vector select_elements(Vector const & vector, Mask const & mask)
    {
      #define CMD "Pipal::Solver::select_elements(...): "

      Integer const size{static_cast<Integer>(vector.size())};
      Integer const count{static_cast<Integer>(mask.count())};

      PIPAL_ASSERT(size == mask.size(),
        CMD "vector and mask must have the same size.");

      if (count == size) {
        return vector;
      } else {
        Vector out(count);
        for (Integer i{0}, j{0}; i < size; ++i) {
          if (mask[i]) {out[j++] = vector[i];}
        }
        return out;
      }

      #undef CMD
    }

    /**
     * \brief Fill the input structure with problem data.
     * \param[in] x0 Initial guess for the optimization variables.
     */
    void fill_input(Vector const & x0)
    {
      #define CMD "Pipal::Solver::fill_input(...): "

      // Create alias for easier access
      Parameter const & params{this->m_params};
      Input & input{this->m_input};

      // Get variable bounds
      Vector bl, bu;
      PIPAL_ASSERT(this->m_problem->primal_lower_bounds(bl),
        CMD "error in evaluating lower bounds on primal variables");
      PIPAL_ASSERT(this->m_problem->primal_upper_bounds(bu),
        CMD "error in evaluating upper bounds on primal variables");

      // Get constraint bounds
      Vector cl, cu;
      PIPAL_ASSERT(this->m_problem->constraints_lower_bounds(cl),
        CMD "error in evaluating lower bounds on constraints");
      PIPAL_ASSERT(this->m_problem->constraints_upper_bounds(cu),
        CMD "error in evaluating upper bounds on constraints");

      // Find index sets
      Real const tolerance{Eigen::NumTraits<Real>::epsilon()};
      Mask const cond_bl(bl.array() <= -params.rhs_bnd);
      Mask const cond_bu(bu.array() >=  params.rhs_bnd);
      Mask const cond_cl(cl.array() <= -params.rhs_bnd);
      Mask const cond_cu(cu.array() >=  params.rhs_bnd);

      input.i1 =  cond_bl && cond_bu;
      input.i2 = (bl.array() - bu.array()).abs() > tolerance;
      input.i3 = !cond_bl &&  cond_bu;
      input.i4 =  cond_bl && !cond_bu;
      input.i5 = !cond_bl && !cond_bu && !input.i2.array();
      input.i6 = (cl.array() - cu.array()).abs() > tolerance;
      input.i7 = !cond_cl &&  cond_cu;
      input.i8 =  cond_cl && !cond_cu;
      input.i9 = !cond_cl && !cond_cu && !input.i6.array();

      // Set right-hand side values
      input.b2 = this->select_elements(bl, input.i2);
      input.l3 = this->select_elements(bl, input.i3);
      input.u4 = this->select_elements(bu, input.i4);
      input.l5 = this->select_elements(bl, input.i5);
      input.u5 = this->select_elements(bu, input.i5);
      input.b6 = this->select_elements(cl, input.i6);
      input.l7 = this->select_elements(cl, input.i7);
      input.u8 = this->select_elements(cu, input.i8);
      input.l9 = this->select_elements(cl, input.i9);
      input.u9 = this->select_elements(cu, input.i9);

      // Set sizes of index sets
      input.n0 = x0.size();
      input.n1 = input.i1.count();
      input.n2 = input.i2.count();
      input.n3 = input.i3.count();
      input.n4 = input.i4.count();
      input.n5 = input.i5.count();
      input.n6 = input.i6.count();
      input.n7 = input.i7.count();
      input.n8 = input.i8.count();
      input.n9 = input.i9.count();

      // Initialize number of invalid bounds
      input.vi = 0;

      // Count invalid bounds
      if (input.n2 > 0) {
        input.vi += (input.b2.array() <= -params.rhs_bnd).count();
        input.vi += (input.b2.array() >=  params.rhs_bnd).count();
      }
      if (input.n3 > 0) {
        input.vi += (input.l3.array() >=  params.rhs_bnd).count();
      }
      if (input.n4 > 0) {
        input.vi += (input.u4.array() <= -params.rhs_bnd).count();
      }
      if (input.n5 > 0) {
        input.vi += (input.l5.array() >=  params.rhs_bnd).count();
        input.vi += (input.u5.array() <= -params.rhs_bnd).count();
        input.vi += (input.l5.array() >   input.u5.array()).count();
      }
      if (input.n6 > 0) {
        input.vi += (input.b6.array() <= -params.rhs_bnd).count();
        input.vi += (input.b6.array() >=  params.rhs_bnd).count();
      }
      if (input.n7 > 0) {
        input.vi += (input.l7.array() >=  params.rhs_bnd).count();
      }
      if (input.n8 > 0) {
        input.vi += (input.u8.array() <= -params.rhs_bnd).count();
      }
      if (input.n9 > 0) {
        input.vi += (input.l9.array() >=  params.rhs_bnd).count();
        input.vi += (input.u9.array() <= -params.rhs_bnd).count();
        input.vi += (input.l9.array() >   input.u9.array()).count();
      }

      // Set number of variables and constraints
      input.nV = input.n1 + input.n3 + input.n4 + input.n5;
      input.nI = input.n3 + input.n4 + 2*input.n5 + input.n7 + input.n8 + 2*input.n9;
      input.nE = input.n6;

      // Set size of primal-dual matrix
      input.nA = input.nV + 3*input.nE + 3*input.nI;

      // Set initial point
      input.x0.resize(input.nV);

      input.x0 <<
        this->select_elements(x0, input.i1), this->select_elements(x0, input.i3),
        this->select_elements(x0, input.i4), this->select_elements(x0, input.i5);
      std::cout << "input.x0 = " << input.x0.transpose() << std::endl;

      #undef CMD
    }

    /**
     * \brief Set the interior-point parameter maximum exponent in increases.
     * \param[in] t_mu_max_exp The interior-point parameter maximum exponent in increases (must be positive).
     */
    void mu_max_exp(Real const t_mu_max_exp)
    {
      PIPAL_ASSERT(t_mu_max_exp > static_cast<Real>(0.0),
        "Pipal::Solver::mu_max_exp(...): input value must be positive");
      this->m_mu_max_exp = t_mu_max_exp;
    }

    /**
     * \brief Get the interior-point parameter maximum exponent in increases.
     * \return The interior-point parameter maximum exponent in increases value.
     */
    Real mu_max_exp() const {return this->m_mu_max_exp;}


  }; // class Pipal

} // namespace Pipal

#endif /* INCLUDE_PIPAL_SOLVER_HH */
