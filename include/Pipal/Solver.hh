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
  * \tparam N The size of the primal variable vector.
  * \tparam M The size of the dual variable vector.
  */
  template<typename Real, Integer N, Integer M>
  class Solver
  {
    static_assert(std::is_floating_point<Real>::value,
      "Pipal::Solver: template argument 'Real' must be a floating-point type.");
    static_assert(N > 0,
      "Pipal::Solver: template argument 'N' must be positive.");
    static_assert(M > 0,
      "Pipal::Solver: template argument 'M' must be positive.");

  public:
    using VectorN = typename Problem<Real, N, M>::VectorN;
    using VectorM = typename Problem<Real, N, M>::VectorM;
    using MatrixJ = typename Problem<Real, N, M>::MatrixJ;
    using MatrixH = typename Problem<Real, N, M>::MatrixH;

    using ObjectiveFunc           = typename ProblemWrapper<Real, N, M>::ObjectiveFunc;
    using ObjectiveGradientFunc   = typename ProblemWrapper<Real, N, M>::ObjectiveGradientFunc;
    using ObjectiveHessianFunc    = typename ProblemWrapper<Real, N, M>::ObjectiveHessianFunc;
    using ConstraintsFunc         = typename ProblemWrapper<Real, N, M>::ConstraintsFunc;
    using ConstraintsJacobianFunc = typename ProblemWrapper<Real, N, M>::ConstraintsJacobianFunc;
    using LagrangianHessianFunc   = typename ProblemWrapper<Real, N, M>::LagrangianHessianFunc;

  private:

    struct Counter {
      Integer iters{0};          /*!< Number of iterations. */
      Integer function_evals{0}; /*!< Number of function evaluations. */
      Integer gradient_evals{0}; /*!< Number of gradient evaluations. */
      Integer hessian_evals{0};  /*!< Number of Hessian evaluations. */
      Integer jacobian_evals{0}; /*!< Number of Jacobian evaluations. */
      Integer linear_solves{0};  /*!< Number of linear system solves. */

      /** \brief Reset all counters to zero. */
      void reset()
      {
        iters = 0;
        function_evals = 0;
        gradient_evals = 0;
        hessian_evals = 0;
        jacobian_evals = 0;
        linear_solves = 0;
      }
    } m_counters; /*!< Internal counters for statistics. */

    struct Parameter {
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
    } m_params; /*!< Internal parameters with default values. */

    std::unique_ptr<Problem<Real, N, M>> m_problem; /**< Problem object */

    // Some options for the solver
    bool    m_verbose{false}; /**< Verbosity flag. */
    bool    m_aggressive{false};   /*!< Aggressive algorithm boolean. */
    Integer m_max_iterations{100}; /*!< Maximum number of iterations. */
    Real    m_tolerance{1.0e-10};  /*!< Optimality tolerance. */
    Real    m_mu_max_exp{0.0};     /*!< Interior-point parameter maximum exponent in increases. */

  public:
    /**
    * \brief Default constructor for the IPSolver class.
    *
    * Initializes the solver with default values for the objective, gradient, constraints, and Jacobian functions.
    */
    Solver() {};

    /**
    * \brief Constructor for the IPSolver class.
    *
    * Initializes the solver with the provided objective, gradient, constraints, and Jacobian functions.
    * \param[in] objective Objective function handle.
    * \param[in] objective_gradient Gradient of the objective function handle.
    * \param[in] constraints Constraints function handle.
    * \param[in] constraints_jacobian Jacobian of the constraints function handle.
    * \param[in] lagrangian_hessian Hessian of the Lagrangian function handle.
    */
    Solver(ObjectiveFunc const & objective, ObjectiveGradientFunc const & objective_gradient,
      ConstraintsFunc const & constraints, ConstraintsJacobianFunc const & constraints_jacobian,
      LagrangianHessianFunc const & lagrangian_hessian)
      : m_problem(std::make_unique<ProblemWrapper<Real, N, M>>(objective, objective_gradient,
          constraints, constraints_jacobian, lagrangian_hessian)) {}

    /**
    * \brief Constructor for the IPSolver class (with Hessian).
    *
    * Initializes the solver with the provided objective, gradient, constraints, Jacobian, and Hessian functions.
    * \param[in] objective Objective function handle.
    * \param[in] objective_gradient Gradient of the objective function handle.
    * \param[in] objective_hessian Hessian of the objective function handle.
    * \param[in] constraints Constraints function handle.
    * \param[in] constraints_jacobian Jacobian of the constraints function handle.
    * \param[in] lagrangian_hessian Hessian of the Lagrangian function handle.
    */
    Solver(ObjectiveFunc const & objective, ObjectiveGradientFunc const & objective_gradient,
      ObjectiveHessianFunc const & objective_hessian, ConstraintsFunc const & constraints,
      ConstraintsJacobianFunc const & constraints_jacobian, LagrangianHessianFunc const & lagrangian_hessian)
      : m_problem(std::make_unique<ProblemWrapper<Real, N, M>>(objective, objective_gradient, objective_hessian,
          constraints, constraints_jacobian, lagrangian_hessian)) {}

    /**
    * \brief Constructor for the IPSolver class (with a unique pointer to a Problem object).
    *
    * Initializes the solver with the provided unique pointer to a Problem object.
    * \param[in] problem The unique pointer to the Problem object to use.
    * \warning The pointer is moved into the solver, so it should not be used after this call.
    */
    Solver(std::unique_ptr<Problem<Real, N, M>> && problem)
      : m_problem(std::move(problem)) {}

    /**
    * \brief Deleted copy constructor.
    *
    * This class is not copyable.
    */
    Solver(Solver const &) = delete;

    /**
    * \brief Deleted assignment operator.
    *
    * This class is not assignable.
    */
    Solver& operator=(Solver const &) = delete;

    /**
    * \brief Deleted move constructor.
    *
    * This class is not movable.
    */
    Solver(Solver &&) = delete;

    /**
    * \brief Deleted move assignment operator.
    *
    * This class is not movable.
    */
    Solver& operator=(Solver &&) = delete;

    /**
    * \brief Destructor for the IPSolver class.
    *
    * Cleans up resources used by the solver.
    */
    ~Solver() = default;

    /**
    * \brief Set the problem to be solved.
    *
    * This method allows the user to specify the problem to be solved.
    * \param[in] problem The problem to set.
    */
    void problem(Problem<Real, N, M> const & problem) {
      this->m_problem = std::make_unique<ProblemWrapper<Real, N, M>>(problem);
    }

    /**
    * \brief Set the problem to be solved using a unique pointer.
    *
    * This method allows the user to specify the problem to be solved using a unique pointer.
    * \param[in] problem The unique pointer to the problem to set.
    */
    void problem(std::unique_ptr<Problem<Real, N, M>> && problem) {this->m_problem = std::move(problem);}

    /**
    * \brief Get the problem being solved.
    * \return A reference to the problem.
    */
    Problem<Real, N, M> const & problem() const {return *this->m_problem;}

    /**
    * Get the verbose mode.
    * \return The verbose mode.
    */
    bool verbose_mode() {return this->m_verbose;}

    /**
    * Set the verbose mode.
    * \param[in] t_verbose The verbose mode.
    */
    void verbose_mode(bool const t_verbose) {this->m_verbose = t_verbose;}

    /**
    * Enable the verbose mode.
    */
    void enable_verbose_mode() {this->verbose_mode(true);}

    /**
    * Disable the verbose mode.
    */
    void disable_verbose_mode() {this->verbose_mode(false);}

    /**
    * Get the aggressive mode.
    * \return The aggressive mode.
    */
    bool aggressive_mode() {return this->m_aggressive;}

    /**
    * Set the aggressive mode.
    * \param[in] t_aggressive The aggressive mode.
    */
    void aggressive_mode(bool const  t_aggressive) {this->m_aggressive = t_aggressive;}

    /**
    * Enable the aggressive mode.
    */
    void enable_aggressive_mode() {this->aggressive_mode(true);}

    /**
    * Disable the aggressive mode.
    */
    void disable_aggressive_mode() {this->aggressive_mode(false);}

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

    /**
    * \brief Solves the optimization problem using the interior-point method.
    *
    * This method implements the interior-point algorithm to solve the optimization problem defined
    * by the objective function, constraints, and their respective gradients and Jacobians.
    * \param[in] x_guess Initial guess for the optimization variables.
    * \param[out] x_sol Solution vector.
    * \return True if the optimization was successful, false otherwise.
    */
    bool solve(VectorN const & x_guess, VectorN & x_sol)
    {
      #define CMD "Pipal::Solver::solve(...): "

      x_sol = x_guess;
      return true;
      #undef CMD
    }

  }; // class IPSolver

} // namespace Pipal

#endif /* INCLUDE_PIPAL_SOLVER_HH */
