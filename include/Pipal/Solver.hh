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

// Standard libraries
#include <iomanip>
#include <algorithm>
#include <type_traits>

// Eigen library
#include <Eigen/Dense>

// Pipal includes
#include "Pipal/Acceptance.hh"
#include "Pipal/Counter.hh"
#include "Pipal/Types.hh"
#include "Pipal/Direction.hh"
#include "Pipal/Input.hh"
#include "Pipal/Iterate.hh"
#include "Pipal/Output.hh"
#include "Pipal/Problem.hh"

namespace Pipal
{

  /**
   * \brief Solver class for the Pipal library.
   *
   * The Solver class provides an interface for solving the optimization problems using Frank E.
   * Curtis Pipal algorithm. It utilizes the Problem class to define the optimization problem and
   * implements various methods for solving it.
   * \tparam Real The floating-point type used for computations (e.g., float, double).
   */
  template <typename Real>
  class Solver
  {
  static_assert(std::is_floating_point_v<Real>,
     "Pipal::Solver<Real>: Real must be a floating-point type.");

  public:
    using ProblemPtr              = typename Problem<Real>::UniquePtr;
    using ObjectiveFunc           = typename ProblemWrapper<Real>::ObjectiveFunc;
    using ObjectiveGradientFunc   = typename ProblemWrapper<Real>::ObjectiveGradientFunc;
    using ConstraintsFunc         = typename ProblemWrapper<Real>::ConstraintsFunc;
    using ConstraintsJacobianFunc = typename ProblemWrapper<Real>::ConstraintsJacobianFunc;
    using LagrangianHessianFunc   = typename ProblemWrapper<Real>::LagrangianHessianFunc;
    using BoundsFunc              = typename ProblemWrapper<Real>::BoundsFunc;

  private:
    Counter          m_counter;    /*!< Internal counters for solver statistics. */
    Acceptance<Real> m_acceptance; /*!< Acceptance criteria for trial points. */
    Input<Real>      m_input;      /*!< Input structure for the solver. */
    Direction<Real>  m_direction;  /*!< Current search direction of the solver. */
    Iterate<Real>    m_iterate;    /*!< Current iterate of the solver. */
    Output<Real>     m_output;     /*!< Output class for managing solver output. */
    Parameter<Real>  m_parameter;  /*!< Internal parameters for the solver algorithm. */
    ProblemPtr       m_problem;    /*!< Problem object pointer. */

    // Some options for the solver
    bool m_verbose{false}; /*!< Verbosity flag. */

  public:
    /**
     * \brief Default constructor for the Pipal class.
     *
     * Initializes the solver with default values for the objective, gradient, constraints, and Jacobian
     * functions.
     */
    Solver() = default;

    /**
     * \brief Constructor for the Pipal class.
     *
     * Initializes the solver with the provided objective, gradient, constraints, and Jacobian functions.
     * \param[in] name Name of the optimization problem.
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
    Solver(std::string const & name, ObjectiveFunc const & objective, ObjectiveGradientFunc const &
      objective_gradient, ConstraintsFunc const & constraints, ConstraintsJacobianFunc const &
      constraints_jacobian, LagrangianHessianFunc const & lagrangian_hessian, BoundsFunc const &
      primal_lower_bounds, BoundsFunc const & primal_upper_bounds, BoundsFunc const &
      constraints_lower_bounds, BoundsFunc const & constraints_upper_bounds)
      : m_problem(std::make_unique<ProblemWrapper<Real>>(name, objective, objective_gradient,
        constraints, constraints_jacobian, lagrangian_hessian, primal_lower_bounds, primal_upper_bounds,
        constraints_lower_bounds, constraints_upper_bounds)) {}

    /**
     * \brief Constructor for the Pipal class (with a unique pointer to a Problem object).
     *
     * Initializes the solver with the provided unique pointer to a Problem object.
     * \param[in] problem The unique pointer to the Problem object to use.
     * \warning The pointer is moved into the solver, so it should not be used after this call.
     */
    Solver(ProblemPtr && problem) : m_problem(std::move(problem)) {}

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
     * \brief Set the problem to be solved using a unique pointer.
     *
     * This method allows the user to specify the problem to be solved using a unique pointer.
     * \param[in] problem The unique pointer to the problem to set.
     */
    void problem(ProblemPtr && problem) {this->m_problem = std::move(problem);}

    /**
     * \brief Get the problem being solved.
     * \return A reference to the problem.
     */
    Problem<Real> const & problem() const {return *this->m_problem;}

    /**
     * \brief Get the verbose mode.
     * \return The verbose mode.
     */
    [[nodiscard]] bool verbose_mode() const {return this->m_verbose;}

    /**
     * \brief Set the verbose mode.
     * \param[in] t_verbose The verbose mode.
     */
    void verbose_mode(bool const t_verbose) {this->m_verbose = t_verbose;}

    /**
     * \brief Get the algorithm mode.
     * \return The algorithm mode.
     */
    [[nodiscard]] Algorithm algorithm() const {return this->m_parameter.algorithm;}

    /**
     * \brief Set the algorithm mode.
     *
     * This method allows the user to specify the algorithm mode, which can be either \c CONSERVATIVE
     * or \c ADAPTIVE.
     * \param[in] t_algorithm The algorithm mode.
     */
    void algorithm(Algorithm const t_algorithm) {this->m_parameter.algorithm = t_algorithm;}

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
      PIPAL_ASSERT(t_tolerance > 0.0,
        "Pipal::Solver::tolerance(...): input value must be positive");
      this->m_parameter.opt_err_tol = t_tolerance;
    }

    /**
     * \brief Get the tolerance for convergence.
     * \return The tolerance for convergence.
     */
    Real tolerance() const {return this->m_parameter.opt_err_tol;}

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
      this->m_parameter.iter_max = t_max_iterations;
    }

    /**
     * \brief Get the maximum number of iterations.
     * \return The maximum number of iterations.
     */
    [[nodiscard]] Integer max_iterations() const {return this->m_parameter.iter_max;}

    /**
     * \brief Solves the optimization problem using the interior-point method.
     *
     * This method implements the interior-point algorithm to solve the optimization problem defined
     * by the objective function, constraints, and their respective gradients and Jacobians.
     * \param[in] x_guess Initial guess for the optimization variables.
     * \param[out] x_sol Solution vector.
     * \return True if the optimization was successful, false otherwise.
     */
    bool optimize(Vector<Real> const & x_guess, Vector<Real> & x_sol)
    {
      #define CMD "Pipal::Solver::optimize(...): "

      // Create alias for easier access
      Parameter<Real>  & p{this->m_parameter};
      Counter          & c{this->m_counter};
      Input<Real>      & i{this->m_input};
      Iterate<Real>    & z{this->m_iterate};
      Direction<Real>  & d{this->m_direction};
      Acceptance<Real> & a{this->m_acceptance};

      // Check that the problem is set
      PIPAL_ASSERT(this->m_problem.get() != nullptr,
        CMD "problem not set, use 'problem(...)' method to set it");

      // Get variable bounds
      Vector<Real> bl, bu;
      PIPAL_ASSERT(this->m_problem->primal_lower_bounds(bl),
        CMD "error in evaluating lower bounds on primal variables");
      PIPAL_ASSERT(this->m_problem->primal_upper_bounds(bu),
        CMD "error in evaluating upper bounds on primal variables");

      // Get constraint bounds
      Vector<Real> cl, cu;
      PIPAL_ASSERT(this->m_problem->constraints_lower_bounds(cl),
        CMD "error in evaluating lower bounds on constraints");
      PIPAL_ASSERT(this->m_problem->constraints_upper_bounds(cu),
        CMD "error in evaluating upper bounds on constraints");

      // Reset counters
      resetCounter(c);
      resetDirection<Real>(d, i);

      // Fill input structure
      buildInput<Real>(i, p, this->m_problem->name(), x_guess, bl, bu, cl, cu);
      buildIterate<Real>(z, p, i, c, *this->m_problem);

      // Print header and break line
      if (this->m_verbose) {this->m_output.printHeader(i, z); this->m_output.printBreak(c);}

      // Iterations loop
      while (!checkTermination<Real>(z, p, i, c))
      {
        // Print iterate
        if (this->m_verbose) {this->m_output.printIterate(c, z);}

        // Evaluate the step
        evalStep<Real>(d, p, i, c, z, a, *this->m_problem);

        // Print direction
        if (this->m_verbose) {this->m_output.printDirection(z, d);}

        lineSearch<Real>(a, p, i, c, z, d, *this->m_problem);

        // Print accepted
        if (this->m_verbose) {this->m_output.printAcceptance(a);}

        updateIterate<Real>(z, p, i, c, d, a, *this->m_problem);

        // Increment iteration counter
        incrementIterationCount(c);

        // Print break
        if (this->m_verbose) {this->m_output.printBreak(c);}
      }
      // Print footer and terminate
      if (this->m_verbose) {this->m_output.printFooter(p, i, c, z);}
      x_sol = x_guess;
      return true;

      #undef CMD
    }

  }; // class Solver

} // namespace Solver

#endif /* INCLUDE_PIPAL_SOLVER_HH */
