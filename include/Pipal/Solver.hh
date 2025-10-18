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
#include <iostream>
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
//#include "Pipal/Output.hh"
#include "Pipal/Parameter.hh"
#include "Pipal/Problem.hh"

namespace Pipal
{

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

    using ProblemPtr              = typename Problem<Real>::UniquePtr;
    using ObjectiveFunc           = typename ProblemWrapper<Real>::ObjectiveFunc;
    using ObjectiveGradientFunc   = typename ProblemWrapper<Real>::ObjectiveGradientFunc;
    using ObjectiveHessianFunc    = typename ProblemWrapper<Real>::ObjectiveHessianFunc;
    using ConstraintsFunc         = typename ProblemWrapper<Real>::ConstraintsFunc;
    using ConstraintsJacobianFunc = typename ProblemWrapper<Real>::ConstraintsJacobianFunc;
    using LagrangianHessianFunc   = typename ProblemWrapper<Real>::LagrangianHessianFunc;
    using BoundsFunc              = typename ProblemWrapper<Real>::BoundsFunc;

  private:
    Acceptance m_acceptance; /*!< Acceptance criteria for trial points. */
    Counter         m_counter; /*!< Internal counters for solver statistics. */
    Input     m_input;  /*!< Input structure for the solver. */
    Direction m_direction; /*!< Current search direction of the solver. */
    Iterate   m_iterate; /*!< Current iterate of the solver. */
    //Output<Real>    m_output;  /*!< Output class for managing solver output. */
    Parameter m_parameter;   /*!< Internal parameters for the solver algorithm. */
    ProblemPtr      m_problem; /*!< Problem object pointer. */

    // Some options for the solver
    bool      m_verbose{false};      /*!< Verbosity flag. */
    Integer   m_max_iterations{100}; /*!< Maximum number of iterations. */
    Real      m_tolerance{1.0e-10};  /*!< Optimality tolerance. */

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
    bool verbose_mode() {return this->m_verbose;}

    /**
     * \brief Set the verbose mode.
     * \param[in] t_verbose The verbose mode.
     */
    void verbose_mode(bool const t_verbose) {this->m_verbose = t_verbose;}

    /**
     * \brief Get the algorithm mode.
     * \return The algorithm mode.
     */
    //FIXME: Algorithm algorithm() const {return this->m_parameters.algorithm;}

    /**
     * \brief Set the algorithm mode.
     *
     * This method allows the user to specify the algorithm mode, which can be either \c CONSERVATIVE
     * or \c ADAPTIVE.
     * \param[in] t_algorithm The algorithm mode.
     */
    //FIXME: void algorithm(Algorithm const t_algorithm) {this->m_parameters.algorithm = t_algorithm;}

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

      // Create alias for easier access
      Parameter  & p{this->m_parameter};
      Counter          & c{this->m_counter};
      Input      & i{this->m_input};
      Iterate    & z{this->m_iterate};
      Direction  & d{this->m_direction};
      Acceptance & a{this->m_acceptance};
      Problem<Real>    * problem{this->m_problem.get()};

      // Get variable bounds
      Vector bl(i.nV), bu(i.nV);
      PIPAL_ASSERT(problem->primal_lower_bounds(bl),
        CMD "error in evaluating lower bounds on primal variables");
      PIPAL_ASSERT(problem->primal_upper_bounds(bu),
        CMD "error in evaluating upper bounds on primal variables");

      // Get constraint bounds
      Vector cl(i.nI), cu(i.nI);
      PIPAL_ASSERT(problem->constraints_lower_bounds(cl),
        CMD "error in evaluating lower bounds on constraints");
      PIPAL_ASSERT(problem->constraints_upper_bounds(cu),
        CMD "error in evaluating upper bounds on constraints");

      // Reset counters
      c.buildCounter();

      // Check that the problem is set
      PIPAL_ASSERT(problem != nullptr,
        CMD "problem not set, use 'problem(...)' method to set it");

      // Fill input structure
      std::string name{"Pipal Problem FIXME"};
      i = buildInput(p, name, problem->objective(), problem->constraints(),
        problem->objective_gradient(), problem->constraints_jacobian(),
        problem->lagrangian_hessian(), x_guess, bl, bu, cl, cu);

      // Print header and break line
      //if (this->m_verbose) {output = Output(input) o.print_header(); o.print_break(c);}

      // Iterations loop
      while (!checkTermination(z, p, i, c))
      {
        // Print iterate
        //if (this->m_verbose) {o.print_iterate(c, i);}

        // Evaluate the step
        evalStep(d, p, i, c, z, a);

        // Print direction
        //if (this->m_verbose) {o.print_direction(d);}

        lineSearch(a, p, i, c, z, d);

        // Print accepted
        //if (this->m_verbose) {o.print_acceptance(a);}

        updateIterate(z, p, i, c, d, a);

        // Increment iteration counter
        incrementIterationCount(c);

        // Print break
        //if (this->m_verbose) {o.print_break(c);}
      }

      // Print footer and terminate
      //if (this->m_verbose) {o.print_footer(p, i, c, i);o.terminate();}
      x_sol = x_guess;

      return true;

      #undef CMD
    }

  }; // class Solver

} // namespace Solver

#endif /* INCLUDE_PIPAL_SOLVER_HH */
