/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Pipal project is distributed under the MIT License.                                       *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// STL includes
#include <memory>
#include <limits>

// GTest library
#include <gtest/gtest.h>

// Pipal includes
#include "Pipal.hh"

using Pipal::Integer;
using Pipal::Real;
using Pipal::Vector;
using Pipal::Matrix;

constexpr bool VERBOSE{true};
constexpr Real SOLVER_TOLERANCE{1.0e-9};
constexpr Real APPROX_TOLERANCE{1.0e-6};
constexpr Integer MAX_ITERATIONS{100};

TEST(Test1, ProblemWrapper) {
  Pipal::ProblemWrapper problem_wrapper("test_basic",
    [] (const Vector & x, Real & out) { // Objective function
      out = 100.0*std::pow(x(1) - x(0)*x(0), 2.0) + std::pow(1 - x(0), 2.0);
      return std::isfinite(out);
    },
    [] (const Vector & x, Vector & out) { // Gradient of the objective function
      out.resize(2);
      out << -400.0*x(0)*(x(1)-x(0)*x(0)) - 2.0*(1 - x(0)), 200.0*(x(1) - x(0)*x(0));
      return out.allFinite();
    },
    [] (const Vector &, Matrix & out) { // Hessian of the objective function
      out.resize(0,0);
      return true;
    },
    [] (const Vector &, Vector & out) { // Constraints function
      out.resize(0);
      return true;
    },
    [] (const Vector &, Matrix & out) { // Jacobian of the constraints function
      out.resize(0,0);
      return true;
    },
    [] (const Vector & x, const Vector &, Matrix & out) ->bool { // Hessian of the Lagrangian
      out.resize(2,2);
      out <<
        1200.0*x(0)*x(0) - 400.0*x(1) + 2, -400.0*x(0),
        -400.0*x(0), 200.0;
      return out.allFinite();
    },
    [] (Vector & out) {out.setConstant(2, -1.0); return true;}, // Lower bounds on the primal variables
    [] (Vector & out) {out.setConstant(2, +1.0); return true;}, // Upper bounds on the primal variables
    [] (Vector & out) {out.resize(0); return true;}, // Lower bounds on the constraints
    [] (Vector & out) {out.resize(0); return true;}  // Upper bounds on the constraints
  );

  Pipal::Solver solver(
    problem_wrapper.name(), problem_wrapper.objective(), problem_wrapper.objective_gradient(),
    problem_wrapper.objective_hessian(), problem_wrapper.constraints(), problem_wrapper.constraints_jacobian(),
    problem_wrapper.lagrangian_hessian(), problem_wrapper.primal_lower_bounds(), problem_wrapper.primal_upper_bounds(),
    problem_wrapper.constraints_lower_bounds(), problem_wrapper.constraints_upper_bounds()
  );
  solver.verbose_mode(VERBOSE);
  solver.tolerance(SOLVER_TOLERANCE);
  solver.max_iterations(MAX_ITERATIONS);
  Vector x_sol, x_guess;
  x_guess.setZero(2);
  EXPECT_TRUE(solver.optimize(x_guess, x_sol));
   EXPECT_TRUE(x_sol.isApprox(x_sol, APPROX_TOLERANCE));
}
