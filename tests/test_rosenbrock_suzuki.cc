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
constexpr Real APPROX_TOLERANCE{1.0e-3};
constexpr Integer MAX_ITERATIONS{100};

TEST(Test1, ProblemWrapper) {
  Pipal::Solver solver("test_rosenbrock_suzuki",
    [] (const Vector & x, Real & out) { // Objective function
      out = x(0)*x(0) + x(1)*x(1) + 2.0*x(2)*x(2) + x(3)*x(3) - 5.0*x(0) - 5.0*x(1) - 21.0*x(2) + 7.0*x(3);
      return std::isfinite(out);
    },
    [] (const Vector & x, Vector & out) { // Gradient of the objective function
      out.resize(4);
      out << 2.0*x(0)-5.0, 2.0*x(1)-5.0, 4.0*x(2)-21.0, 2.0*x(3)+7.0;
      return out.allFinite();
    },
    [] (const Vector &, Matrix & out) { // Hessian of the objective function
      out.resize(4,4);
      out.setZero();
      out.diagonal() << 2.0, 2.0, 4.0, 2.0;
      return true;
    },
    [] (const Vector & x, Vector & out) { // Constraints function
      out.resize(4);
      out <<
        (x(0)*x(0) + x(1)*x(1) + x(2)*x(2) + x(3)*x(3) + x(0) - x(1) + x(2) - x(3)) - 8.0,
        (x(0)*x(0) + 2.0*x(1)*x(1) + x(2)*x(2) + 2.0*x(3)*x(3) - x(0) - x(3)) - 10.0,
        (2.0*x(0)*x(0) + x(1)*x(1) + x(2)*x(2) + 2.0*x(0) - x(1) - x(3)) - 5.0,
        (x(0) + x(1) + x(2) + x(3) - 5.0);
      return out.allFinite();
    },
    [] (const Vector & x, Matrix & out) { // Jacobian of the constraints function
      out.resize(4,4);
      out <<
        2.0*x(0) + 1.0, 2.0*x(1) - 1.0, 2.0*x(2) + 1.0, 2.0*x(3) - 1.0,
        2.0*x(0) - 1.0, 4.0*x(1),       2.0*x(2),       4.0*x(3) - 1.0,
        4.0*x(0) + 2.0, 2.0*x(1) - 1.0, 2.0*x(2),       -1.0,
        1.0,            1.0,            1.0,             1.0;
      return out.allFinite();
    },
    [] (const Vector &, const Vector &z, Matrix & out) ->bool { // Hessian of the Lagrangian
      out.resize(4,4);
      Vector diag0(4); diag0 << 2.0, 2.0, 4.0, 2.0;
      Vector diag1(4); diag1 << 2.0, 2.0, 2.0, 2.0;
      Vector diag2(4); diag2 << 2.0, 4.0, 2.0, 4.0;
      Vector diag3(4); diag3 << 4.0, 2.0, 2.0, 0.0;
      out = diag0.asDiagonal() + z(0)*diag1.asDiagonal() + z(1)*diag2.asDiagonal() + z(2)*diag3.asDiagonal();
      return out.allFinite();
    },
    [] (Vector & out) {out.setConstant(4, -INFINITY); out(0) = -3000.0; return true;}, // Lower bounds on the primal variables
    [] (Vector & out) {out.setConstant(4, INFINITY); return true;}, // Upper bounds on the primal variables
    [] (Vector & out) {out.setConstant(4, -INFINITY); out(3) = 0.0; return true;}, // Lower bounds on the constraints
    [] (Vector & out) {out.setConstant(4, 0.0); return true;}  // Upper bounds on the constraints
  );
  solver.algorithm(Pipal::Algorithm::ADAPTIVE);
  solver.verbose_mode(VERBOSE);
  solver.tolerance(SOLVER_TOLERANCE);
  solver.max_iterations(MAX_ITERATIONS);
  Vector x_sol(4), x_guess(4), x_opt(4);
  x_guess << -4000.0, 1.0, 1.0, 1.0;
  x_opt << 1.0, 1.0, 1.0, 1.0;
  EXPECT_TRUE(solver.optimize(x_guess, x_sol));
  EXPECT_TRUE(x_sol.isApprox(x_opt, APPROX_TOLERANCE));
}
