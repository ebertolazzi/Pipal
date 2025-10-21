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

// Quadratic program class
class QuadraticProgram : public Pipal::Problem
{
private:
  Matrix m_Q;
  Vector m_c;
  Matrix m_A;
  Vector m_b;

public:
  // Constructor
  QuadraticProgram(Matrix const & Q, Vector const & c, Matrix const & A, Vector const & b)
    : Problem("test_linear"), m_Q(Q), m_c(c), m_A(A), m_b(b) {}

  // Objective function
  bool objective(Vector const & x, Real & out) const override
  {
    out = 0.5 * x.dot(this->m_Q * x) + this->m_c.dot(x);
    return std::isfinite(out);
  }

  // Gradient of the objective function
  bool objective_gradient(Vector const & x, Vector & out) const override
  {
    out = this->m_Q * x + this->m_c;
    return out.allFinite();
  }

  // Hessian of the objective function
  bool objective_hessian(Vector const & /*x*/, Matrix & out) const override
  {
    out = this->m_Q;
    return out.allFinite();
  }

  // Constraints function
  bool constraints(Vector const & x, Vector & out) const override
  {
    out = this->m_A * x - this->m_b;
    return out.allFinite();
  }

  // Jacobian of the constraints function
  bool constraints_jacobian(Vector const & /*x*/, Matrix & out) const override
  {
    out = this->m_A;
    return out.allFinite();
  }

  // Hessian of the Lagrangian function
  bool lagrangian_hessian(Vector const & /*x*/, Vector const & /*z*/, Matrix & out) const override
  {
    out = this->m_Q;
    return out.allFinite();
  }

  // Lower bounds on the primal variables
  bool primal_lower_bounds(Vector & out) const override
  {
    out = Vector::Constant(this->m_c.size(), -1e+10);
    return out.allFinite();
  }

  // Upper bounds on the primal variables
  bool primal_upper_bounds(Vector & out) const override
  {
    out = Vector::Constant(this->m_c.size(), 1e+10);
    return out.allFinite();
  }

  // Lower bounds on the constraints
  bool constraints_lower_bounds(Vector & out) const override
  {
    out = Vector::Constant(this->m_b.size(), -1e+10);
    return out.allFinite();
  }

  // Upper bounds on the constraints
  bool constraints_upper_bounds(Vector & out) const override
  {
    out = Vector::Zero(this->m_b.size());
    return out.allFinite();
  }
};

// Test fixture for Pipal with QuadraticProgram
class LinearRegression : public testing::Test {
protected:
  static constexpr Integer N{2};
  static constexpr Integer M{5};

  std::unique_ptr<QuadraticProgram> problem;
  Vector x_guess, sol;

  void SetUp() override {
    Matrix Q(2.0 * Matrix::Identity(N, N));
    Vector c(N); c << -2.0, -5.0;
    Matrix A(M, N);
    A <<
       1.0,  2.0,
      -1.0,  2.0,
      -1.0, -2.0,
       1.0,  0.0,
       0.0,  1.0;
    Vector b(M); b << 6.0, 2.0, 2.0, 3.0, 2.0;

    x_guess.resize(N); x_guess << 0.5, 0.5;
    sol.resize(N); sol << 1.4, 1.7;

    problem = std::make_unique<QuadraticProgram>(Q, c, A, b);
  }
};

TEST_F(LinearRegression, Problem) {
  Pipal::Solver solver(std::move(problem));
  solver.verbose_mode(VERBOSE);
  solver.tolerance(SOLVER_TOLERANCE);
  solver.max_iterations(MAX_ITERATIONS);

  Vector x_sol;
  EXPECT_TRUE(solver.optimize(x_guess, x_sol));
  EXPECT_TRUE(x_sol.isApprox(sol, APPROX_TOLERANCE));
}

TEST_F(LinearRegression, ProblemWrapper) {
  Pipal::ProblemWrapper problem_wrapper("test_linear",
    [this] (const Vector & x, Real & out) {return this->problem->objective(x, out);},
    [this] (const Vector & x, Vector & out) {return this->problem->objective_gradient(x, out);},
    [this] (const Vector & x, Matrix & out) {return this->problem->objective_hessian(x, out);},
    [this] (const Vector & x, Vector & out) {return this->problem->constraints(x, out);},
    [this] (const Vector & x, Matrix & out) {return this->problem->constraints_jacobian(x, out);},
    [this] (const Vector & x, const Vector & z, Matrix & out) {return this->problem->lagrangian_hessian(x, z, out);},
    [this] (Vector & out) {return this->problem->primal_lower_bounds(out);},
    [this] (Vector & out) {return this->problem->primal_upper_bounds(out);},
    [this] (Vector & out) {return this->problem->constraints_lower_bounds(out);},
    [this] (Vector & out) {return this->problem->constraints_upper_bounds(out);}
  );

 Pipal::Solver solver(problem_wrapper.name(),
   problem_wrapper.objective(), problem_wrapper.objective_gradient(), problem_wrapper.objective_hessian(),
   problem_wrapper.constraints(), problem_wrapper.constraints_jacobian(), problem_wrapper.lagrangian_hessian(),
   problem_wrapper.primal_lower_bounds(), problem_wrapper.primal_upper_bounds(),
   problem_wrapper.constraints_lower_bounds(), problem_wrapper.constraints_upper_bounds()
 );
 solver.verbose_mode(VERBOSE);
 solver.tolerance(SOLVER_TOLERANCE);
 solver.max_iterations(MAX_ITERATIONS);
 Vector x_sol;
 EXPECT_TRUE(solver.optimize(x_guess, x_sol));
 EXPECT_TRUE(x_sol.isApprox(sol, APPROX_TOLERANCE));
}
