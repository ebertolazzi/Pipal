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

// GTest library
#include <gtest/gtest.h>

// IPsolPipalver includes
#include "Pipal.hh"

using Pipal::Integer;

constexpr bool VERBOSE{true};
constexpr double SOLVER_TOLERANCE{5.0e-5};
constexpr double APPROX_TOLERANCE{1.0e-4};
constexpr Integer MAX_ITERATIONS{100};

// Quadratic program class
template<typename Real, Integer N, Integer M>
class QuadraticProgram : public Pipal::Problem<Real, N, M>
{
public:
  using typename Pipal::Problem<Real, N, M>::VectorN;
  using typename Pipal::Problem<Real, N, M>::VectorM;
  using typename Pipal::Problem<Real, N, M>::MatrixH;
  using typename Pipal::Problem<Real, N, M>::MatrixJ;

private:
  MatrixH m_Q;
  VectorN m_c;
  MatrixJ m_A;
  VectorM m_b;

public:
  // Constructor
  QuadraticProgram(MatrixH const & Q, VectorN const & c, MatrixJ const & A, VectorM const & b)
    : m_Q(Q), m_c(c), m_A(A), m_b(b) {}

  // Objective function
  Real objective(VectorN const & x) const override
  {
    return 0.5 * x.dot(m_Q * x) + m_c.dot(x);
  }

  // Gradient of the objective function
  VectorN objective_gradient(VectorN const & x) const override
  {
    return m_Q * x + m_c;
  }

  // Hessian of the objective function
  MatrixH objective_hessian(VectorN const & /*x*/) const override
  {
    return m_Q;
  }

  // Constraints function
  VectorM constraints(VectorN const & x) const override
  {
    return m_A * x - m_b;
  }

  // Jacobian of the constraints function
  MatrixJ constraints_jacobian(VectorN const & /*x*/, VectorM const & /*z*/) const override
  {
    return m_A;
  }

  // Hessian of the Lagrangian function
  MatrixH lagrangian_hessian(VectorN const & /*x*/, VectorM const & /*z*/) const override
  {
    return m_Q;
  }
};

// Test fixture for Pipal with QuadraticProgram
class LinearRegression : public testing::Test {
protected:
  using TestType = double;

  static constexpr Integer N{2};
  static constexpr Integer M{5};

  using VectorN = typename Pipal::Problem<TestType, N, M>::VectorN;
  using VectorM = typename Pipal::Problem<TestType, N, M>::VectorM;
  using MatrixJ = typename Pipal::Problem<TestType, N, M>::MatrixJ;
  using MatrixH = typename Pipal::Problem<TestType, N, M>::MatrixH;

  std::unique_ptr<QuadraticProgram<TestType, N, M>> problem;
  VectorN x_guess, sol;

  void SetUp() override {
    MatrixH Q(2.0 * MatrixH::Identity(N, N));
    VectorN c(N); c << -2.0, -5.0;
    MatrixJ A(M, N);
    A <<  1.0,  2.0,
       -1.0,  2.0,
       -1.0, -2.0,
        1.0,  0.0,
        0.0,  1.0;
    VectorM b(M); b << 6.0, 2.0, 2.0, 3.0, 2.0;

    x_guess.resize(N); x_guess << 0.5, 0.5;
    sol.resize(N); sol << 1.4, 1.7;

    problem = std::make_unique<QuadraticProgram<TestType, N, M>>(Q, c, A, b);
  }
};

TEST_F(LinearRegression, DISABLED_Problem_BFGS) {
  Pipal::Solver<TestType, N, M> solver(std::move(problem));
  solver.verbose_mode(VERBOSE);
  solver.tolerance(SOLVER_TOLERANCE);
  solver.max_iterations(MAX_ITERATIONS);

  VectorN x_sol;
  EXPECT_TRUE(solver.solve(x_guess, x_sol));
  EXPECT_TRUE(x_sol.isApprox(sol, APPROX_TOLERANCE));
}

TEST_F(LinearRegression, DISABLED_Problem_Newton) {
  Pipal::Solver<TestType, N, M> solver(std::move(problem));
  solver.verbose_mode(VERBOSE);
  solver.tolerance(SOLVER_TOLERANCE);
  solver.max_iterations(MAX_ITERATIONS);

  VectorN x_sol;
  EXPECT_TRUE(solver.solve(x_guess, x_sol));
  EXPECT_TRUE(x_sol.isApprox(sol, APPROX_TOLERANCE));
}

TEST_F(LinearRegression, Problem_Steepest) {
  Pipal::Solver<TestType, N, M> solver(std::move(problem));
  solver.verbose_mode(VERBOSE);
  solver.tolerance(SOLVER_TOLERANCE);
  solver.max_iterations(MAX_ITERATIONS);

  VectorN x_sol;
  EXPECT_TRUE(solver.solve(x_guess, x_sol));
  EXPECT_TRUE(x_sol.isApprox(sol, APPROX_TOLERANCE));
}

TEST_F(LinearRegression, DISABLED_ProblemWrapper_BFGS) {
  Pipal::ProblemWrapper<TestType, N, M> problem_wrapper(
    [this] (const VectorN & x) {return this->problem->objective(x);},
    [this] (const VectorN & x) {return this->problem->objective_gradient(x);},
    [this] (const VectorN & x) {return this->problem->objective_hessian(x);},
    [this] (const VectorN & x) {return this->problem->constraints(x);},
    [this] (const VectorN & x, const VectorM & z) {return this->problem->constraints_jacobian(x, z);},
    [this] (const VectorN & x, const VectorM & z) {return this->problem->lagrangian_hessian(x, z);}
  );

  Pipal::Solver<TestType, N, M> solver(
    problem_wrapper.objective(), problem_wrapper.objective_gradient(), problem_wrapper.objective_hessian(),
    problem_wrapper.constraints(), problem_wrapper.constraints_jacobian(), problem_wrapper.lagrangian_hessian()
  );
  solver.verbose_mode(VERBOSE);
  solver.tolerance(SOLVER_TOLERANCE);
  solver.max_iterations(MAX_ITERATIONS);
  VectorN x_sol;
  EXPECT_TRUE(solver.solve(x_guess, x_sol));
  EXPECT_TRUE(x_sol.isApprox(sol, APPROX_TOLERANCE));
}

TEST_F(LinearRegression, DISABLED_ProblemWrapper_Newton) {
  Pipal::ProblemWrapper<TestType, N, M> problem_wrapper(
    [this] (const VectorN & x) {return problem->objective(x);},
    [this] (const VectorN & x) {return problem->objective_gradient(x);},
    [this] (const VectorN & x) {return problem->objective_hessian(x);},
    [this] (const VectorN & x) {return problem->constraints(x);},
    [this] (const VectorN & x, const VectorM & z) {return problem->constraints_jacobian(x, z);},
    [this] (const VectorN & x, const VectorM & z) {return problem->lagrangian_hessian(x, z);}
  );

  Pipal::Solver<TestType, N, M> solver(
    problem_wrapper.objective(), problem_wrapper.objective_gradient(), problem_wrapper.objective_hessian(),
    problem_wrapper.constraints(), problem_wrapper.constraints_jacobian(), problem_wrapper.lagrangian_hessian()
  );
  solver.verbose_mode(VERBOSE);
  solver.tolerance(SOLVER_TOLERANCE);
  solver.max_iterations(MAX_ITERATIONS);
  VectorN x_sol;
  EXPECT_TRUE(solver.solve(x_guess, x_sol));
  EXPECT_TRUE(x_sol.isApprox(sol, APPROX_TOLERANCE));
}

TEST_F(LinearRegression, ProblemWrapper_Steepest) {
  Pipal::ProblemWrapper<TestType, N, M> problem_wrapper(
    [this] (const VectorN & x) {return problem->objective(x);},
    [this] (const VectorN & x) {return problem->objective_gradient(x);},
    [this] (const VectorN & x) {return problem->objective_hessian(x);},
    [this] (const VectorN & x) {return problem->constraints(x);},
    [this] (const VectorN & x, const VectorM & z) {return problem->constraints_jacobian(x, z);},
    [this] (const VectorN & x, const VectorM & z) {return problem->lagrangian_hessian(x, z);}
  );

  Pipal::Solver<TestType, N, M> solver(
    problem_wrapper.objective(), problem_wrapper.objective_gradient(), problem_wrapper.objective_hessian(),
    problem_wrapper.constraints(), problem_wrapper.constraints_jacobian(), problem_wrapper.lagrangian_hessian()
  );
  solver.verbose_mode(VERBOSE);
  solver.tolerance(SOLVER_TOLERANCE);
  solver.max_iterations(MAX_ITERATIONS);
  VectorN x_sol;
  EXPECT_TRUE(solver.solve(x_guess, x_sol));
  EXPECT_TRUE(x_sol.isApprox(sol, APPROX_TOLERANCE));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
