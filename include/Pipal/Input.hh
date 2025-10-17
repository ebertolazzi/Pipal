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

#ifndef INCLUDE_PIPAL_INPUT_HH
#define INCLUDE_PIPAL_INPUT_HH

// Pipal includes
#include "Pipal/Defines.hh"
#include "Pipal/Problem.hh"
#include "Pipal/Parameter.hh"

// Standard libraries
#include <type_traits>
#include <string>

// Eigen library
#include <Eigen/Dense>

namespace Pipal
{
  /**
   * \brief Input structure holding all the data defining the optimization problem.
   * \tparam Real The floating-point type.
   * \tparam Integer The integer type.
   */
  template<typename Real, typename Integer>
  struct Input
  {
    static_assert(std::is_floating_point<Real>::value,
      "Pipal::Input: template argument 'Real' must be a floating-point type.");
    static_assert(std::is_integral<Integer>::value,
      "Pipal::Input: template argument 'Integer' must be an integer type.");

    using std::string;
    using Vector  = Eigen::Vector<Real, Eigen::Dynamic>;
    using Indices = Eigen::Array<Integer, Eigen::Dynamic, 1>;

    using ObjectiveFunc           = std::function<bool(Vector const &, Real &)>;
    using ObjectiveGradientFunc   = std::function<bool(Vector const &, Vector &)>;
    using ObjectiveHessianFunc    = std::function<bool(Vector const &, Matrix &)>;
    using ConstraintsFunc         = std::function<bool(Vector const &, Vector &)>;
    using ConstraintsJacobianFunc = std::function<bool(Vector const &, Vector const &, Matrix &)>;
    using LagrangianHessianFunc   = std::function<bool(Vector const &, Vector const &, Matrix &)>;
    using BoundsFunc              = std::function<bool(Vector &)>;

    string  id; // Problem identity
    Integer n0; // Number of original formulation variables
    Indices I1; // Indices of free variables
    Indices I2; // Indices of fixed variables
    Indices I3; // Indices of lower bounded variables
    Indices I4; // Indices of upper bounded variables
    Indices I5; // Indices of lower and upper bounded variables
    Indices I6; // Indices of equality constraints
    Indices I7; // Indices of lower bounded constraints
    Indices I8; // Indices of upper bounded constraints
    Indices I9; // Indices of lower and upper bounded constraints
    Vector  x0; // Initial guess for the primal variables
    Vector  b2; // Right-hand side of fixed variables
    Vector  l3; // Right-hand side of lower bounded variables
    Vector  u4; // Right-hand side of upper bounded variables
    Vector  l5; // Right-hand side of lower half of lower and upper bounded variables
    Vector  u5; // Right-hand side of upper half of lower and upper bounded variables
    Vector  b6; // Right-hand side of equality constraints
    Vector  l7; // Right-hand side of lower bounded constraints
    Vector  u8; // Right-hand side of upper bounded constraints
    Vector  l9; // Right-hand side of lower half of lower and upper bounded constraints
    Vector  u9; // Right-hand side of upper half of lower and upper bounded constraints
    Integer n1; // Number of free variables
    Integer n2; // Number of fixed variables
    Integer n3; // Number of lower bounded variables
    Integer n4; // Number of upper bounded variables
    Integer n5; // Number of lower and upper bounded variables
    Integer n6; // Number of equality constraints
    Integer n7; // Number of lower bounded constraints
    Integer n8; // Number of upper bounded constraints
    Integer n9; // Number of lower and upper bounded constraints
    Integer nV; // Number of variables
    Integer nI; // Number of inequality constraints
    Integer nE; // Number of equality constraints
    Integer nA; // Size of primal-dual matrix
    Integer vi; // Counter for invalid bounds

    ObjectiveFunc           f_orig; // original objective
    ConstraintsFunc         c_orig; // original constraints
    ObjectiveGradientFunc   g_orig; // original gradient of objective
    ConstraintsJacobianFunc J_orig; // original jacobian of constraints
    LagrangianHessianFunc   H_orig; // original hessian of the lagrangian

    // Constructor
    void Input(Parameter<Real> const & p, string const & name, ObjectiveFunc const & f_orig,
      ConstraintsFunc const & c_orig, ObjectiveGradientFunc const & g_orig,
      ConstraintsJacobianFunc const & J_orig, LagrangianHessianFunc const & H_orig,
      Vector const & x0, Vector const & bl, Vector const & bu, Vector const & cl, Vector const & cu)
    {
      #define CMD "Pipal::Solver::fill_input(...): "

      using Mask = Eigen::Array<bool, Eigen::Dynamic, 1>;

      // Get variable bounds
      Vector bl, bu;
      PIPAL_ASSERT(problem.primal_lower_bounds(bl),
        CMD "error in evaluating lower bounds on primal variables");
      PIPAL_ASSERT(problem.primal_upper_bounds(bu),
        CMD "error in evaluating upper bounds on primal variables");

      // Get constraint bounds
      Vector cl, cu;
      PIPAL_ASSERT(problem.constraints_lower_bounds(cl),
        CMD "error in evaluating lower bounds on constraints");
      PIPAL_ASSERT(problem.constraints_upper_bounds(cu),
        CMD "error in evaluating upper bounds on constraints");

      // Set problem identity
      this->id = name;

      // store function pointers to original problem functions
      this->f_orig = f_orig;
      this->c_orig = c_orig;
      this->g_orig = g_orig;
      this->J_orig = J_orig;
      this->H_orig = H_orig;

      // Set number of original formulation variables
      this->n0 = x0.size();

      // Find indices sets
      Real const tolerance{Eigen::NumTraits<Real>::epsilon()};
      Mask const cond_bl(bl.array() <= -params.rhs_bnd);
      Mask const cond_bu(bu.array() >=  params.rhs_bnd);
      Mask const cond_cl(cl.array() <= -params.rhs_bnd);
      Mask const cond_cu(cu.array() >=  params.rhs_bnd);
      Mask const cond_bq((bl.array() - bu.array()).abs() > tolerance);
      Mask const cond_cq((cl.array() - cu.array()).abs() > tolerance);

      this->I1 = this->find(cond_bl && cond_bu);
      this->I2 = this->find(cond_bq);
      this->I3 = this->find(!cond_bl &&  cond_bu);
      this->I4 = this->find(cond_bl && !cond_bu);
      this->I5 = this->find(!cond_bl && !cond_bu && !cond_bq);
      this->I6 = this->find(cond_cq);
      this->I7 = this->find(!cond_cl &&  cond_cu);
      this->I8 = this->find(cond_cl && !cond_cu);
      this->I9 = this->find(!cond_cl && !cond_cu && !cond_cq);

      // Set right-hand side values
      this->b2 = bl(this->I2);
      this->l3 = bl(this->I3);
      this->u4 = bu(this->I4);
      this->l5 = bl(this->I5);
      this->u5 = bu(this->I5);
      this->b6 = cl(this->I6);
      this->l7 = cl(this->I7);
      this->u8 = cu(this->I8);
      this->l9 = cl(this->I9);
      this->u9 = cu(this->I9);

      // Set sizes of indices sets
      this->n1 = this->I1.count();
      this->n2 = this->I2.count();
      this->n3 = this->I3.count();
      this->n4 = this->I4.count();
      this->n5 = this->I5.count();
      this->n6 = this->I6.count();
      this->n7 = this->I7.count();
      this->n8 = this->I8.count();
      this->n9 = this->I9.count();

      // Initialize number of invalid bounds
      this->vi = 0;

      // Count invalid bounds
      if (this->n2 > 0) {
        this->vi += (this->b2.array() <= -params.rhs_bnd).count();
        this->vi += (this->b2.array() >= params.rhs_bnd).count();
      }
      if (this->n3 > 0) {
        this->vi += (this->l3.array() >= params.rhs_bnd).count();
      }
      if (this->n4 > 0) {
        this->vi += (this->u4.array() <= -params.rhs_bnd).count();
      }
      if (this->n5 > 0) {
        this->vi += (this->l5.array() >= params.rhs_bnd).count();
        this->vi += (this->u5.array() <= -params.rhs_bnd).count();
        this->vi += (this->l5.array() > this->u5.array()).count();
      }
      if (this->n6 > 0) {
        this->vi += (this->b6.array() <= -params.rhs_bnd).count();
        this->vi += (this->b6.array() >= params.rhs_bnd).count();
      }
      if (this->n7 > 0) {
        this->vi += (this->l7.array() >= params.rhs_bnd).count();
      }
      if (this->n8 > 0) {
        this->vi += (this->u8.array() <= -params.rhs_bnd).count();
      }
      if (this->n9 > 0) {
        this->vi += (this->l9.array() >= params.rhs_bnd).count();
        this->vi += (this->u9.array() <= -params.rhs_bnd).count();
        this->vi += (this->l9.array() > this->u9.array()).count();
      }

      // Set number of variables and constraints
      this->nV = this->n1 + this->n3 + this->n4 + this->n5;
      this->nI = this->n3 + this->n4 + 2*this->n5 + this->n7 + this->n8 + 2*this->n9;
      this->nE = this->n6;

      // Set size of primal-dual matrix
      this->nA = this->nV + 3*this->nE + 3*this->nI;

      // Set initial point
      this->x0.resize(this->nV);
      this->x0 << x0(this->I1), x0(this->I3), x0(this->I4), x0(this->I5);

      #undef CMD
    }

  }; // struct Input

} // namespace Pipal

#endif /* INCLUDE_PIPAL_INPUT_HH */
