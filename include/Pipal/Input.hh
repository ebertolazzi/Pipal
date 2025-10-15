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
#include "Pipal.hh"
#include "Pipal/Problem.hh"
#include "Pipal/Parameter.hh"

// Standard libraries
#include <type_traits>

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

    using Vector  = Eigen::Vector<Real, Eigen::Dynamic>;
    using Indices = Eigen::Array<Integer, Eigen::Dynamic, 1>;

    Indices i1; /*!< Indices of free variables. */
    Indices i2; /*!< Indices of fixed variables. */
    Indices i3; /*!< Indices of lower bounded variables. */
    Indices i4; /*!< Indices of upper bounded variables. */
    Indices i5; /*!< Indices of lower and upper bounded variables. */
    Indices i6; /*!< Indices of equality constraints. */
    Indices i7; /*!< Indices of lower bounded constraints. */
    Indices i8; /*!< Indices of upper bounded constraints. */
    Indices i9; /*!< Indices of lower and upper bounded constraints. */
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

    /**
     * \brief Fill the input structure with problem data.
     * \param[in] x0 Initial guess for the optimization variables.
     */
    void fill_input(Problem<Real> const & problem, Parameter<Real> const & params, Vector const & x0)
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

      // Find indices sets
      Real const tolerance{Eigen::NumTraits<Real>::epsilon()};
      Mask const cond_bl(bl.array() <= -params.rhs_bnd);
      Mask const cond_bu(bu.array() >=  params.rhs_bnd);
      Mask const cond_cl(cl.array() <= -params.rhs_bnd);
      Mask const cond_cu(cu.array() >=  params.rhs_bnd);
      Mask const cond_bq((bl.array() - bu.array()).abs() > tolerance);
      Mask const cond_cq((cl.array() - cu.array()).abs() > tolerance);

      this->i1 = this->find(cond_bl && cond_bu);
      this->i2 = this->find(cond_bq);
      this->i3 = this->find(!cond_bl &&  cond_bu);
      this->i4 = this->find(cond_bl && !cond_bu);
      this->i5 = this->find(!cond_bl && !cond_bu && !cond_bq);
      this->i6 = this->find(cond_cq);
      this->i7 = this->find(!cond_cl &&  cond_cu);
      this->i8 = this->find(cond_cl && !cond_cu);
      this->i9 = this->find(!cond_cl && !cond_cu && !cond_cq);

      // Set right-hand side values
      this->b2 = bl(this->i2);
      this->l3 = bl(this->i3);
      this->u4 = bu(this->i4);
      this->l5 = bl(this->i5);
      this->u5 = bu(this->i5);
      this->b6 = cl(this->i6);
      this->l7 = cl(this->i7);
      this->u8 = cu(this->i8);
      this->l9 = cl(this->i9);
      this->u9 = cu(this->i9);

      // Set sizes of indices sets
      this->n0 = x0.size();
      this->n1 = this->i1.count();
      this->n2 = this->i2.count();
      this->n3 = this->i3.count();
      this->n4 = this->i4.count();
      this->n5 = this->i5.count();
      this->n6 = this->i6.count();
      this->n7 = this->i7.count();
      this->n8 = this->i8.count();
      this->n9 = this->i9.count();

      // Initialize number of invalid bounds
      this->vi = 0;

      // Count invalid bounds
      if (this->n2 > 0) {
        this->vi += (this->b2.array() <= -params.rhs_bnd).count();
        this->vi += (this->b2.array() >=  params.rhs_bnd).count();
      }
      if (this->n3 > 0) {
        this->vi += (this->l3.array() >=  params.rhs_bnd).count();
      }
      if (this->n4 > 0) {
        this->vi += (this->u4.array() <= -params.rhs_bnd).count();
      }
      if (this->n5 > 0) {
        this->vi += (this->l5.array() >=  params.rhs_bnd).count();
        this->vi += (this->u5.array() <= -params.rhs_bnd).count();
        this->vi += (this->l5.array() >   this->u5.array()).count();
      }
      if (this->n6 > 0) {
        this->vi += (this->b6.array() <= -params.rhs_bnd).count();
        this->vi += (this->b6.array() >=  params.rhs_bnd).count();
      }
      if (this->n7 > 0) {
        this->vi += (this->l7.array() >=  params.rhs_bnd).count();
      }
      if (this->n8 > 0) {
        this->vi += (this->u8.array() <= -params.rhs_bnd).count();
      }
      if (this->n9 > 0) {
        this->vi += (this->l9.array() >=  params.rhs_bnd).count();
        this->vi += (this->u9.array() <= -params.rhs_bnd).count();
        this->vi += (this->l9.array() >   this->u9.array()).count();
      }

      // Set number of variables and constraints
      this->nV = this->n1 + this->n3 + this->n4 + this->n5;
      this->nI = this->n3 + this->n4 + 2*this->n5 + this->n7 + this->n8 + 2*this->n9;
      this->nE = this->n6;

      // Set size of primal-dual matrix
      this->nA = this->nV + 3*this->nE + 3*this->nI;

      // Set initial point
      this->x0.resize(this->nV);

      this->x0 << x0(this->i1), x0(this->i3), x0(this->i4), x0(this->i5);
      std::cout << "this->x0 = " << this->x0.transpose() << std::endl;

      #undef CMD
    }

  private:
    /**
     * \brief Select elements from a vector based on a boolean mask.
     * \param[in] vector The input vector.
     * \param[in] mask The boolean mask.
     * \return The selected elements from the input vector.
     */
    template<typename MaskType>
    static Indices find(MaskType const & mask)
    {
      Indices out(mask.count());
      for (Integer i{0}, j{0}; i < mask.size(); ++i) {
        if (mask[i]) {out[j++] = i;}
      }
      return out;
    }

  }; // struct Input

} // namespace Pipal

#endif /* INCLUDE_PIPAL_INPUT_HH */
