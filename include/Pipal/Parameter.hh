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

#ifndef INCLUDE_PIPAL_PARAMETER_HH
#define INCLUDE_PIPAL_PARAMETER_HH

// Pipal includes
#include "Pipal/Defines.hh"

// Standard libraries
#include <type_traits>

namespace Pipal
{

  /**
   * \brief Internal parameters for the solver algorithm.
   * \tparam Real The floating-point type.
   */
  template<typename Real>
  struct Parameter
  {
    static_assert(std::is_floating_point<Real>::value,
      "Pipal::Parameter: template argument 'Real' must be a floating-point type.");

    static constexpr Real    rhs_bnd{1.0e+18};       /*!< Maximum absolute value allowed for constraint right-hand side. */
    static constexpr Real    grad_max{1.0e+02};      /*!< Gradient norm limit for scaling. */
    static constexpr Real    infeas_max{1.0e+02};    /*!< Infeasibility limit for penalty parameter update. */
    static constexpr Real    nnz_max{2.0e+04};       /*!< Maximum non-zeros in (upper triangle of) Newton matrix. */
    static constexpr Integer opt_err_mem{6};         /*!< Optimality error history length. */
    static constexpr Real    ls_factor{5.0e-01};     /*!< Line search reduction factor. */
    static constexpr Real    ls_thresh{1.0e-08};     /*!< Line search threshold value. */
    static constexpr Real    ls_frac{1.0e-02};       /*!< Line search fraction-to-boundary constant. */
    static constexpr Real    pivot_thresh{5.0e-01};  /*!< Pivot threshold for LDL factorization. */
    static constexpr Real    slack_min{1.0e-20};     /*!< Slack variable bound. */
    static constexpr Real    shift_min{1.0e-12};     /*!< Hessian shift (non-zero) minimum value. */
    static constexpr Real    shift_factor1{5.0e-01}; /*!< Hessian shift update value (for decreases). */
    static constexpr Real    shift_factor2{6.0e-01}; /*!< Hessian shift update value (for increases). */
    static constexpr Real    shift_max{1.0e+08};     /*!< Hessian shift maximum value. */
    static constexpr Real    rho_init{1.0e-01};      /*!< Penalty parameter initial value. */
    static constexpr Real    rho_min{1.0e-12};       /*!< Penalty parameter minimum value. */
    static constexpr Real    rho_factor{5.0e-01};    /*!< Penalty parameter reduction factor. */
    static constexpr Integer rho_trials{8};          /*!< Penalty parameter number of trial values per iteration. */
    static constexpr Real    mu_init{1.0e-01};       /*!< Interior-point parameter initial value. */
    static constexpr Real    mu_min{1.0e-12};        /*!< Interior-point parameter minimum value. */
    static constexpr Real    mu_factor{1.0e-01};     /*!< Interior-point parameter reduction factor. */
    static constexpr Real    mu_factor_exp{1.5};     /*!< Interior-point parameter reduction exponent. */
    static constexpr Integer mu_trials{4};           /*!< Interior-point parameter number of trial values per iteration. */
    static constexpr Real    mu_max{1.0e-01};        /*!< Interior-point parameter maximum value. */
    static constexpr Real    mu_max_exp_0{0.0};      /*!< Interior-point parameter maximum exponent in increases (default). */
    static constexpr Real    update_con_1{1.0e-02};  /*!< Steering rule constant 1. */
    static constexpr Real    update_con_2{1.0e-02};  /*!< Steering rule constant 2. */
    static constexpr Real    update_con_3{1.01};     /*!< Adaptive interior-point rule constant. */

    Real m_mu_max_exp{0.0}; /*!< Interior-point parameter maximum exponent in increases. */

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

  }; // struct Parameter

} // namespace Pipal

#endif /* INCLUDE_PIPAL_PARAMETER_HH */
