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
#include "Pipal/Types.hh"

namespace Pipal {
  /**
   * \brief Initialize algorithm parameters.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[out] p Parameter object to initialize.
   * \param[in] a Algorithm selection enumerator.
   */
  template <typename Real>
  inline
  void
  buildParameter(Parameter<Real> & p, Algorithm a) {p.algorithm = a;}

  /**
   * \brief Reset maximum exponent used for mu increases to its default.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] p Parameter object whose mu exponent limit is reset.
   */
  template <typename Real>
  inline
  void
  resetMuMaxExp(Parameter<Real> & p) {p.mu_max_exp = p.mu_max_exp0;}

  /**
   * \brief Force mu exponent increases to use zero as maximum exponent.
   * \tparam Real Floating-point type used by the algorithm.
   * \param[in] p Parameter object to modify.
   */
  template <typename Real>
  inline
  void
  setMuMaxExpZero(Parameter<Real> & p) {p.mu_max_exp = 0.0;}

} // namespace Pipal

#endif /* INCLUDE_PIPAL_PARAMETER_HH */
