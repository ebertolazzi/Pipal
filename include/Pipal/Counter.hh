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

#ifndef INCLUDE_PIPAL_COUNTER_HH
#define INCLUDE_PIPAL_COUNTER_HH

// Pipal includes
#include "Pipal/Types.hh"

namespace Pipal
{

  /**
   * \brief Reset all internal counters to zero.
   * \param[out] c Counter object to reset.
   */
  inline void resetCounter(Counter & c) {c.f = c.g = c.H = c.k = c.M = 0;}

  /**
   * \brief Increment the matrix factorization counter.
   * \param[in] c Counter whose matrix-factorization count is incremented.
   */
  inline void incrementFactorizationCount(Counter & c) {++c.M;}

  /**
   * \brief Increment the function evaluation counter.
   * \param[in] c Counter whose function-evaluation count is incremented.
   */
  inline void incrementFunctionCount(Counter & c) {++c.f;}

  /**
   * \brief Increment the gradient evaluation counter.
   * \param[in] c Counter whose gradient-evaluation count is incremented.
   */
  inline void incrementGradientCount(Counter & c) {++c.g;}

  /**
   * \brief Increment the Hessian evaluation counter.
   * \param[in] c Counter whose Hessian-evaluation count is incremented.
   */
  inline void incrementHessianCount(Counter & c) {++c.H;}

  /**
   * \brief Increment the iteration counter.
   * \param[in] c Counter whose iteration count is incremented.
   */
  inline void incrementIterationCount(Counter & c) {++c.k;}

} // namespace Pipal

#endif /* INCLUDE_PIPAL_COUNTER_HH */
