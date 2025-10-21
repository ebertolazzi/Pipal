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

  // Reset all internal counters to zero
  inline void buildCounter(Counter & c) {c.f = c.g = c.H = c.k = c.M = 0;}

  // Matrix factorization counter incrementor
  inline void incrementFactorizationCount(Counter & c) {++c.M;}

  // Function evaluation counter incrementor
  inline void incrementFunctionCount(Counter & c) {++c.f;}

  // Gradient evaluation counter incrementor
  inline void incrementGradientCount(Counter & c) {++c.g;}

  // Hessian evaluation counter incrementor
  inline void incrementHessianCount(Counter & c) {++c.H;}

  // Iteration counter incrementor
  inline void incrementIterationCount(Counter & c) {++c.k;}

} // namespace Pipal

#endif /* INCLUDE_PIPAL_COUNTER_HH */
