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
  void buildCounter(Counter & c) {c.f = c.g = c.H = c.k = c.M = 0;}

  // Matrix factorization counter incrementor
  void incrementFactorizationCount(Counter & c) {++c.M;}

  // Function evaluation counter incrementor
  void incrementFunctionCount(Counter & c) {++c.f;}

  // Gradient evaluation counter incrementor
  void incrementGradientCount(Counter & c) {++c.g;}

  // Hessian evaluation counter incrementor
  void incrementHessianCount(Counter & c) {++c.H;}

  // Iteration counter incrementor
  void incrementIterationCount(Counter & c) {++c.k;}

} // namespace Pipal

#endif /* INCLUDE_PIPAL_COUNTER_HH */
