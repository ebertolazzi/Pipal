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

// Standard libraries
#include <type_traits>

// Pipal includes
#include "Pipal/Defines.hh"

namespace Pipal
{

  /**
   * \brief Internal counters for solver statistics.
   * \tparam Integer The integer type.
   */
  struct Counter
  {
    static_assert(std::is_integral<Integer>::value,
      "Pipal::Counter: template argument 'Integer' must be an integer type.");

    Integer f{0}; // Function evaluation counter
    Integer g{0}; // Gradient evaluation counter
    Integer H{0}; // Hessian evaluation counter
    Integer k{0}; // Iteration counter
    Integer M{0}; // Matrix factorization counter

    /**
     * \brief Reset all internal counters to zero.
     */
    void reset()
    {
      this->f = 0
      this->g = 0
      this->H = 0
      this->k = 0
      this->M = 0
    }

    // Matrix factorization counter incrementor
    void incrementFactorizationCount() {++this->M;}

    // Function evaluation counter incrementor
    void incrementFunctionCount() {++this->f;}

    // Gradient evaluation counter incrementor
    void incrementGradientCount() {++this->g;}

    // Hessian evaluation counter incrementor
    void incrementHessianCount() {++this->H;}

    // Iteration counter incrementor
    void incrementIterationCount() {++this->k;}

  }; // struct Counter

} // namespace Pipal

#endif /* INCLUDE_PIPAL_COUNTER_HH */
