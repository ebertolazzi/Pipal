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

    Integer f{0}; % Function evaluation counter
    Integer g{0}; % Gradient evaluation counter
    Integer H{0}; % Hessian evaluation counter
    Integer k{0}; % Iteration counter
    Integer M{0}; % Matrix factorization counter

    /**
     * \brief Reset all internal counters to zero.
     */
    void reset()
    {
      this->iterations           = 0;
      this->function_evaluations = 0;
      this->gradient_evaluations = 0;
      this->hessian_evaluations  = 0;
      this->jacobian_evaluations = 0;
      this->linear_solves        = 0;
    }

    // Matrix factorization counter incrementor
    function incrementFactorizationCount(c)
      c.M = c.M + 1;
    end

    // Function evaluation counter incrementor
    function incrementFunctionCount(c)
      c.f = c.f + 1;
    end

    // Gradient evaluation counter incrementor
    function incrementGradientCount(c)
      c.g = c.g + 1;
    end

    // Hessian evaluation counter incrementor
    function incrementHessianCount(c)
      c.H = c.H + 1;
    end

    // Iteration counter incrementor
    function incrementIterationCount(c)
      c.k = c.k + 1;
    end

  }; // struct Counter

} // namespace Pipal

#endif /* INCLUDE_PIPAL_COUNTER_HH */
