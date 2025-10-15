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

    Integer iterations{0};           /*!< Number of iterations. */
    Integer function_evaluations{0}; /*!< Number of function evaluations. */
    Integer gradient_evaluations{0}; /*!< Number of gradient evaluations. */
    Integer hessian_evaluations{0};  /*!< Number of Hessian evaluations. */
    Integer jacobian_evaluations{0}; /*!< Number of Jacobian evaluations. */
    Integer linear_solves{0};        /*!< Number of linear system solves. */

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

  }; // struct Counter

} // namespace Pipal

#endif /* INCLUDE_PIPAL_COUNTER_HH */
