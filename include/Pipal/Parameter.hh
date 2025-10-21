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

namespace Pipal
{

  // Constructor
  inline void buildParameter(Parameter & p, Algorithm a) {p.algorithm = a;}

  // Reset interior-point parameter maximum exponent in increases to default
  inline void resetMuMaxExp(Parameter & p) {p.mu_max_exp = p.mu_max_exp0;}

  // Set interior-point parameter maximum exponent in increases to zero
  inline void setMuMaxExpZero(Parameter & p) {p.mu_max_exp = 0.0;}

} // namespace Pipal

#endif /* INCLUDE_PIPAL_PARAMETER_HH */
