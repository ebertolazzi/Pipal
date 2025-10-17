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

#ifndef INCLUDE_PIPAL_DEFINES_HH
#define INCLUDE_PIPAL_DEFINES_HH

// Standard libraries
#include <iostream>
#include <sstream>
#include <stdexcept>

// Print Pipal errors
#ifndef PIPAL_ERROR
#define PIPAL_ERROR(MSG)                \
  {                                     \
    std::ostringstream os;              \
    os << MSG;                          \
    throw std::runtime_error(os.str()); \
  }
#endif

// Assert for Pipal
#ifndef PIPAL_ASSERT
#define PIPAL_ASSERT(COND, MSG) \
  if (!(COND))                  \
  {                             \
    PIPAL_ERROR(MSG);           \
  }
#endif

// Warning for Pipal
#ifndef PIPAL_WARNING
#define PIPAL_WARNING(MSG)         \
  {                                \
    std::cout << MSG << std::endl; \
  }
#endif

// Warning assert for Pipal
#ifndef PIPAL_ASSERT_WARNING
#define PIPAL_ASSERT_WARNING(COND, MSG) \
  if (!(COND))                          \
  {                                     \
    PIPAL_WARNING(MSG);                 \
  }
#endif

#ifndef PIPAL_DEFAULT_INTEGER_TYPE
#define PIPAL_DEFAULT_INTEGER_TYPE int
#endif

/**
 * \brief Namespace for the Pipal library.
 *
 * Penalty-interior-point algorithm (Pipal) is a library for nonlinear constrained optimization with
 * inequality constraints (it does not explicitly handle equality constraints). Precisely speaking,
 * it will compute the solution to the optimization problem
 * \f[
 *  \begin{array}{l}
 *    \text{minimize} ~ f(\mathbf{x}) \\
 *    \text{subject to} ~ \mathbf{c}(\mathbf{x}) \leq \mathbf{0}
 *  \end{array} \text{,}
 * \f]
 * where \f$\mathbf{x} \in \mathbb{R}^n\f$ is the vector of optimization variables, \f$f: \mathbb{R}^n
 * \to \mathbb{R}\f$ is the objective function, and \f$\mathbf{c}: \mathbb{R}^n \to \mathbb{R}^m\f$
 * are the constraints.
 *
 * \note To create an equality constraint \f$ h(\mathbf{x}) = 0 \f$, one can define two inequality
 * constraints \f$ h(\mathbf{x}) \leq 0 \f$ and \f$ -h(\mathbf{x}) \leq 0 \f$.
 *
 * This code is mostly based on the descriptions provided in this reference:
 *
 * - Frank E. Curtis. A penalty-interior-point algorithm for nonlinear constrained optimization.
 *   Mathematical Programming Computation (2012) 4:181-209. DOI: 10.1007/s12532-012-0041-4.
 */
namespace Pipal
{
  /**
   * \brief The Integer type as used for the API.
   *
   * The Integer type, \c \#define the preprocessor symbol \c PIPAL_DEFAULT_INTEGER_TYPE. The default
   * value is \c int.
   */
  using Integer = PIPAL_DEFAULT_INTEGER_TYPE;

  /**
   * \brief Select elements from a vector based on a boolean mask.
   * \tparam Indices The type of the output indices.
   * \tparam Mask The type of the boolean mask.
   * \param[in] vector The input vector.
   * \param[in] mask The boolean mask.
   * \return The selected elements from the input vector.
   */
  template<typename Indices, typename Mask>
  static Indices find(Mask const & mask)
  {
    Indices out(mask.count());
    for (Integer i{0}, j{0}; i < mask.size(); ++i) {
      if (mask[i]) {out[j++] = i;}
    }
    return out;
  }

  /**
   * \brief Enumeration for matrix view types.
   */
  using MatrixView = enum class MatrixView : Integer {
    FULL = 0, /*!< Full matrix view. */
    TRIL = 1, /*!< Lower triangular matrix view. */
    TRIU = 2  /*!< Upper triangular matrix view. */
  };

  /**
   * \brief Count the number of non-zero elements in a matrix.
   * \tparam Real The floating-point type.
   * \tparam Matrix The matrix type.
   * \tparam View The matrix view type (FULL, TRIL, TRIU).
   * \param[in] mat The input matrix.
   * \param[in] tol The tolerance for considering an element as non-zero.
   * \return The number of non-zero elements in the matrix.
   */
  template<typename Real, typename Matrix, MatrixView View = MatrixView::FULL>
  static Integer nnz(Matrix const & mat, Integer Real const tol = Eigen::NumTraits<Real>::epsilon())
  {
    Integer count{0};
    const Integer rows{static_cast<Integer>(mat.rows())};
    const Integer cols{static_cast<Integer>(mat.cols())};
    for (Integer j{0}; j < cols; ++j) {
      const Integer i_start{(View == MatrixView::TRIU) ? j : 0};
      const Integer i_end{(View == MatrixView::TRIL) ? std::min(j + 1, rows) : rows};
      for (Integer i{i_start}; i < i_end; ++i) {
        if (std::abs(mat(i, j)) > tol) {++count;}
      }
    }
    return count;

  }


} // namespace Pipal

#endif // INCLUDE_PIPAL_DEFINES_HH
