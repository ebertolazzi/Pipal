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

#ifndef INCLUDE_PIPAL_INPUT_HH
#define INCLUDE_PIPAL_INPUT_HH

// Pipal includes
#include "Pipal/Types.hh"
#include "Pipal/Parameter.hh"

namespace Pipal
{

  // Constructor
  void buildInput(struct Input & i, struct Parameter & p, String const & name, ObjectiveFunc const & f_orig,
    ConstraintsFunc const & c_orig, ObjectiveGradientFunc const & g_orig,
    ConstraintsJacobianFunc const & J_orig, LagrangianHessianFunc const & H_orig,
    Vector const & x0, Vector const & bl, Vector const & bu, Vector const & cl, Vector const & cu)
  {
    #define CMD "Pipal::resetInput(...): "

    // Set problem identity
    i.id = name;

    // store function pointers to original problem functions
    i.f_fun = f_orig;
    i.c_fun = c_orig;
    i.g_fun = g_orig;
    i.J_fun = J_orig;
    i.H_fun = H_orig;

    // Set number of original formulation variables
    i.n0 = x0.size();

    // Find indices sets
    Real const tolerance{Eigen::NumTraits<Real>::epsilon()};
    Mask const cond_bl(bl.array() <= -p.rhs_bnd);
    Mask const cond_bu(bu.array() >=  p.rhs_bnd);
    Mask const cond_cl(cl.array() <= -p.rhs_bnd);
    Mask const cond_cu(cu.array() >=  p.rhs_bnd);
    Mask const cond_bq((bl.array() - bu.array()).abs() > tolerance);
    Mask const cond_cq((cl.array() - cu.array()).abs() > tolerance);

    i.I1 = find(cond_bl && cond_bu);
    i.I2 = find(cond_bq);
    i.I3 = find(!cond_bl &&  cond_bu);
    i.I4 = find(cond_bl && !cond_bu);
    i.I5 = find(!cond_bl && !cond_bu && !cond_bq);
    i.I6 = find(cond_cq);
    i.I7 = find(!cond_cl &&  cond_cu);
    i.I8 = find(cond_cl && !cond_cu);
    i.I9 = find(!cond_cl && !cond_cu && !cond_cq);

    // Set right-hand side values
    i.b2 = bl(i.I2);
    i.l3 = bl(i.I3);
    i.u4 = bu(i.I4);
    i.l5 = bl(i.I5);
    i.u5 = bu(i.I5);
    i.b6 = cl(i.I6);
    i.l7 = cl(i.I7);
    i.u8 = cu(i.I8);
    i.l9 = cl(i.I9);
    i.u9 = cu(i.I9);

    // Set sizes of indices sets
    i.n1 = i.I1.count();
    i.n2 = i.I2.count();
    i.n3 = i.I3.count();
    i.n4 = i.I4.count();
    i.n5 = i.I5.count();
    i.n6 = i.I6.count();
    i.n7 = i.I7.count();
    i.n8 = i.I8.count();
    i.n9 = i.I9.count();

    // Initialize number of invalid bounds
    i.vi = 0;

    // Count invalid bounds
    if (i.n2 > 0) {
      i.vi += (i.b2.array() <= -p.rhs_bnd).count();
      i.vi += (i.b2.array() >= p.rhs_bnd).count();
    }
    if (i.n3 > 0) {
      i.vi += (i.l3.array() >= p.rhs_bnd).count();
    }
    if (i.n4 > 0) {
      i.vi += (i.u4.array() <= -p.rhs_bnd).count();
    }
    if (i.n5 > 0) {
      i.vi += (i.l5.array() >= p.rhs_bnd).count();
      i.vi += (i.u5.array() <= -p.rhs_bnd).count();
      i.vi += (i.l5.array() > i.u5.array()).count();
    }
    if (i.n6 > 0) {
      i.vi += (i.b6.array() <= -p.rhs_bnd).count();
      i.vi += (i.b6.array() >= p.rhs_bnd).count();
    }
    if (i.n7 > 0) {
      i.vi += (i.l7.array() >= p.rhs_bnd).count();
    }
    if (i.n8 > 0) {
      i.vi += (i.u8.array() <= -p.rhs_bnd).count();
    }
    if (i.n9 > 0) {
      i.vi += (i.l9.array() >= p.rhs_bnd).count();
      i.vi += (i.u9.array() <= -p.rhs_bnd).count();
      i.vi += (i.l9.array() > i.u9.array()).count();
    }

    // Set number of variables and constraints
    i.nV = i.n1 + i.n3 + i.n4 + i.n5;
    i.nI = i.n3 + i.n4 + 2*i.n5 + i.n7 + i.n8 + 2*i.n9;
    i.nE = i.n6;

    // Set size of primal-dual matrix
    i.nA = i.nV + 3*i.nE + 3*i.nI;

    // Set initial point
    i.x0.resize(i.nV);
    i.x0 << x0(i.I1), x0(i.I3), x0(i.I4), x0(i.I5);

    #undef CMD
  }

} // namespace Pipal

#endif /* INCLUDE_PIPAL_INPUT_HH */
