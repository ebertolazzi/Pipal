/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (counter) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Pipal project is distributed under the MIT License.                                       *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_PIPAL_ACCEPTANCE_HH
#define INCLUDE_PIPAL_ACCEPTANCE_HH

// Pipal includes
#include "Pipal/Types.hh"
#include "Pipal/Iterate.hh"
#include "Pipal/Input.hh"
#include "Pipal/Parameter.hh"
#include "Pipal/Direction.hh"
#include "Pipal/Counter.hh"

namespace Pipal
{

  // Backtracking line search
  void backtracking(Acceptance & a, Parameter & p, Input & i, Counter & c,
    Iterate & z, Direction & d)
  {
    // Store current values
    Vector x(z.x);
    Real f{z.f};
    Vector cE(z.cE), r1(z.r1), r2(z.r2), lE(z.lE);
    Vector cI(z.cI), s1(z.s1), s2(z.s2), lI(z.lI);
    Real phi{z.phi};

    // Backtracking loop
    while (a.p >= EPSILON)
    {
      // Set trial point
      updatePoint(z, i, d, a);
      evalFunctions(z, i, c);

      // Check for function evaluation error
      if (z.err == 0)
      {
        // Set remaining trial values
        evalSlacks(z, p, i);
        evalMerit(z, i);

        // Check for nonlinear fraction-to-boundary violation
        Integer ftb{0};
        if (i.nE > 0) {
          ftb += (z.r1.array() < (std::min(p.ls_frac, z.mu)*r1).array()).count() + (z.r2.array() < (std::min(p.ls_frac, z.mu)*r2).array()).count();
        }
        if (i.nI > 0) {
          ftb += (z.s1.array() < (std::min(p.ls_frac, z.mu)*s1).array()).count() + (z.s2.array() < (std::min(p.ls_frac, z.mu)*s2).array()).count();
        }

        // Check Armijo condition
        if (ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*std::max(d.qtred, 0.0))
        {
          // Reset variables and return
          setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
          return;
        }
        else
        {
          // Reset variables
          setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }
      }
      else
      {
        // Reset variables
        setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      }

      // Reduce steplength
      a.p *= p.ls_factor;
    }
  }
  // Fraction-to-boundary line search
  void fractionToBoundary(Acceptance & a, Parameter & p, Input & i, Iterate & z,
    Direction & d)
  {

    // Initialize primal fraction-to-boundary
    a.p0 = 1.0;

    // Update primal fraction-to-boundary for constraint slacks
    Real frac{std::min(p.ls_frac, z.mu)};
    if (i.nE > 0) {
      Indices const idx_r1(find(d.r1.array() < 0.0));
      Indices const idx_r2(find(d.r1.array() < 0.0));
      Real min_1{INFINITY}, min_2{INFINITY};
      if (idx_r1.size() > 0) {min_1 = (((frac - 1.0) * z.r1(idx_r1).array()) / d.r1(idx_r1).array()).minCoeff();}
      if (idx_r2.size() > 0) {min_2 = (((frac - 1.0) * z.r2(idx_r2).array()) / d.r2(idx_r2).array()).minCoeff();}
      a.p0 = std::min({a.p0, min_1, min_2});
    }
    if (i.nI > 0) {
      Indices const idx_s1(find(d.s1.array() < 0.0));
      Indices const idx_s2(find(d.s1.array() < 0.0));
      Real min_1{INFINITY}, min_2{INFINITY};
      if (idx_s1.size() > 0) {min_1 = (((frac - 1.0) * z.s1(idx_s1).array()) / d.s1(idx_s1).array()).minCoeff();}
      if (idx_s2.size() > 0) {min_2 = (((frac - 1.0) * z.s2(idx_s2).array()) / d.s2(idx_s2).array()).minCoeff();}
      a.p0 = std::min({a.p0, min_1, min_2});
    }

    // Initialize primal steplength
    a.p = a.p0;

    // Initialize dual fraction-to-boundary
    a.d = 1.0;

    // Update dual fraction-to-boundary for constraint multipliers
    if (i.nE > 0) {
      Indices const idx_l(find(d.lE.array() < 0.0));
      Indices const idx_g(find(d.lE.array() > 0.0));
      Real min_1{INFINITY}, min_2{INFINITY};
      if (idx_l.size() > 0) {min_1 = (((frac - 1.0) * (1.0 + z.lE(idx_l).array())) / d.lE(idx_l).array()).minCoeff();}
      if (idx_g.size() > 0) {min_2 = (((1.0 - frac) * (1.0 - z.lE(idx_g).array())) / d.lE(idx_g).array()).minCoeff();}
      a.d = std::min({a.d, min_1, min_2});
    }
    if (i.nI > 0) {
      Indices const idx_l(find(d.lI.array() < 0.0));
      Indices const idx_g(find(d.lI.array() > 0.0));
      Real min_1{INFINITY}, min_2{INFINITY};
      if (idx_l.size() > 0) {min_1 = (((frac - 1.0) * (0.0 + z.lI(idx_l).array())) / d.lI(idx_l).array()).minCoeff();}
      if (idx_g.size() > 0) {min_2 = (((1.0 - frac) * (1.0 - z.lI(idx_g).array())) / d.lI(idx_g).array()).minCoeff();}
      a.d = std::min({a.d, min_1, min_2});
    }
  }

  // Full step search for trial penalty parameters
  Integer fullStepCheck(Acceptance & a, Parameter & p, Input & i, Counter & c, Iterate & z, Direction & d)
  {
    // Initialize boolean
    Integer b{0};

    // Set current and last penalty parameters
    Real rho{z.rho}, rho_temp{z.rho_};

    // Loop through last penalty parameters
    while (rho < rho_temp)
    {
      // Set penalty parameter
      setRho(z, rho_temp);

      // Evaluate merit
      evalMerit(z, i);

      // Store current values
      Vector x(z.x);
      Real f{z.f};
      Vector cE(z.cE), r1(z.r1), r2(z.r2), lE(z.lE);
      Vector cI(z.cI), s1(z.s1), s2(z.s2), lI(z.lI);
      Real phi{z.phi};

      // Set trial point
      updatePoint(z, i, d, a);
      evalFunctions(z, i, c);

      // Check for function evaluation error
      if (z.err == 0)
      {
        // Set remaining trial values
        evalSlacks(z, p, i);
        evalMerit(z, i);

        // Check for nonlinear fraction-to-boundary violation
        Integer ftb{0};
        if (i.nE > 0) {
          ftb += (z.r1.array() < (std::min(p.ls_frac, z.mu)*r1).array()).count() + (z.r2.array() < (std::min(p.ls_frac, z.mu)*r2).array()).count();
        }
        if (i.nI > 0) {
          ftb += (z.s1.array() < (std::min(p.ls_frac, z.mu)*s1).array()).count() + (z.s2.array() < (std::min(p.ls_frac, z.mu)*s2).array()).count();
        }

        // Check Armijo condition
        if (ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*std::max(d.qtred, 0.0))
        {
          // Reset variables, set boolean, and return
          setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi); b = 1;
          return b;
        }
        else
        {
          // Reset variables
          setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }
      }
      else
      {
        // Reset variables
        setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      }

      // Decrease rho
      rho_temp *= p.rho_factor;
    }

    // Set rho
    setRho(z, rho);

    // Evaluate merit
    evalMerit(z, i);

    return b;
  }

  // Line search
  void lineSearch(Acceptance & a, Parameter & p, Input & i, Counter & c,
    Iterate & z, Direction & d)
  {
    // Check fraction-to-boundary rule
    fractionToBoundary(a, p, i, z, d);

    // Check for full step for trial penalty parameters
    Integer b{fullStepCheck(a, p, i, c, z, d)};
    // Run second-order correction
    a.s = 0;
    if (b == 0) {
      b = secondOrderCorrection(a, p, i, c, z, d);
      if (b == 2) {a.s = 1;}
    }

    // Run backtracking line search
    if (b == 0) {backtracking(a, p, i, c, z, d);}
  }

  // Second-order Correction
  Integer secondOrderCorrection(Acceptance & a, Parameter & p, Input & i, Counter & c,
    Iterate & z, Direction & d)
  {
    // Initialize flag
    Integer b{0};

    // Store current iterate values
    Vector x(z.x);
    Real f{z.f};
    Vector cE(z.cE), r1(z.r1), r2(z.r2), lE(z.lE);
    Vector cI(z.cI), s1(z.s1), s2(z.s2), lI(z.lI);
    Real phi{z.phi}, v{z.v};

    // Set trial point
    updatePoint(z, i, d, a);
    evalFunctions(z, i, c);

    // Check for function evaluation error
    if (z.err == 0)
    {
      // Set remaining trial values
      evalSlacks(z, p, i);
      evalMerit(z, i);

      // Check for nonlinear fraction-to-boundary violation
      Integer ftb{0};
      if (i.nE > 0) {
        ftb += (z.r1.array() < std::min(p.ls_frac, z.mu)*r1.array()).count() + (z.r2.array() < std::min(p.ls_frac, z.mu)*r2.array()).count();
      }
      if (i.nI > 0) {
        ftb += (z.s1.array() < std::min(p.ls_frac, z.mu)*s1.array()).count() + (z.s2.array() < std::min(p.ls_frac, z.mu)*s2.array()).count();
      }

      // Check Armijo condition
      if (ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*std::max(d.qtred, 0.0))
      {
        // Reset variables, set flag, and return
        setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
        b = 1; return b;
      }
      else if (evalViolation(i, z.cE, z.cI) < v)
      {
        // Reset variables and return
        setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        return b;
      }
      else
      {
        // Reset variables (but leave constraint values for second-order correction)
        setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, f, z.cE, z.cI, phi);
      }
    }
    else
    {
      // Reset variables and return
      setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      return b;
    }

    // Recompute slacks for second order correction
    evalSlacks(z, p, i);

    // Evaluate trial primal-dual right-hand side vector
    evalNewtonRhs(z, i);

    // Store current direction values
    Vector dx(d.x);
    Vector dr1(d.r1), dr2(d.r2), dlE(d.lE);
    Vector ds1(d.s1), ds2(d.s2), dlI(d.lI);
    Real dx_norm{d.x_norm};
    Real dl_norm{d.l_norm};

    // Evaluate search direction
    evalNewtonStep(d, i, z);

    // Set trial direction
    setDirection(d, i, a.p*dx+d.x, a.p*dr1+d.r1, a.p*dr2+d.r2, a.p*ds1+d.s1,
      a.p*ds2+d.s2, a.d*dlE+d.lE, a.d*dlI+d.lI, (a.p*dx+d.x).norm(),
      std::sqrt((a.d*dlE+d.lE).squaredNorm() + (a.d*dlI+d.lI).squaredNorm()));

    // Set trial point
    updatePoint(z, i, d, a);
    evalFunctions(z, i, c);

    // Check for function evaluation error
    if (z.err == 0)
    {
      // Set remaining trial values
      evalSlacks(z, p, i);
      evalMerit(z, i);

      // Check for nonlinear fraction-to-boundary violation
      Integer ftb{0};
      if (i.nE > 0) {
        ftb += (z.r1.array() < std::min(p.ls_frac, z.mu)*r1.array()).count() + (z.r2.array() < std::min(p.ls_frac, z.mu)*r2.array()).count();
      }
      if (i.nI > 0) {
        ftb += (z.s1.array() < std::min(p.ls_frac, z.mu)*s1.array()).count() + (z.s2.array() < std::min(p.ls_frac, z.mu)*s2.array()).count();
      }

      // Check Armijo condition
      if (ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*std::max(d.qtred, 0.0))
      {
        // Reset variables, set flag, and return
        setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
        b = 2; return b;
      }
      else
      {
        // Reset variables
        setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      }
    }
    else
    {
      // Reset variables
      setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
    }

    // Reset direction
    setDirection(d, i, dx, dr1, dr2, ds1, ds2, dlE, dlI, dx_norm, dl_norm);

    // Reduce steplength
    a.p *= p.ls_factor;

    return b;
  }

} // namespace Pipal

#endif /* INCLUDE_PIPAL_ACCEPTANCE_HH */
