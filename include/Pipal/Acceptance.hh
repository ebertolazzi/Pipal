/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (counter) 2025, Davide Stocco and Enrico Bertolazzi.                                *
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
#include "Pipal/Direction.hh"

namespace Pipal
{

  // Backtracking line search
  inline void backtracking(Acceptance & a, Parameter & p, Input & i, Counter & c, Iterate & z, Direction const & d)
  {
    // Store current values
    Vector x(z.x);
    Real f{z.f};
    Array cE(z.cE), r1(z.r1), r2(z.r2), lE(z.lE);
    Array cI(z.cI), s1(z.s1), s2(z.s2), lI(z.lI);
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
          ftb += (z.r1 < (std::min(p.ls_frac, z.mu)*r1)).count() + (z.r2 < (std::min(p.ls_frac, z.mu)*r2)).count();
        }
        if (i.nI > 0) {
          ftb += (z.s1 < (std::min(p.ls_frac, z.mu)*s1)).count() + (z.s2 < (std::min(p.ls_frac, z.mu)*s2)).count();
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

      // Reduce step length
      a.p *= p.ls_factor;
    }
  }
  // Fraction-to-boundary line search
  inline void fractionToBoundary(Acceptance & a, Parameter & p, Input const & i, Iterate & z, Direction & d)
  {

    // Initialize primal fraction-to-boundary
    a.p0 = 1.0;

    // Update primal fraction-to-boundary for constraint slacks
    Real frac{std::min(p.ls_frac, z.mu)};
    if (i.nE > 0) {
      Indices const idx_r1(find(d.r1 < 0.0));
      Indices const idx_r2(find(d.r2 < 0.0));
      Real min_1{INFINITY}, min_2{INFINITY};
      if (idx_r1.size() > 0) {min_1 = (((frac - 1.0) * z.r1(idx_r1)) / d.r1(idx_r1)).minCoeff();}
      if (idx_r2.size() > 0) {min_2 = (((frac - 1.0) * z.r2(idx_r2)) / d.r2(idx_r2)).minCoeff();}
      a.p0 = std::min(a.p0, std::min(min_1, min_2));
    }
    if (i.nI > 0) {
      Indices const idx_s1(find(d.s1 < 0.0));
      Indices const idx_s2(find(d.s2 < 0.0));
      Real min_1{INFINITY}, min_2{INFINITY};
      if (idx_s1.size() > 0) {min_1 = (((frac - 1.0) * z.s1(idx_s1)) / d.s1(idx_s1)).minCoeff();}
      if (idx_s2.size() > 0) {min_2 = (((frac - 1.0) * z.s2(idx_s2)) / d.s2(idx_s2)).minCoeff();}
      a.p0 = std::min(a.p0, std::min(min_1, min_2));
    }

    // Initialize primal step length
    a.p = a.p0;

    // Initialize dual fraction-to-boundary
    a.d = 1.0;

    // Update dual fraction-to-boundary for constraint multipliers
    if (i.nE > 0) {
      Indices const idx_l(find(d.lE < 0.0));
      Indices const idx_g(find(d.lE > 0.0));
      Real min_1{INFINITY}, min_2{INFINITY};
      if (idx_l.size() > 0) {min_1 = (((frac - 1.0) * (1.0 + z.lE(idx_l))) / d.lE(idx_l)).minCoeff();}
      if (idx_g.size() > 0) {min_2 = (((1.0 - frac) * (1.0 - z.lE(idx_g))) / d.lE(idx_g)).minCoeff();}
      a.d = std::min({a.d, min_1, min_2});
    }
    if (i.nI > 0) {
      Indices const idx_l(find(d.lI < 0.0));
      Indices const idx_g(find(d.lI > 0.0));
      Real min_1{INFINITY}, min_2{INFINITY};
      if (idx_l.size() > 0) {min_1 = (((frac - 1.0) * (0.0 + z.lI(idx_l))) / d.lI(idx_l)).minCoeff();}
      if (idx_g.size() > 0) {min_2 = (((1.0 - frac) * (1.0 - z.lI(idx_g))) / d.lI(idx_g)).minCoeff();}
      a.d = std::min({a.d, min_1, min_2});
    }
  }

  // Full step search for trial penalty parameters
  inline Integer fullStepCheck(Acceptance const & a, Parameter & p, Input & i, Counter & c, Iterate & z, Direction const & d)
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
      Array cE(z.cE), r1(z.r1), r2(z.r2), lE(z.lE);
      Array cI(z.cI), s1(z.s1), s2(z.s2), lI(z.lI);
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
          ftb += (z.r1 < (std::min(p.ls_frac, z.mu)*r1)).count() + (z.r2 < (std::min(p.ls_frac, z.mu)*r2)).count();
        }
        if (i.nI > 0) {
          ftb += (z.s1 < (std::min(p.ls_frac, z.mu)*s1)).count() + (z.s2 < (std::min(p.ls_frac, z.mu)*s2)).count();
        }

        // Check Armijo condition
        if (ftb == 0 && z.phi - phi <= -p.ls_thresh*a.p*std::max(d.qtred, 0.0))
        {
          // Reset variables, set boolean, and return
          setPrimals(z, i, x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
          b = 1; return b;
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
  inline void lineSearch(Acceptance & a, Parameter & p, Input & i, Counter & c,
    Iterate & z, Direction & d)
  {
    // Check fraction-to-boundary rule
    fractionToBoundary(a, p, i, z, d);

    // Check for full step for trial penalty parameters
    Integer b{fullStepCheck(a, p, i, c, z, d)};
    // Run second-order correction
    a.s = false;
    if (b == 0) {
      b = secondOrderCorrection(a, p, i, c, z, d);
      if (b == 2) {a.s = true;}
    }

    // Run backtracking line search
    if (b == 0) {backtracking(a, p, i, c, z, d);}
  }

  // Second-order Correction
  inline Integer secondOrderCorrection(Acceptance & a, Parameter & p, Input & i, Counter & c, Iterate & z, Direction & d)
  {
    // Initialize flag
    Integer b{0};

    // Store current iterate values
    Vector x(z.x);
    Real f{z.f};
    Array cE(z.cE), r1(z.r1), r2(z.r2), lE(z.lE);
    Array cI(z.cI), s1(z.s1), s2(z.s2), lI(z.lI);
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
        ftb += (z.r1 < std::min(p.ls_frac, z.mu)*r1).count() + (z.r2 < std::min(p.ls_frac, z.mu)*r2).count();
      }
      if (i.nI > 0) {
        ftb += (z.s1 < std::min(p.ls_frac, z.mu)*s1).count() + (z.s2 < std::min(p.ls_frac, z.mu)*s2).count();
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
    setDirection(d, i, a.p*dx+d.x, a.p*dr1+d.r1.matrix(), a.p*dr2+d.r2.matrix(), a.p*ds1+d.s1.matrix(),
      a.p*ds2+d.s2.matrix(), a.d*dlE+d.lE.matrix(), a.d*dlI+d.lI.matrix(), (a.p*dx+d.x).norm(),
      std::sqrt((a.d*dlE+d.lE.matrix()).matrix().squaredNorm() + (a.d*dlI+d.lI.matrix()).squaredNorm()));

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
        ftb += (z.r1 < std::min(p.ls_frac, z.mu)*r1).count() + (z.r2 < std::min(p.ls_frac, z.mu)*r2).count();
      }
      if (i.nI > 0) {
        ftb += (z.s1 < std::min(p.ls_frac, z.mu)*s1).count() + (z.s2 < std::min(p.ls_frac, z.mu)*s2).count();
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

    // Reduce step length
    a.p *= p.ls_factor;

    return b;
  }

} // namespace Pipal

#endif /* INCLUDE_PIPAL_ACCEPTANCE_HH */
