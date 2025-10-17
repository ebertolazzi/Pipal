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
#include "Pipal/Defines.hh"

namespace Pipal
{

  struct Acceptance
  {

    Real p0{0.0};  // Fraction-to-the-boundary steplength
    Real p{0.0};   // Primal steplength
    Real d{0.0};   // Dual steplength
    bool s{false}; // Bool for second-order correction

    // Backtracking line search
    void backtracking(Parameter<Real> const & p, Input<Real> const & i, Counter const & c,
      Iterate<Real> const & z, Direction<Real> const & d)
    {
      // Store current values
      x = z.x; f = z.f;
      cE = z.cE; r1 = z.r1; r2 = z.r2; lE = z.lE;
      cI = z.cI; s1 = z.s1; s2 = z.s2; lI = z.lI;
      phi = z.phi;

      // Backtracking loop
      while (this->p >= eps)
      {
        // Set trial point
        z.updatePoint(i, d, a);
        z.evalFunctions(i, c);

        // Check for function evaluation error
        if (z.err == 0)
        {
          // Set remaining trial values
          z.evalSlacks(p, i);
          z.evalMerit(i);

          // Check for nonlinear fraction-to-boundary violation
          Integer ftb{0};
          if (i.nE > 0) {
            ftb += (z.r1 < std::min(p.ls_frac, z.mu)*r1).count() + (z.r2 < std::min(p.ls_frac, z.mu)*r2).count();
          }
          if (i.nI > 0) {
            ftb += (z.s1 < std::min(p.ls_frac, z.mu)*s1).count() + (z.s2 < std::min(p.ls_frac, z.mu)*s2).count();
          }

          // Check Armijo condition
          if (ftb == 0 && z.phi - phi <= -p.ls_thresh*this->p*std::max(d.qtred, 0))
          {
            // Reset variables and return
            z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
            return;
          }
          else
          {
            // Reset variables
            z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
          }
        }
        else
        {
          // Reset variables
          z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }

        // Reduce steplength
        this->p *= p.ls_factor;
      }
    }

    // Fraction-to-boundary line search
    void fractionToBoundary(Parameter<Real> const & p, Input<Real> const & i, Iterate<Real> const & z,
      Direction<Real> const & d)
    {

      // Initialize primal fraction-to-boundary
      this->p0 = 1;

      // Update primal fraction-to-boundary for constraint slacks
      if (i.nE > 0) {
        this->p0 = std::min({this->p0,
          (((std::min(p.ls_frac, z.mu)-1)*z.r1(d.r1 < 0)).array() / (d.r1(d.r1 < 0)).array()).minCoeff(),
          (((std::min(p.ls_frac, z.mu)-1)*z.r2(d.r2 < 0)).array() / (d.r2(d.r2 < 0)).array()).minCoeff()
        });
      }
      if (i.nI > 0) {
        this->p0 = std::min({this->p0,
          (((std::min(p.ls_frac, z.mu)-1)*z.s1(d.s1 < 0)).array() / (d.s1(d.s1 < 0)).array()).minCoeff(),
          (((std::min(p.ls_frac, z.mu)-1)*z.s2(d.s2 < 0)).array() / (d.s2(d.s2 < 0)).array()).minCoeff()
        });
      }

      // Initialize primal steplength
      this->p = this->p0;

      // Initialize dual fraction-to-boundary
      this->d = 1;

      // Update dual fraction-to-boundary for constraint multipliers
      if (i.nE > 0) {
        this->d = std::min({this->d,
          (((std::min(p.ls_frac, z.mu)-1)*(1+z.lE(d.lE < 0))).array() / (d.lE(d.lE < 0)).array()).minCoeff(),
          (((1-std::min(p.ls_frac, z.mu))*(1-z.lE(d.lE > 0))).array() / (d.lE(d.lE > 0)).array()).minCoeff()
        });
      } // FIXME
      if (i.nI > 0) {
        this->d = std::min({this->d;
          (((std::min(p.ls_frac, z.mu)-1)*(0+z.lI(d.lI < 0))).array() / (d.lI(d.lI < 0)).array()).minCoeff(),
          (((1-std::min(p.ls_frac, z.mu))*(1-z.lI(d.lI > 0))).array() / (d.lI(d.lI > 0)).array()).minCoeff()
        });
      } // FIXME
    }

    // Full step search for trial penalty parameters
    Integer fullStepCheck(Parameter<Real> const & p, Input<Real> const & i, Counter const & c,
      Iterate<Real> const & z, Direction<Real> const & d)
    {
      // Initialize boolean
      b = 0;

      // Set current and last penalty parameters
      rho      = z.rho;
      rho_temp = z.rho_;

      // Loop through last penalty parameters
      while (rho < rho_temp)
      {
        // Set penalty parameter
        z.setRho(rho_temp);

        // Evaluate merit
        z.evalMerit(i);

        // Store current values
        x = z.x; f = z.f;
        cE = z.cE; r1 = z.r1; r2 = z.r2; lE = z.lE;
        cI = z.cI; s1 = z.s1; s2 = z.s2; lI = z.lI;
        phi = z.phi;

        // Set trial point
        z.updatePoint(i, d, a);
        z.evalFunctions(i, c);

        // Check for function evaluation error
        if (z.err == 0)
        {
          // Set remaining trial values
          z.evalSlacks(p, i);
          z.evalMerit(i);

          // Check for nonlinear fraction-to-boundary violation
          Integer ftb{0};
          if (i.nE > 0) {
            ftb += (z.r1 < std::min(p.ls_frac, z.mu)*r1).count() + (z.r2 < std::min(p.ls_frac, z.mu)*r2).count();
          }
          if (i.nI > 0) {
            ftb += (z.s1 < std::min(p.ls_frac, z.mu)*s1).count() + (z.s2 < std::min(p.ls_frac, z.mu)*s2).count();
          }

          // Check Armijo condition
          if (ftb == 0 && z.phi - phi <= -p.ls_thresh*this->p*std::max(d.qtred, 0))
          {
            // Reset variables, set boolean, and return
            z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi); b = 1;
            return b;
          }
          else
          {
            // Reset variables
            z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
          }
        }
        else
        {
          // Reset variables
          z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }

        // Decrease rho
        rho_temp = p.rho_factor*rho_temp;

        return b;
      }

      // Set rho
      z.setRho(rho);

      // Evaluate merit
      z.evalMerit(i);

      return b;
    }

    // Line search
    void lineSearch(Parameter<Real> const & p, Input<Real> const & i, Counter const & c,
      Iterate<Real> & z, Direction<Real> & d)
    {
      // Check fraction-to-boundary rule
      this->fractionToBoundary(p, i, z, d);

      // Check for full step for trial penalty parameters
      b = this->fullStepCheck(p, i, c, z, d);

      // Run second-order correction
      this->s = 0;
      if (b == 0) {
        b = this->secondOrderCorrection(p, i, c, z, d);
        if (b == 2) {this->s = 1;}
      }

      // Run backtracking line search
      if (b == 0) {this->backtracking(p, i, c, z, d);}
    }

    // Second-order Correction
    void secondOrderCorrection(Parameter<Real> const & p, Input<Real> const & i, Counter const & c,
      Iterate<Real> & z, Direction<Real> & d)
    {
      // Initialize flag
      b = 0;

      // Store current iterate values
      x = z.x; f = z.f;
      cE = z.cE; r1 = z.r1; r2 = z.r2; lE = z.lE;
      cI = z.cI; s1 = z.s1; s2 = z.s2; lI = z.lI;
      phi = z.phi; v = z.v;

      // Set trial point
      z.updatePoint(i, d, a);
      z.evalFunctions(i, c);

      // Check for function evaluation error
      if (z.err == 0)
      {
        // Set remaining trial values
        z.evalSlacks(p, i);
        z.evalMerit(i);

        // Check for nonlinear fraction-to-boundary violation
        Integer ftb{0};
        if (i.nE > 0) {
          ftb += (z.r1 < std::min(p.ls_frac, z.mu)*r1).count() + (z.r2 < std::min(p.ls_frac, z.mu)*r2).count();
        }
        if (i.nI > 0) {
          ftb += (z.s1 < std::min(p.ls_frac, z.mu)*s1).count() + (z.s2 < std::min(p.ls_frac, z.mu)*s2).count();
        }

        // Check Armijo condition
        if (ftb == 0 && z.phi - phi <= -p.ls_thresh*this->p*std::max(d.qtred, 0))
        {
          // Reset variables, set flag, and return
          z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
          b = 1; return;
        }
        else if (z.evalViolation(i,z.cE,z.cI) < v)
        {
          // Reset variables and return
          z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
          return;
        }
        else
        {
          // Reset variables (but leave constraint values for second-order correction)
          z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, f, z.cE, z.cI, phi);
        }
      }
      else
      {
        // Reset variables and return
        z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        return;
      }

      // Recompute slacks for second order correction
      z.evalSlacks(p, i);

      // Evaluate trial primal-dual right-hand side vector
      z.evalNewtonRhs(i);

      // Store current direction values
      dx  = d.x ;
      dr1 = d.r1; dr2 = d.r2; dlE = d.lE;
      ds1 = d.s1; ds2 = d.s2; dlI = d.lI;
      dx_norm = d.x_norm;
      dl_norm = d.l_norm;

      // Evaluate search direction
      d.evalNewtonStep(i, z);

      // Set trial direction
      d.setDirection(i, this->p*dx+d.x, this->p*dr1+d.r1, this->p*dr2+d.r2, this->p*ds1+d.s1,
        this->p*ds2+d.s2, this->d*dlE+d.lE, this->d*dlI+d.lI, this->p*dx+d.x.norm(),
        Vector(this->d*dlE+d.lE, this->d*dlI+d.lI).norm()); // OPTIMIZE

      // Set trial point
      z.updatePoint(i, d, a);
      z.evalFunctions(i, c);

      // Check for function evaluation error
      if (z.err == 0)
      {
        // Set remaining trial values
        z.evalSlacks(p, i);
        z.evalMerit(i);

        // Check for nonlinear fraction-to-boundary violation
        Integer ftb{0};
        if (i.nE > 0) {
          ftb += (z.r1 < std::min(p.ls_frac, z.mu)*r1).count() + (z.r2 < std::min(p.ls_frac, z.mu)*r2).count();
        }
        if (i.nI > 0) {
          ftb += (z.s1 < std::min(p.ls_frac, z.mu)*s1).count() + (z.s2 < std::min(p.ls_frac, z.mu)*s2).count();
        }

        // Check Armijo condition
        if (ftb == 0 && z.phi - phi <= -p.ls_thresh*this->p*max(d.qtred,0))
        {
          // Reset variables, set flag, and return
          z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, z.f, z.cE, z.cI, z.phi);
          b = 2; return;
        {
        else
        }
          // Reset variables
          z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }
      }
      else
      {
        // Reset variables
        z.setPrimals(i, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      }

      // Reset direction
      d.setDirection(i, dx, dr1, dr2, ds1, ds2, dlE, dlI, dx_norm, dl_norm);

      // Reduce steplength
      this->p *= p.ls_factor;
    }

  }; // struct Acceptance

} // namespace Pipal

#endif /* INCLUDE_PIPAL_ACCEPTANCE_HH */
