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

    Real p0{0.0};  /*!< Fraction-to-the-boundary steplength. */
    Real parameter{0.0};   /*!< Primal steplength. */
    Real direction{0.0};   /*!< Dual steplength. */
    bool s{false}; /*!< Bool for second-order correction. */

    // Backtracking line search
    void backtracking(Parameter<Real> const & parameter, Input<Real> const & input,
      Iterate<Real> const & iterate, Direction<Real> const & direction)
    {
      // Store current values
      Vector x(iterate.x); Real f{iterate.f};
      cE = iterate.cE; r1 = iterate.r1; r2 = iterate.r2; lE = iterate.lE;
      cI = iterate.cI; s1 = iterate.s1; s2 = iterate.s2; lI = iterate.lI;
      phi = iterate.phi;

      // Backtracking loop
      while (this->parameter >= eps)
      {
        // Set trial point
        iterate.update_point  (input,  direction, *this);
        iterate.eval_functions(input,counter    );

        // Check for function evaluation error
        if (iterate.err == 0)
        {
          // Set remaining trial values
          iterate.eval_slacks(parameter, input);
          iterate.eval_merit(input);

          // Check for nonlinear fraction-to-boundary violation
          Integer ftb{0};
          if (input.nE > 0) {
            ftb += (iterate.r1 < std::min(parameter.ls_frac,iterate.mu)*r1).count() +
                   (iterate.r2 < std::min(parameter.ls_frac,iterate.mu)*r2).count();
          }
          if (input.nI > 0) {
            ftb += (iterate.s1 < std::min(parameter.ls_frac,iterate.mu)*s1).count() +
                   (iterate.s2 < std::min(parameter.ls_frac,iterate.mu)*s2).count();
          }

          // Check Armijo condition
          if (ftb == 0 && iterate.phi - phi <= -parameter.ls_thresh * this->parameter * std::max(direction.qtred, 0.0))
          {
            // Reset variables and return
            iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, iterate.f, iterate.cE, iterate.cI, iterate.phi);
            return;
          }
          else
          {
            // Reset variables
            iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
          }
        }
        else
        {
          // Reset variables
          iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }

        // Reduce steplength
        this->parameter *= parameter.ls_factor;
      }
    }

    // Fraction-to-boundary line search
    void fraction_to_boundary(Parameter<Real> const & parameter, Input<Real> const & input,
      Iterate<Real> const & iterate, Direction<Real> const & direction)
    {
      // Initialize primal fraction-to-boundary
      this->p0 = 1;

      // Update primal fraction-to-boundary for constraint slacks
      if (input.nE > 0) {
        this->p0 = std::min([
          this->p0;
          (min(parameter.ls_frac,iterate.mu)-1)*iterate.r1(direction.r1<0)./direction.r1(direction.r1<0);
          (std::min(parameter.ls_frac,iterate.mu) - 1) * iterate.r2(direction.r2 < 0) ./ direction.r2(direction.r2 < 0)
        ]);
      }
      if (input.nI > 0) {
        this->p0 = std::min([
          this->p0;
          (min(parameter.ls_frac,iterate.mu)-1)*iterate.s1(direction.s1<0)./direction.s1(direction.s1<0);
          (std::min(parameter.ls_frac,iterate.mu) - 1) * iterate.s2(direction.s2 < 0) ./ direction.s2(direction.s2 < 0)
        ]);
      }

      // Initialize primal steplength
      this->parameter = this->p0;

      // Initialize dual fraction-to-boundary
      this->direction = 1;

      // Update dual fraction-to-boundary for constraint multipliers
      if (input.nE > 0) {
        this->direction = std::min([
          this->direction;
          (std::min(parameter.ls_frac,iterate.mu)-1)*(1+iterate.lE(direction.lE < 0))./direction.lE(direction.lE < 0);
          (1-std::min(parameter.ls_frac,iterate.mu))*(1-iterate.lE(direction.lE > 0))./direction.lE(direction.lE > 0)
        ]);
      }
      if (input.nI > 0) {
        this->direction = std::min([
          this->direction;
          (std::min(parameter.ls_frac,iterate.mu)-1)*(0+iterate.lI(direction.lI < 0))./direction.lI(direction.lI < 0);
          (1-std::min(parameter.ls_frac,iterate.mu))*(1-iterate.lI(direction.lI > 0))./direction.lI(direction.lI > 0)
        ]);
      }
    }

    // Full step search for trial penalty parameters
    Integer full_step_check(Parameter<Real> const & parameter, Input<Real> const & input,
       Counter<Integer> & counter, Iterate<Real> const & iterate, Direction<Real> const & direction)
    {
      // Set current and last penalty parameters
      rho      = iterate.rho;
      rho_temp = iterate.rho_;

      // Loop through last penalty parameters
      while (rho < rho_temp)
      {
        // Set penalty parameter
        iterate.set_rho(rho_temp);

        // Evaluate merit
        iterate.eval_merit(input);

        // Store current values
        x = iterate.x; f = iterate.f;
        cE = iterate.cE; r1 = iterate.r1; r2 = iterate.r2; lE = iterate.lE;
        cI = iterate.cI; s1 = iterate.s1; s2 = iterate.s2; lI = iterate.lI;
        phi = iterate.phi;

        // Set trial point
        iterate.update_point(input, direction, *this);
        iterate.eval_functions(input, counter);

        // Check for function evaluation error
        if (iterate.err == 0)
        {
          // Set remaining trial values
          iterate.eval_slacks(parameter, input);
          iterate.eval_merit(input);

          // Check for nonlinear fraction-to-boundary violation
          Integer ftb{0};
          if (input.nE > 0) {
            ftb += (iterate.r1 < std::min(parameter.ls_frac, iterate.mu)*r1).count() +
                  (iterate.r2 < std::min(parameter.ls_frac,iterate.mu)*r2).count();
          }
          if (input.nI > 0) {
            ftb += (iterate.s1 < std::min(parameter.ls_frac, iterate.mu)*s1).count() +
                  (iterate.s2 < std::min(parameter.ls_frac,iterate.mu)*s2).count();
          }

          // Check Armijo condition
          if (ftb == 0 && iterate.phi - phi <= -parameter.ls_thresh*this->parameter * std::max(direction.qtred, 0.0))
          {
            // Reset variables, set boolean, and return
            iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, iterate.f, iterate.cE, iterate.cI, iterate.phi);
            return 1;
          }
          else
          {
            // Reset variables
            iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
          }
        }
        else
        {
          // Reset variables
          iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }

        // Decrease rho
        rho_temp = parameter.rho_factor*rho_temp;

        return 0;
      }

      // Set rho
      iterate.set_rho(rho);

      // Evaluate merit
      iterate.eval_merit(input);
    }

    // Line search
    void line_search(Parameter<Real> const & parameter, Input<Real> const & input,
       Counter<Integer> & counter, Iterate<Real> const & iterate, Direction<Real> const & direction)
    {
      // Check fraction-to-boundary rule
      this->fraction_to_boundary(parameter, input, iterate, direction);

      // Check for full step for trial penalty parameters
      Integer b{this->full_step_check(parameter, input, counter, iterate, direction)};

      // Run second-order correction
      this->s = 0;
      if (b == 0) {
        b = this->second_order_correction(parameter, input, counter, iterate, direction);
        if (b == 2) {this->s = 1}
      }

      // Run backtracking line search
      if (b == 0) {this->backtracking(parameter, input, counter, iterate, direction);}
    }

    // Second-order Correction
    Integer second_order_correction(Parameter<Real> const & parameter, Input<Real> const & input,
       Counter<Integer> & counter, Iterate<Real> const & iterate, Direction<Real> const & direction)
    {
      // Store current iterate values
      x = iterate.x; f = iterate.f;
      cE = iterate.cE; r1 = iterate.r1; r2 = iterate.r2; lE = iterate.lE;
      cI = iterate.cI; s1 = iterate.s1; s2 = iterate.s2; lI = iterate.lI;
      phi = iterate.phi; v = iterate.v;

      // Set trial point
      iterate.update_point(input, direction, *this);
      iterate.eval_functions(input, counter);

      // Check for function evaluation error
      if (iterate.err == 0)
      {
        // Set remaining trial values
        iterate.eval_slacks(parameter, input);
        iterate.eval_merit(input);

        // Check for nonlinear fraction-to-boundary violation
        Integer ftb{0};
        if (input.nE > 0) {
          ftb += (iterate.r1 < std::min(parameter.ls_frac,iterate.mu)*r1).count() +
                 (iterate.r2 < std::min(parameter.ls_frac,iterate.mu)*r2).count();
        }
        if (input.nI > 0) {
          ftb += (iterate.s1 < std::min(parameter.ls_frac,iterate.mu)*s1).count() +
                 (iterate.s2 < std::min(parameter.ls_frac,iterate.mu)*s2).count();
        }

        // Check Armijo condition
        if (ftb == 0 && iterate.phi - phi <= -parameter.ls_thresh*this->parameter * std::max(direction.qtred, 0.0))
        {
          // Reset variables, set flag, and return
          iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, iterate.f, iterate.cE, iterate.cI, iterate.phi);
          return 1;
        }
        else if (iterate.eval_violation(input,iterate.cE,iterate.cI) < v)
        {
          // Reset variables and return
          iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
          return 0;
        }
        else
        {
          // Reset variables (but leave constraint values for second-order correction)
          iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, f,iterate.cE,iterate.cI,phi);
        }
      }
      else
      {
        // Reset variables and return
        iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        return 0;
      }

      // Recompute slacks for second order correction
      iterate.eval_slacks(parameter, input);

      // Evaluate trial primal-dual right-hand side vector
      iterate.eval_newtonr_rhs(input);

      // Store current direction values
      dx  = direction.x ;
      dr1 = direction.r1; dr2 = direction.r2; dlE = direction.lE;
      ds1 = direction.s1; ds2 = direction.s2; dlI = direction.lI;
      dx_norm = direction.x_norm;
      dl_norm = direction.l_norm;

      // Evaluate search direction
      direction.eval_newton_step(input,iterate);

      // Set trial direction
      direction.set_direction(input, this->parameter*dx+direction.x, this->parameter*dr1+direction.r1,
        this->parameter*dr2+direction.r2, this->parameter*ds1+direction.s1, this->parameter*ds2+direction.s2,
        this->direction*dlE+direction.lE, this->direction*dlI+direction.lI,
        norm(this->parameter*dx+direction.x),norm([this->direction*dlE+direction.lE;this->direction*dlI+direction.lI]));

      // Set trial point
      iterate.update_point(input,direction,a);
      iterate.eval_functions(input,counter);

      // Check for function evaluation error
      if (iterate.err == 0)
      {
        // Set remaining trial values
        iterate.eval_slacks(parameter, input);
        iterate.eval_merit (input);

        // Check for nonlinear fraction-to-boundary violation
        Integer ftb{0};
        if (input.nE > 0) {
          ftb += (iterate.r1 < std::min(parameter.ls_frac, iterate.mu)*r1).count() +
                (iterate.r2 < std::min(parameter.ls_frac,iterate.mu)*r2).count();
        }
        if (input.nI > 0) {
          ftb += (iterate.s1 < std::min(parameter.ls_frac, iterate.mu)*s1).count() +
                (iterate.s2 < std::min(parameter.ls_frac,iterate.mu)*s2).count();
        }

        // Check Armijo condition
        if (ftb == 0 && iterate.phi - phi <= -parameter.ls_thresh*this->parameter*std::max(direction.qtred,0))
        {
          // Reset variables, set flag, and return
          iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, iterate.f, iterate.cE, iterate.cI, iterate.phi);
          return 2;
        }
        else
        {
          // Reset variables
          iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
        }
      }
      else
      {
        // Reset variables
        iterate.set_primals(input, x, r1, r2, s1, s2, lE, lI, f, cE, cI, phi);
      }

      // Reset direction
      direction.set_direction(input,dx,dr1,dr2,ds1,ds2,dlE,dlI,dx_norm,dl_norm);

      // Reduce steplength
      this->parameter *= parameter.ls_factor;

      return 0;
    }

  }; // struct Acceptance

} // namespace Pipal

#endif /* INCLUDE_PIPAL_ACCEPTANCE_HH */
