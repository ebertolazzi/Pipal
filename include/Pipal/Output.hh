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

#ifndef INCLUDE_PIPAL_OUTPUT_HH
#define INCLUDE_PIPAL_OUTPUT_HH

// Standard libraries
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
#include <ostream>

// Pipal includes
#include "Pipal/Types.hh"
#include "Pipal/Iterate.hh"

namespace Pipal
{

  class Output
  {
    static_assert(std::is_floating_point<Real>::value,
      "Pipal::Output: template argument 'Real' must be a floating-point type.");

  private:
    using String      = std::string;
    using Ostream     = std::ostream;
    using Seconds     = std::chrono::seconds;
    using SteadyClock = std::chrono::steady_clock;
    using TimePoint   = SteadyClock::time_point;

    Ostream & s{std::cout}; // Output stream (reference to avoid copying std::cout)
    String    l; // Line break string
    String    q; // Quantities header
    String    n; // Footer line
    TimePoint t; // Timer

  public:
    // Constructor
    Output()
    {
      this->t = SteadyClock::now();
      this->l = "======+=========================+====================================+=========================+===========================================================================+=======================";
      this->q = "Iter. |  Objective     Infeas.  |  Pen. Par.   I.P. Par.  Opt. Error |    Merit     P.I.P. Err.|    Shift    ||P.Step||  ||D.Step||   Lin. Red.    Quad. Red.    Quality   | Pri. Step.  Dual Step.";
      this->n = "-----------  ---------- | ----------  ----------  ----------  -----------  -----------  ----------- | ----------  ----------";
    }

    // Constructor
    // FIXME Output(Ostream & stream) : Output() {this->s = &stream;}

    // Destructor closes stream
    ~Output() {}

    // Header printing
    void printHeader(Input & i, Iterate & z)
    {
      this->s << "Problem name" << std::endl << "============" << std::endl << "  " << i.id << std::endl << std::endl;

      this->s << "Problem size" << std::endl << "============" << std::endl;
      this->s << "  Number of variables....................... : " << std::setw(8) << i.nV << std::endl;
      this->s << "  Number of equality constraints............ : " << std::setw(8) << i.nE << std::endl;
      this->s << "  Number of inequality constraints.......... : " << std::setw(8) << i.nI << std::endl << std::endl;

      this->s << "Problem sparsity" << std::endl << "================" << std::endl;
      this->s << "  Nonzeros in Hessian of Lagrangian......... : " << std::setw(8) << z.Hnnz << std::endl;
      this->s << "  Nonzeros in equality constraint Jacobian.. : " << std::setw(8) << z.JEnnz << std::endl;
      this->s << "  Nonzeros in inequality constraint Jacobian : " << std::setw(8) << z.JInnz << std::endl << std::endl;
    }

    // Break printing (every 20 iterations)
    void printBreak(Counter & c)
    {
      if (c.k % 20 == 0) {
        this->s << this->l << std::endl << this->q << std::endl << this->l << std::endl;
      }
    }

    // Direction info
    void printDirection(Iterate & z, Direction & d)
    {
      this->s << std::showpos << std::scientific << std::setprecision(4)
        << z.phi << "  " << z.kkt[2] << " | " << z.shift << "  " << d.x_norm << "  " << d.l_norm << "  "
        << d.ltred << "  " << d.qtred << "  " << d.m << " | ";
    }

    // Iterate info
    void printIterate(Counter & c, Iterate & z)
    {
      this->s << std::setw(5) << c.k << " | " << std::showpos << std::scientific << std::setprecision(4)
        << z.f << "  " << z.v << " | " << z.rho << "  " << z.mu << "  " << z.kkt[1] << " | ";
    }

    // Acceptance info
    void printAcceptance(Acceptance & a)
    {
      this->s << std::scientific << std::setprecision(4) << a.p << "  " << a.d;
        if (a.s == 1) {this->s << " SOC";}
      this->s << std::endl;
    }

    // Footer printing
    void printFooter(Parameter & p, Input & i, Counter & c, Iterate & z)
    {
      this->printIterate(c, z);
      this->s << this->n << std::endl << this->l << std::endl << std::endl;

      const Integer b{checkTermination(z, p, i, c)};

      this->s << "Final result" << std::endl << "============" << std::endl;
      switch (b) {
        case 0:  this->s << "  EXIT: No termination message set" << std::endl; break;
        case 1:  this->s << "  EXIT: Optimal solution found" << std::endl; break;
        case 2:  this->s << "  EXIT: Infeasible stationary point found" << std::endl; break;
        case 3:  this->s << "  EXIT: Iteration limit reached" << std::endl; break;
        case 4:  this->s << "  EXIT: Invalid bounds" << std::endl; break;
        case 5:  this->s << "  EXIT: Function evaluation error" << std::endl; break;
        default: this->s << "  EXIT: Unknown termination" << std::endl; break;
      }

      this->s << "\nFinal values" << std::endl << "============" << std::endl;
      this->s << "  Objective function........................ : " << z.fu << std::endl;
      this->s << "  Feasibility violation..................... : " << z.vu << std::endl;
      this->s << "  Optimality error (feasibility)............ : " << z.kkt(0) << std::endl;
      this->s << "  Optimality error (penalty)................ : " << z.kkt(1) << std::endl;
      this->s << "  Optimality error (penalty-interior-point). : " << z.kkt(2) << std::endl;
      this->s << "  Penalty parameter......................... : " << z.rho << std::endl;
      this->s << "  Interior-point parameter.................. : " << z.mu << std::endl << std::endl;

      this->s << "Final counters" << std::endl << "==============" << std::endl;
      this->s << "  Iterations................................ : " << c.k << std::endl;
      this->s << "  Function evaluations...................... : " << c.f << std::endl;
      this->s << "  Gradient evaluations...................... : " << c.g << std::endl;
      this->s << "  Hessian evaluations....................... : " << c.H << std::endl;
      this->s << "  Matrix factorizations..................... : " << c.M << std::endl;
      this->s << "  CPU seconds............................... : " <<
          std::chrono::duration_cast<Seconds>(SteadyClock::now() - this->t).count() << std::endl;
    }

  }; // class Output

} // namespace Pipal

#endif // INCLUDE_PIPAL_OUTPUT_HH
