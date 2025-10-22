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
#include <iostream>
#include <iomanip>
#include <string>
#include <ostream>

// Pipal includes
#include "Pipal/Types.hh"

namespace Pipal
{

  class Output
  {
    using Ostream      = std::ostream;
    using MicroSeconds = std::chrono::microseconds;
    using SteadyClock  = std::chrono::steady_clock;
    using TimePoint    = SteadyClock::time_point;

    Ostream &   s{std::cout}; // Output stream (reference to avoid copying std::cout)
    std::string l; // Line break std::string
    std::string q; // Quantities header
    std::string n; // Footer line
    TimePoint   t; // Timer

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
    ~Output() = default;

    // Header printing
    void printHeader(Input const & i, Iterate const & z) const {
      this->s
        << "Problem name" << std::endl
        << "============" << std::endl
        << "  " << i.id << std::endl
        << std::endl;

      this->s
        << "Problem size" << std::endl << "============" << std::endl
        << "  Number of variables....................... : " << i.nV << std::endl
        << "  Number of equality constraints............ : " << i.nE << std::endl
        << "  Number of inequality constraints.......... : " << i.nI << std::endl
        << std::endl;

      this->s
        << "Problem sparsity" << std::endl << "================" << std::endl
        << "  Nonzeros in Hessian of Lagrangian......... : " << z.Hnnz << std::endl
        << "  Nonzeros in equality constraint Jacobian.. : " << z.JEnnz << std::endl
        << "  Nonzeros in inequality constraint Jacobian : " << z.JInnz << std::endl
        << std::endl;
    }

    // Break printing (every 20 iterations)
    void printBreak(Counter const & c) const {
      if (c.k % 20 == 0) {
        this->s << this->l << std::endl << this->q << std::endl << this->l << std::endl;
      }
    }

    // Direction info
    void printDirection(Iterate const & z, Direction const & d) const {
      this->s
        << std::scientific << std::setprecision(4)
        << std::showpos << z.phi << "  " << std::noshowpos << z.kkt[2] << " | " << z.shift << "  " << d.x_norm << "  " << d.l_norm << "  "
        << std::showpos << d.ltred << "  " << d.qtred << "  " << d.m << std::noshowpos << " | ";
    }

    // Iterate info
    void printIterate(Counter const & c, Iterate const & z) const {
      this->s << std::setw(5) << c.k << " | " << std::scientific << std::setprecision(4)
        << std::showpos << z.f << std::noshowpos << "  " << z.v << " | " << z.rho << "  " << z.mu << "  " << z.kkt[1] << " | ";
    }

    // Acceptance info
    void printAcceptance(Acceptance const & a) const {
      this->s << std::scientific << std::setprecision(4) << a.p << "  " << a.d;
        if (a.s == 1) {this->s << " SOC";}
      this->s << std::endl;
    }

    // Footer printing
    void printFooter(Parameter const & p, Input const & i, Counter const & c, Iterate const & z) const {
      this->printIterate(c, z);
      this->s
        << this->n << std::endl << this->l << std::endl
        << std::endl;

      const Integer b{checkTermination(z, p, i, c)};

      this->s
        << "Final result" << std::endl
        << "============" << std::endl;
      switch (b) {
        case 0:  this->s << "  EXIT: No termination message set" << std::endl; break;
        case 1:  this->s << "  EXIT: Optimal solution found" << std::endl; break;
        case 2:  this->s << "  EXIT: Infeasible stationary point found" << std::endl; break;
        case 3:  this->s << "  EXIT: Iteration limit reached" << std::endl; break;
        case 4:  this->s << "  EXIT: Invalid bounds" << std::endl; break;
        case 5:  this->s << "  EXIT: Function evaluation error" << std::endl; break;
        default: this->s << "  EXIT: Unknown termination" << std::endl; break;
      }

      this->s
        << "\nFinal values" << std::endl
        << "============" << std::showpos << std::endl
        << "  Objective function........................ : " << z.fu << std::endl
        << "  Feasibility violation..................... : " << z.vu << std::endl
        << "  Optimality error (feasibility)............ : " << z.kkt(0) << std::endl
        << "  Optimality error (penalty)................ : " << z.kkt(1) << std::endl
        << "  Optimality error (penalty-interior-point). : " << z.kkt(2) << std::endl
        << "  Penalty parameter......................... : " << z.rho << std::endl
        << "  Interior-point parameter.................. : " << z.mu << std::noshowpos << std::endl
        << std::endl;

      this->s
        << "Final counters" << std::endl
        << "==============" << std::endl
        << "  Iterations................................ : " << c.k << std::endl
        << "  Function evaluations...................... : " << c.f << std::endl
        << "  Gradient evaluations...................... : " << c.g << std::endl
        << "  Hessian evaluations....................... : " << c.H << std::endl
        << "  Matrix factorizations..................... : " << c.M << std::endl
        << "  CPU millseconds........................... : " << std::scientific << std::setprecision(4)
        << std::chrono::duration_cast<MicroSeconds>(SteadyClock::now() - this->t).count()/1.0e3 << std::endl;
    }

  }; // class Output

} // namespace Pipal

#endif // INCLUDE_PIPAL_OUTPUT_HH
