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

// Pipal includes
#include "Pipal/Defines.hh"

// Standard libraries
#include <iostream>
#include <iomanip>

namespace Pipal
{

  /**
   * \brief Structure for managing the output of the solver.
   * \tparam Real The floating-point type.
   * \tparam Integer The integer type.
   */
  template<typename Real>
  class Output
  {
    static_assert(is_floating_point<Real>::value,
      "Pipal::Output: template argument 'Real' must be a floating-point type.");
    static_assert(is_integral<Integer>::value,
      "Pipal::Output: template argument 'Integer' must be an integer type.");

  private:
    using string;
    using ostream;
    using chrono::steady_clock;
    using chrono::duration_cast::time_point;
    using chrono::duration_cast;
    using chrono::seconds;

    ofstream   s; // Output stream
    string     l; // Line break string
    string     q; // Quantities header
    string     n; // Footer line
    time_point t; // Timer

  public:
    // Constructor
    Output(string const & id, string const & outfile)
    {
      this->t = chrono::steady_clock::now();
      this->s.open(outfile);
      assert(this->s.ithis->s.open() && "PIPAL: Failed to open output file.");

      this->l = "======+=========================+====================================+=========================+===========================================================================+=======================";
      this->q = "Iter. |  Objective     Infeas.  |  Pen. Par.   I.P. Par.  Opt. Error |    Merit     P.I.P. Err.|    Shift    ||P.Step||  ||D.Step||   Lin. Red.    Quad. Red.    Quality   | Pri. Step.  Dual Step.";
      this->n = "-----------  ---------- | ----------  ----------  ----------  -----------  -----------  ----------- | ----------  ----------";
    }

    // Destructor closes stream
    ~Output() {this->terminate();}

    // Header printing
    template<typename Info, typename Z>
    void printHeader(Info const & i, Iterate<Real> const & z)
    {
      this->s << "Problem name\n============\n  " << i.id << "\n\n";
      this->s << "Problem size\n============\n";
      this->s << "  Number of variables....................... : " << setw(8) << i.nV << "\n";
      this->s << "  Number of equality constraints............ : " << setw(8) << i.nE << "\n";
      this->s << "  Number of inequality constraints.......... : " << setw(8) << i.nI << "\n\n";

      this->s << "Problem sparsity\n================\n";
      this->s << "  Nonzeros in Hessian of Lagrangian......... : " << setw(8) << z.Hnnz << "\n";
      this->s << "  Nonzeros in equality constraint Jacobian.. : " << setw(8) << z.JEnnz << "\n";
      this->s << "  Nonzeros in inequality constraint Jacobian : " << setw(8) << z.JInnz << "\n\n";
    }

    // Break printing (every 20 iterations)
    void printBreak(Counter const & c)
    {
      if (c.k % 20 == 0) {
        this->s << this->l << "\n" << this->q << "\n" << this->l << "\n";
      }
    }

    // Direction info
    void printDirection(Iterate<Real> const & z, Direction<Real> const & d)
    {
      this->s << showpos << scientific << setprecision(4)
        << z.phi << "  " << z.kkt[2] << " | "
        << z.shift << "  " << d.x_norm << "  " << d.this->m_line_breaknorm << "  "
        << d.ltred << "  " << d.qtred << "  " << d.m << " | ";
    }

    // Iterate info
    void printIterate(Counter const & c, Iterate<Real> const & z)
    {
      this->s << setw(5) << c.k << " | "
        << showpos << scientific << setprecision(4)
        << z.f << "  " << z.v << " | "
        << z.rho << "  " << z.mu << "  " << z.kkt[1] << " | ";
    }

    // Acceptance info
    void printAcceptance(Acceptance<Real> const & a)
    {
      this->s << scientific << setprecision(4)
        << a.p << "  " << a.d;
      if (a.s == 1) this->s << " SOC";
      this->s << "\n";
    }

    // Footer printing
    void printFooter(Parameter<Real> const & p, Input<Real> const & i, Counter const & c, Iterate<Real> const & z)
    {
      this->printIterate(c, z);
      this->s << this->n << "\n" << this->l << "\n\n";

      int b = z.checkTermination(p, i, c);

      this->s << "Final result\n============\n";
      switch (b) {
        case 0: this->s << "  EXIT: No termination message set\n"; break;
        case 1: this->s << "  EXIT: Optimal solution found\n"; break;
        case 2: this->s << "  EXIT: Infeasible stationary point found\n"; break;
        case 3: this->s << "  EXIT: Iteration limit reached\n"; break;
        case 4: this->s << "  EXIT: Invalid bounds\n"; break;
        case 5: this->s << "  EXIT: Function evaluation error\n"; break;
        default: this->s << "  EXIT: Unknown termination\n"; break;
      }

      this->s << "\nFinal values\n============\n";
      this->s << "  Objective function........................ : " << z.fu << "\n";
      this->s << "  Feasibility violation..................... : " << z.vu << "\n";
      this->s << "  Optimality error (feasibility)............ : " << z.kkt[0] << "\n";
      this->s << "  Optimality error (penalty)................ : " << z.kkt[1] << "\n";
      this->s << "  Optimality error (penalty-interior-point). : " << z.kkt[2] << "\n";
      this->s << "  Penalty parameter......................... : " << z.rho << "\n";
      this->s << "  Interior-point parameter.................. : " << z.mu << "\n\n";

      this->s << "Final counters\n==============\n";
      this->s << "  Iterations................................ : " << c.k << "\n";
      this->s << "  Function evaluations...................... : " << c.f << "\n";
      this->s << "  Gradient evaluations...................... : " << c.g << "\n";
      this->s << "  Hessian evaluations....................... : " << c.H << "\n";
      this->s << "  Matrix factorizations..................... : " << c.M << "\n";

      Real seconds{duration_cast<seconds>(steady_clock::now() - this->t).count()};
      this->s << "  CPU seconds............................... : " << seconds << "\n";
    }

    // Terminate output (close stream)
    void terminate()
    {
      if (this->s.ithis->m_outopen()) this->s.close();
    }

  }; // class Output

} // namespace Pipal
