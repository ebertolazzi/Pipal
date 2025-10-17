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

    ofstream   m_out;      // Output stream
    string     m_line_break;        // Line break string
    string     m_quantities;        // Quantities header
    string     m_footer;        // Footer line
    time_point m_start; // Timer

  public:
    // Constructor
    Output(string const & id, string const & outfile)
    {
      this->m_start = chrono::steady_clock::now();
      this->m_out.open(outfile);
      assert(this->m_out.ithis->m_out.open() && "PIPAL: Failed to open output file.");

      this->m_line_break = "======+=========================+====================================+=========================+===========================================================================+=======================";
      this->m_quantities = "Iter. |  Objective     Infeas.  |  Pen. Par.   I.P. Par.  Opt. Error |    Merit     P.I.P. Err.|    Shift    ||P.Step||  ||D.Step||   Lin. Red.    Quad. Red.    Quality   | Pri. Step.  Dual Step.";
      this->m_footer     = "-----------  ---------- | ----------  ----------  ----------  -----------  -----------  ----------- | ----------  ----------";
    }

    // Destructor closes stream
    ~Output() {this->terminate();}

    // ------------------------------------------------------------------------
    // Header printing
    // ------------------------------------------------------------------------
    template<typename Info, typename Z>
    void printHeader(Info const & i, Z const & z)
    {
      this->m_out << "Problem name\n============\n  " << i.id << "\n\n";
      this->m_out << "Problem size\n============\n";
      this->m_out << "  Number of variables....................... : " << setw(8) << i.nV << "\n";
      this->m_out << "  Number of equality constraints............ : " << setw(8) << i.nE << "\n";
      this->m_out << "  Number of inequality constraints.......... : " << setw(8) << i.nI << "\n\n";

      this->m_out << "Problem sparsity\n================\n";
      this->m_out << "  Nonzeros in Hessian of Lagrangian......... : " << setw(8) << z.Hnnz << "\n";
      this->m_out << "  Nonzeros in equality constraint Jacobian.. : " << setw(8) << z.JEnnz << "\n";
      this->m_out << "  Nonzeros in inequality constraint Jacobian : " << setw(8) << z.JInnz << "\n\n";
    }

    // ------------------------------------------------------------------------
    // Break printing (every 20 iterations)
    // ------------------------------------------------------------------------
    void print_break(Counter<Real> const & c)
    {
      if (c.k % 20 == 0)
        this->m_out << this->m_line_break << "\n" << this->m_quantities << "\n" << this->m_line_break << "\n";
    }

    // ------------------------------------------------------------------------
    // Direction info
    // ------------------------------------------------------------------------
    template<typename Z, typename D>
    void print_direction(Z const & z, D const & d)
    {
      this->m_out << showpos << scientific << setprecision(4)
        << z.phi << "  " << z.kkt[2] << " | "
        << z.shift << "  " << d.x_norm << "  " << d.this->m_line_breaknorm << "  "
        << d.ltred << "  " << d.qtred << "  " << d.m << " | ";
    }

    // ------------------------------------------------------------------------
    // Iterate info
    // ------------------------------------------------------------------------
    template<typename Counters, typename Z>
    void print_iterate(Counters const & c, Z const & z)
    {
      this->m_out << setw(5) << c.k << " | "
        << showpos << scientific << setprecision(4)
        << z.f << "  " << z.v << " | "
        << z.rho << "  " << z.mu << "  " << z.kkt[1] << " | ";
    }

    // ------------------------------------------------------------------------
    // Acceptance info
    // ------------------------------------------------------------------------
    template<typename A>
    void print_acceptance(A const & a)
    {
      this->m_out << scientific << setprecision(4)
        << a.p << "  " << a.d;
      if (a.s == 1) this->m_out << " SOC";
      this->m_out << "\n";
    }

    // ------------------------------------------------------------------------
    // Footer printing
    // ------------------------------------------------------------------------
    template<typename P, typename I, typename C, typename Z>
    void print_footer(P const & p, I const & i, C const & c, Z const & z)
    {
      print_iterate(c, z);
      this->m_out << this->m_footer << "\n" << this->m_line_break << "\n\n";

      int b = z.check_termination(p, i, c);

      this->m_out << "Final result\n============\n";
      switch (b) {
        case 0: this->m_out << "  EXIT: No termination message set\n"; break;
        case 1: this->m_out << "  EXIT: Optimal solution found\n"; break;
        case 2: this->m_out << "  EXIT: Infeasible stationary point found\n"; break;
        case 3: this->m_out << "  EXIT: Iteration limit reached\n"; break;
        case 4: this->m_out << "  EXIT: Invalid bounds\n"; break;
        case 5: this->m_out << "  EXIT: Function evaluation error\n"; break;
        default: this->m_out << "  EXIT: Unknown termination\n"; break;
      }

      this->m_out << "\nFinal values\n============\n";
      this->m_out << "  Objective function........................ : " << z.fu << "\n";
      this->m_out << "  Feasibility violation..................... : " << z.vu << "\n";
      this->m_out << "  Optimality error (feasibility)............ : " << z.kkt[0] << "\n";
      this->m_out << "  Optimality error (penalty)................ : " << z.kkt[1] << "\n";
      this->m_out << "  Optimality error (penalty-interior-point). : " << z.kkt[2] << "\n";
      this->m_out << "  Penalty parameter......................... : " << z.rho << "\n";
      this->m_out << "  Interior-point parameter.................. : " << z.mu << "\n\n";

      this->m_out << "Final counters\n==============\n";
      this->m_out << "  Iterations................................ : " << c.k << "\n";
      this->m_out << "  Function evaluations...................... : " << c.f << "\n";
      this->m_out << "  Gradient evaluations...................... : " << c.g << "\n";
      this->m_out << "  Hessian evaluations....................... : " << c.H << "\n";
      this->m_out << "  Matrix factorizations..................... : " << c.M << "\n";

      Real seconds{duration_cast<seconds>(steady_clock::now() - this->m_start).count()};
      this->m_out << "  CPU seconds............................... : " << seconds << "\n";
    }

    // ------------------------------------------------------------------------
    // Terminate output (close stream)
    // ------------------------------------------------------------------------
    void terminate()
    {
      if (this->m_out.ithis->m_outopen()) this->m_out.close();
    }
  };

} // namespace Pipal
