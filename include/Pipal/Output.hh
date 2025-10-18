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
#include "Pipal/Types.hh"

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
    static_assert(std::is_floating_point<Real>::value,
      "Pipal::Output: template argument 'Real' must be a floating-point type.");

  private:
    using String       = std::string;
    using Ostream      = std::ostream;
    using SteadyClock  = std::chrono::steady_clock;
    using DurationCast = std::chrono::duration_cast;
    using TimePoint    = std::chrono::duration_cast::time_point;
    using Seconds      = std::chrono::seconds;

    Ostream   s; // Output stream
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

    // Destructor closes stream
    ~Output() {}

    // Header printing
    void printHeader(Input const & i, Iterate const & z)
    {
      this->s << "Problem name\n============\n  " << i.id << "\n\n";
      this->s << "Problem size\n============\n";
      this->s << "  Number of variables....................... : " << std::setw(8) << i.nV << "\n";
      this->s << "  Number of equality constraints............ : " << std::setw(8) << i.nE << "\n";
      this->s << "  Number of inequality constraints.......... : " << std::setw(8) << i.nI << "\n\n";

      this->s << "Problem sparsity\n================\n";
      this->s << "  Nonzeros in Hessian of Lagrangian......... : " << std::setw(8) << z.Hnnz << "\n";
      this->s << "  Nonzeros in equality constraint Jacobian.. : " << std::setw(8) << z.JEnnz << "\n";
      this->s << "  Nonzeros in inequality constraint Jacobian : " << std::setw(8) << z.JInnz << "\n\n";
    }

    // Break printing (every 20 iterations)
    void printBreak(Counter const & c)
    {
      if (c.k % 20 == 0) {
        this->s << this->l << "\n" << this->q << "\n" << this->l << "\n";
      }
    }

    // Direction info
    void printDirection(Iterate const & z, Direction const & d)
    {
      this->s << std::showpos << std::scientific << std::setprecision(4)
        << z.phi << "  " << z.kkt[2] << " | "
        << z.shift << "  " << d.x_norm << "  " << d.l_norm << "  "
        << d.ltred << "  " << d.qtred << "  " << d.m << " | ";
    }

    // Iterate info
    void printIterate(Counter const & c, Iterate const & z)
    {
      this->s << std::setw(5) << c.k << " | "
        << std::showpos << std::scientific << std::setprecision(4)
        << z.f << "  " << z.v << " | "
        << z.rho << "  " << z.mu << "  " << z.kkt[1] << " | ";
    }

    // Acceptance info
    void printAcceptance(Acceptance const & a)
    {
      this->s << std::scientific << std::setprecision(4)
        << a.p << "  " << a.d;
      if (a.s == 1) this->s << " SOC";
      this->s << "\n";
    }

    // Footer printing
    void printFooter(Parameter const & p, Input const & i, Counter const & c, Iterate const & z)
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

      Real seconds{std::chrono::duration_cast<Seconds>(SteadyClock::now() - this->t).count()};
      this->s << "  CPU seconds............................... : " << seconds << "\n";
    }

  }; // class Output

} // namespace Pipal

#endif // INCLUDE_PIPAL_OUTPUT_HH
