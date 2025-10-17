/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzinput.                                      *
 *                                                                                               *
 * The Pipal project is distributed under the MIT License.                                       *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_PIPAL_DIRECTION_HH
#define INCLUDE_PIPAL_DIRECTION_HH

// Pipal includes
#include "Pipal/Defines.hh"

namespace Pipal
{

  template<typename Real>
  struct Direction
  {
    using Vector = Eigen::Vector<Real, Eigen::Dynamic>;
    using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    using Array  = Eigen::Array<Real, Eigen::Dynamic, 1>;

    static constexpr Real INF{std::numeric_limits<Real>::infinity()}; /*!< Infinity value. */

    Vector x;       // Primal direction
    Real   x_norm;  // Primal direction norm value
    Real   x_norm_; // Primal direction norm last value
    Vector r1;      // Equality constraint slack direction
    Vector r2;      // Equality constraint slack direction
    Vector lE;      // Equality constraint multiplier direction
    Vector s1;      // Inequality constraint slack direction
    Vector s2;      // Inequality constraint slack direction
    Vector lI;      // Inequality constraint multiplier direction
    Real   l_norm;  // Constraint multiplier direction norm
    Real   lred0;   // Penalty-interior-point linear model value for zero penalty parameter
    Real   ltred0;  // Penalty-interior-point linear model reduction value for zero penalty parameter
    Real   ltred;   // Penalty-interior-point linear model reduction value
    Real   qtred;   // Penalty-interior-point quadratic model reduction value
    Real   m;       // Quality function value

    // Constructor
    Direction() {
      // Initialize last direction norm
      this->x_norm_ = INF;
    }

    // Evaluate linear combination of directions
    void evalLinearCombination(Input<Real> const & i, Direction<Real> const & d1, Direction<Real> const & d2,
      Direction<Real> const & d3, Real const a0, Real const a1, Real const a2)
    {
      // Evaluate linear combinations
      d.x = a0*d1.x  + a1*d2.x  + a2*d3.x;
      if (i.nE > 0) {
        d.r1 = a0*d1.r1 + a1*d2.r1 + a2*d3.r1;
        d.r2 = a0*d1.r2 + a1*d2.r2 + a2*d3.r2;
      }
      if (i.nI > 0) {
        d.s1 = a0*d1.s1 + a1*d2.s1 + a2*d3.s1;
        d.s2 = a0*d1.s2 + a1*d2.s2 + a2*d3.s2;
      }
      if (i.nE > 0) {
        d.lE = a0*d1.lE + a1*d2.lE + a2*d3.lE;
      }
      if (i.nI > 0) {
        d.lI = a0*d1.lI + a1*d2.lI + a2*d3.lI;
      }

      // Evaluate primal direction norm
      d.x_norm = d.x.norm();

      // Evaluate dual direction norm
      d.l_norm = Vector(d.lE, d.lI).norm(); // OPTIMIZE
    }

    // Evaluate model and model reductions
    void evalModels(Input<Real> const & i, Iterate<Real> const & z)
    {
      // Evaluate reduction in linear model of penalty-interior-point objective for zero penalty parameter
      d.lred0 = 0;
      if (i.nE > 0) {d.lred0 = d.lred0 - (Vector(1-z.mu/z.r1, 1-z.mu/z.r2).array()*Vector(d.r1, d.r2).array()).sum();}
      if (i.nI > 0) {d.lred0 = d.lred0 - (Vector(0-z.mu/z.s1, 1-z.mu/z.s2).array()*Vector(d.s1, d.s2).array()).sum();}

      // Evaluate remaining quantities only for nonzero penalty parameter
      if (z.rho > 0)
      {
        // Evaluate reduction in linear model of merit function for zero penalty parameter
        d.ltred0 = 0;
        if (i.nE > 0) {
          d.ltred0 = d.ltred0 - 0.5*full(sum(((1-z.mu/z.r1).*(-1+z.cE./(std::sqrt(z.cE.array().square().matrix() + z.mu*z.mu)))+(1-z.mu./z.r2).*(1+z.cE./(std::sqrt(z.cE.array().square().matrix() + z.mu*z.mu)))).*(z.JE*d.x)));
        }
        if (i.nI > 0) {
          d.ltred0 = d.ltred0 - 0.5*full(sum(((0-z.mu/z.s1).*(-1+z.cI./(std::sqrt(z.cI.array().square().matrix() + 4*z.mu*z.mu)))+(1-z.mu./z.s2).*(1+z.cI./(std::sqrt(z.cI.array().square().matrix() + 4*z.mu*z.mu)))).*(z.JI*d.x)));
        }

        // Evaluate reduction in linear model of merit function
        d.ltred = -z.rho*z.g.transpose()*d.x + d.ltred0;

        // Evaluate reduction in quadratic model of merit function
        d.qtred = d.ltred - 0.5*d.x.transpose()*z.H*d.x;
        if (i.nE > 0) {Jd = z.JE*d.x; Dinv = z.r1./(1+z.lE)+z.r2./(1-z.lE); d.qtred = d.qtred - 0.5*Jd.transpose()*(Jd./Dinv);}
        if (i.nI > 0) {Jd = z.JI*d.x; Dinv = z.s1./(0+z.lI)+z.s2./(1-z.lI); d.qtred = d.qtred - 0.5*Jd.transpose()*(Jd./Dinv);}

        // Initialize quality function vector
        vec.setZero(i.nV+2*i.nE+2*i.nI);

        // Set gradient of objective
        vec(Eigen::seq(0, i.nV-1)) = z.rho*z.g;

        // Set gradient of Lagrangian for constraints
        if (i.nE > 0) {vec(Eigen::seq(0, i.nV-1)) = vec(Eigen::seq(0, i.nV-1)) + ((z.lE+d.lE).transpose()*z.JE).transpose();}
        if (i.nI > 0) {vec(Eigen::seq(0, i.nV-1)) = vec(Eigen::seq(0, i.nV-1)) + ((z.lI+d.lI).transpose()*z.JI).transpose();}

        // Set complementarity for constraint slacks
        if (i.nE > 0) {vec(Eigen::seq(i.nV, i.nV+2*i.nE-1)) <<
          ((z.r1+d.r1).array() * (1 + (z.lE+d.lE)).array()).matrix(),
          ((z.r2+d.r2).array() * (1 - (z.lE+d.lE)).array()).matrix();
        }
        if (i.nI > 0) {vec(Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1)) <<
          ((z.s1+d.s1).array() * (0 + (z.lI+d.lI)).array()).matrix(),
          ((z.s2+d.s2).array() * (1 - (z.lI+d.lI)).array()).matrix();
        }

        // Evaluate quality function
        d.m = vec.array().abs().maxCoeff();
      }
    }

    // Evaluate Newton step
    void evalNewtonStep(Input<Real> const & i, Iterate<Real> const & z)
    {
      // Evaluate direction
      // dir = z.AS(:,z.AP)*(z.AL'\(z.AD\(z.AL\(z.AS(z.AP,:)*(-z.b))))); // FIXME

      // Parse direction
      d.x = dir(Eigen::seq(0, i.nV-1));
      if (i.nE > 0) {
        d.r1 = dir(Eigen::seq(i.nV, i.nV+i.nE-1));
        d.r2 = dir(Eigen::seq(i.nV+i.nE, i.nV+i.nE+i.nE-1));
      }
      if (i.nI > 0) {
        d.s1 = dir(Eigen::seq(i.nV+i.nE+i.nE, i.nV+i.nE+i.nE+i.nI-1));
        d.s2 = dir(Eigen::seq(i.nV+i.nE+i.nE+i.nI, i.nV+i.nE+i.nE+i.nI+i.nI-1));
      }
      if (i.nE > 0) {
        d.lE = dir(Eigen::seq(i.nV+i.nE+i.nE+i.nI+i.nI, i.nV+i.nE+i.nE+i.nI+i.nI+i.nE-1));
      }
      if (i.nI > 0) {
        d.lI = dir(Eigen::seq(i.nV+i.nE+i.nE+i.nI+i.nI+i.nE, i.nV+i.nE+i.nE+i.nI+i.nI+i.nE+i.nI-1));
      }

      // Evaluate primal direction norm
      d.x_norm = d.x.norm();

      // Evaluate dual direction norm
      d.l_norm = Vector(d.lE, d.lI).norm(); // OPTIMIZE
    }

    // Evaluate search direction quantities
    void evalStep(Parameter<Real> const & p, Input<Real> const & i, Counter const & c,
      Iterate<Real> const & z, Acceptance<Real> const & a)
    {
      // Reset maximum exponent for interior-point parameter increases
      p.resetMuMaxExp();

      // Update penalty-interior-point parameters based on KKT errors
      z.updateParameters(p, i);

      // Evaluate matrices
      z.evalMatrices(p, i, c);

      // Set last penalty parameter
      z.setRhoLast(z.rho);

      // Check for aggressive algorithm
      if (p.algorithm == 1)
      {
        // Check KKT memory for potential mu increase limit
        if (z.kkt(1) > std::max(z.kkt_)) {p.setMuMaxExpZero();}

        // Store current penalty and interior-point parameters
        rho_curr = z.rho; mu_curr = z.mu;

        // Evaluate trial steps
        d.evalTrialSteps(i, z, d1, d2, d3);

        // Set trial interior-point parameter values
        Eigen::ArrayXd exponents(p.mu_trials);
        for (Integer i{0}; i < p.mu_trials; ++i) {exponents[i] = (p.mu_trials - 1 - i) - p.mu_max_exp;}

        Array Mu(mu_curr * exponents.unaryExpr([p.mu_factor](Real e){return std::pow(p.mu_factor, e);}));
        Mu = Mu.min(p.mu_max).max(p.mu_min);

        // Initialize feasibility direction data
        lred0_0_mu.setZero(p.mu_trials);

        // Loop through interior-point parameter values
        for (Integer j{0}, j < p.mu_trials, ++j)
        {
          // Set penalty and interior-point parameters
          z.setRho(0); z.setMu(Mu(j));

          // Evaluate direction
          d.evalLinearCombination(i, d1, d2, d3, (z.rho/rho_curr+z.mu/mu_curr-1), (1-z.mu/mu_curr), (1-z.rho/rho_curr));

          // Cut length
          d.x = std::min(d.x_norm_/std::max(d.x_norm, 1), 1)*d.x;

          // Run fraction-to-boundary
          a.fractionToBoundary(p, i, z, d);

          // Cut length
          d.evalTrialStepCut(i, a);

          // Evaluate models
          d.evalModels(i, z);

          // Set feasibility direction data
          lred0_0_mu(j) = d.lred0;
        }

        // Initialize updating data
        ltred0_rho_mu.setZero(p.mu_trials);
        qtred_rho_mu.setZero(p.mu_trials);
        m_rho_mu.setZero(p.mu_trials);

        // Initialize check
        check = 0;

        // Loop through penalty parameter values
        for (Integer k{0}, k < p.rho_trials, ++k)
        {
          // Set penalty parameter
          z.setRho(std::max(p.rho_min, std::pow(p.rho_factor, k-1)*rho_curr));

          // Set last penalty parameter
          if (rho_curr > z.kkt(0)*z.kkt(0)) {z.setRhoLast(z.rho);}

          // Loop through interior-point parameter values
          for (Integer j{0}, j < p.mu_trials, ++j)
          {
            // Set interior-point parameter
            z.setMu(Mu(j));

            // Evaluate direction
            d.evalLinearCombination(i, d1, d2, d3, (z.rho/rho_curr+z.mu/mu_curr-1), (1-z.mu/mu_curr), (1-z.rho/rho_curr));

            // Run fraction-to-boundary
            a.fractionToBoundary(p, i, z, d);

            // Cut steps
            d.evalTrialStepCut(i, a);

            // Evaluate models
            d.evalModels(i, z);

            // Set updating data
            ltred0_rho_mu(j) = d.ltred0;
            qtred_rho_mu(j)  = d.qtred;
            m_rho_mu(j)      = d.m;

            // Check updating conditions for infeasible points
            if (z.v > p.opt_err_tol && (ltred0_rho_mu(j) < p.update_con_1*lred0_0_mu(j) || qtred_rho_mu(j) < p.update_con_2*lred0_0_mu(j) || z.rho > z.kkt(0)*z.kkt(0))) {m_rho_mu(j) = INF;}

            // Check updating conditions for feasible points
            if (z.v <= p.opt_err_tol && qtred_rho_mu(j) < 0) {m_rho_mu(j) = INF;}

          }

          // Find minimum m for current rho
          m_min = std::min(m_rho_mu);

          // Check for finite minimum
          if (m_min < INF)
          {
            // Loop through mu values
            for (Integer j{0}, j < p.mu_trials, ++j)
            {
              // Set mu
              mu = Mu(j);

              // Check condition
              if (m_rho_mu(j) <= p.update_con_3*m_min) {z.setMu(mu);}
            }

            // Set condition check
            check = 1;

            // Break loop
            break;
          }
        }

        // Check conditions
        if (check == 0) {z.setRho(rho_curr); z.setMu(mu_curr);}

        // Evaluate merit
        z.evalMerit(i);
      }

      // Evaluate primal-dual right-hand side vector
      z.evalNewtonRhs(i);

      // Evaluate search direction
      d.evalNewtonStep(i, z);

      // Evaluate models
      d.evalModels(i, z);

      // Store last direction norm
      d.x_norm_ = d.x_norm;
    }

    // Evaluate and store trial step
    void evalTrialStep(Input<Real> const & i, Direction<Real> const & v)
    {
      // Set direction components
      v.x = d.x;
      if (i.nE > 0) {v.r1 = d.r1; v.r2 = d.r2;}
      if (i.nI > 0) {v.s1 = d.s1; v.s2 = d.s2;}
      if (i.nE > 0) {v.lE = d.lE;}
      if (i.nI > 0) {v.lI = d.lI;}
    }

    // Evaluate trial step cut by fraction-to-boundary rule
    function evalTrialStepCut(Input<Real> const & i, Acceptance<Real> const & a)
    {
      // Set direction components
      d.x = a.p*d.x ;
      if (i.nE > 0) {d.r1 = a.p*d.r1; d.r2 = a.p*d.r2;}
      if (i.nI > 0) {d.s1 = a.p*d.s1; d.s2 = a.p*d.s2;}
      if (i.nE > 0) {d.lE = a.d*d.lE;}
      if (i.nI > 0) {d.lI = a.d*d.lI;}
    }

    // Evaluate and store directions for parameter combinations
    void evalTrialSteps(i, z, d1, d2, d3)
    {
      // Store current penalty and interior-point parameters
      rho_curr = z.rho;
      mu_curr  = z.mu;

      // Evaluate direction for current penalty and interior-point parameters
      z.setRho(rho_curr);
      z.setMu(mu_curr);
      z.evalNewtonRhs(i);
      d.evalNewtonStep(i, z);
      d.evalTrialStep(i, d1);

      // Evaluate direction for zero interior-point parameter
      z.setRho(rho_curr);
      z.setMu(0);
      z.evalNewtonRhs(i);
      d.evalNewtonStep(i, z);
      d.evalTrialStep(i, d2);

      // Evaluate direction for zero penalty parameter
      z.setRho(0);
      z.setMu(mu_curr);
      z.evalNewtonRhs(i);
      d.evalNewtonStep(i, z);
      d.evalTrialStep(i, d3);
    }

    // Set direction
    void setDirection(i,dx,dr1,dr2,ds1,ds2,dlE,dlI,dx_norm,dl_norm)
    {
      // Set primal variables
      d.x = dx;
      if (i.nE > 0) {d.r1 = dr1; d.r2 = dr2; d.lE = dlE;}
      if (i.nI > 0) {d.s1 = ds1; d.s2 = ds2; d.lI = dlI;}
      d.x_norm = dx_norm;
      d.l_norm = dl_norm;
    }

  }; // struct Direction

} // namespace Pipal

#endif /* INCLUDE_PIPAL_DIRECTION_HH */
