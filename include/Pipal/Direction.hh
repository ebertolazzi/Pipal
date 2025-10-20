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
#include "Pipal/Types.hh"
#include "Pipal/Acceptance.hh"
#include "Pipal/Input.hh"
#include "Pipal/Iterate.hh"

namespace Pipal
{

  // Evaluate linear combination of directions
  void evalLinearCombination(Direction & d, Input & i, Direction & d1, Direction & d2,
    Direction & d3, Real const a0, Real const a1, Real const a2)
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
    d.l_norm = std::sqrt(d.lE.squaredNorm() + d.lI.squaredNorm());
  }

  // Evaluate model and model reductions
  void evalModels(Direction & d, Input & i, Iterate & z)
  {
    // Evaluate reduction in linear model of penalty-interior-point objective for zero penalty parameter
    d.lred0 = 0;
    if (i.nE > 0) {
      Vector tmp1(1.0 - z.mu / z.r1.array());
      Vector tmp2(1.0 - z.mu / z.r2.array());
      Vector tmp3(tmp1.size() + tmp2.size());
      tmp3 << tmp1, tmp2;
      Vector tmp4(d.r1.size() + d.r2.size());
      tmp4 << d.r1, d.r2;
      d.lred0 -= (tmp3.array()*tmp4.array()).sum();
    }
    if (i.nI > 0) {
      Vector tmp1(0.0 - z.mu / z.s1.array());
      Vector tmp2(1.0 - z.mu / z.s2.array());
      Vector tmp3(tmp1.size() + tmp2.size());
      tmp3 << tmp1, tmp2;
      Vector tmp4(d.s1.size() + d.s2.size());
      tmp4 << d.s1, d.s2;
      d.lred0 -= (tmp3.array()*tmp4.array()).sum();
    }

    // Evaluate remaining quantities only for nonzero penalty parameter
    if (z.rho > 0)
    {
      // Evaluate reduction in linear model of merit function for zero penalty parameter
      d.ltred0 = 0;
      if (i.nE > 0) {
        Array sqrt_term((z.cE.array().square() + z.mu*z.mu).sqrt());
        Array expr = ((1.0 - z.mu / z.r1.array()) * (-1.0 + z.cE.array() / sqrt_term) +
          (1.0 - z.mu / z.r2.array()) * (1.0 + z.cE.array() / sqrt_term)) * (z.JE * d.x).array();
        d.ltred0 -= 0.5 * expr.sum();
      }
      if (i.nI > 0) {
        Array sqrt_term((z.cI.array().square() + z.mu*z.mu).sqrt());
        Array expr = ((1.0 - z.mu / z.s1.array()) * (-1.0 + z.cI.array() / sqrt_term) +
          (1.0 - z.mu / z.s2.array()) * (1.0 + z.cI.array() / sqrt_term)) * (z.JI * d.x).array();
        d.ltred0 -= 0.5 * expr.sum();
      }

      // Evaluate reduction in linear model of merit function
      d.ltred = -z.rho*z.g.transpose()*d.x + d.ltred0;

      // Evaluate reduction in quadratic model of merit function
      d.qtred = d.ltred - 0.5*d.x.transpose()*z.H*d.x;
      if (i.nE > 0) {
        Vector Jd(z.JE*d.x);
        Matrix Dinv((z.r1.array()/(1.0+z.lE.array()).array() + z.r2.array()/(1.0-z.lE.array()).array()).matrix());
        d.qtred -= 0.5*Jd.transpose()*((Jd.array()/Dinv.array()).matrix());
      }
      if (i.nI > 0) {
        Vector Jd(z.JI*d.x);
        Matrix Dinv((z.s1.array()/(0.0+z.lI.array()).array() + z.s2.array()/(1.0-z.lI.array()).array()).matrix());
        d.qtred -= 0.5*Jd.transpose()*((Jd.array()/Dinv.array()).matrix());
      }

      // Initialize quality function vector
      Vector vec(i.nV+2*i.nE+2*i.nI);
      vec.setZero();

      // Set gradient of objective
      vec(Eigen::seq(0, i.nV-1)) = z.rho*z.g;

      // Set gradient of Lagrangian for constraints
      if (i.nE > 0) {vec(Eigen::seq(0, i.nV-1)) = vec(Eigen::seq(0, i.nV-1)) + ((z.lE+d.lE).transpose()*z.JE).transpose();}
      if (i.nI > 0) {vec(Eigen::seq(0, i.nV-1)) = vec(Eigen::seq(0, i.nV-1)) + ((z.lI+d.lI).transpose()*z.JI).transpose();}

      // Set complementarity for constraint slacks
      if (i.nE > 0) {vec(Eigen::seq(i.nV, i.nV+2*i.nE-1)) <<
        ((z.r1+d.r1).array() * (1.0 + (z.lE+d.lE).array()).array()).matrix(),
        ((z.r2+d.r2).array() * (1.0 - (z.lE+d.lE).array()).array()).matrix();
      }
      if (i.nI > 0) {vec(Eigen::seq(i.nV+2*i.nE, i.nV+2*i.nE+2*i.nI-1)) <<
        ((z.s1+d.s1).array() * (0.0 + (z.lI+d.lI).array()).array()).matrix(),
        ((z.s2+d.s2).array() * (1.0 - (z.lI+d.lI).array()).array()).matrix();
      }

      // Evaluate quality function
      d.m = vec.array().abs().maxCoeff();
    }
  }

  // Evaluate Newton step
  void evalNewtonStep(Direction & d, Input & i, Iterate & z)
  {
    // Evaluate direction
    Vector dir(z.ldlt.solve(-z.b));

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
    d.l_norm = std::sqrt(d.lE.squaredNorm() + d.lI.squaredNorm());
  }

  // Evaluate search direction quantities
  void evalStep(Direction & d, Parameter & p, Input & i, Counter & c,
    Iterate & z, Acceptance & a)
  {
    // Reset maximum exponent for interior-point parameter increases
    resetMuMaxExp(p);

    // Update penalty-interior-point parameters based on KKT errors
    updateParameters(z, p, i);

    // Evaluate matrices
    evalMatrices(z, p, i, c);

    // Set last penalty parameter
    setRhoLast(z, z.rho);

    // Check for aggressive algorithm
    if (p.algorithm == Algorithm::ADAPTIVE)
    {
      // Check KKT memory for potential mu increase limit
      if (z.kkt(1) > z.kkt_.maxCoeff()) {setMuMaxExpZero(p);}

      // Store current penalty and interior-point parameters
      Real rho_curr{z.rho}, mu_curr{z.mu};

      // Evaluate trial steps
      Direction d1, d2, d3;
      evalTrialSteps(d, i, z, d1, d2, d3);

      // Set trial interior-point parameter values
      Array exponents(p.mu_trials);
      for (Integer i{0}; i < p.mu_trials; ++i) {exponents[i] = (p.mu_trials - 1.0 - i) - p.mu_max_exp;}

      Array Mu(mu_curr * exponents.unaryExpr([&p](Real e){return std::pow(p.mu_factor, e);}));
      Mu = Mu.min(p.mu_max).max(p.mu_min);

      // Initialize feasibility direction data
      Vector lred0_0_mu(p.mu_trials);
      lred0_0_mu.setZero();

      // Loop through interior-point parameter values
      for (Integer j{0}; j < p.mu_trials; ++j)
      {
        // Set penalty and interior-point parameters
        setRho(z, 0.0);
        setMu(z, Mu(j));

        // Evaluate direction
        evalLinearCombination(d, i, d1, d2, d3, (z.rho/rho_curr+z.mu/mu_curr-1.0), (1.0-z.mu/mu_curr), (1.0-z.rho/rho_curr));

        // Cut length
        d.x = std::min(d.x_norm_/std::max(d.x_norm, 1.0), 1.0)*d.x;

        // Run fraction-to-boundary
        fractionToBoundary(a, p, i, z, d);

        // Cut length
        evalTrialStepCut(d, i, a);

        // Evaluate models
        evalModels(d, i, z);

        // Set feasibility direction data
        lred0_0_mu(j) = d.lred0;
      }

      // Initialize updating data
      Vector ltred0_rho_mu(p.mu_trials);
      ltred0_rho_mu.setZero();
      Vector qtred_rho_mu(p.mu_trials);
      qtred_rho_mu.setZero();
      Vector m_rho_mu(p.mu_trials);
      m_rho_mu.setZero();

      // Initialize check
      bool check{false};

      // Loop through penalty parameter values
      for (Integer k{0}; k < p.rho_trials; ++k)
      {
        // Set penalty parameter
        setRho(z, std::max(p.rho_min, std::pow(p.rho_factor, k-1)*rho_curr));

        // Set last penalty parameter
        if (rho_curr > z.kkt(0)*z.kkt(0)) {setRhoLast(z, z.rho);}

        // Loop through interior-point parameter values
        for (Integer j{0}; j < p.mu_trials; ++j)
        {
          // Set interior-point parameter
          setMu(z, Mu(j));

          // Evaluate direction
          evalLinearCombination(d, i, d1, d2, d3, (z.rho/rho_curr+z.mu/mu_curr-1), (1-z.mu/mu_curr), (1-z.rho/rho_curr));

          // Run fraction-to-boundary
          fractionToBoundary(a, p, i, z, d);

          // Cut steps
          evalTrialStepCut(d, i, a);

          // Evaluate models
          evalModels(d, i, z);

          // Set updating data
          ltred0_rho_mu(j) = d.ltred0;
          qtred_rho_mu(j)  = d.qtred;
          m_rho_mu(j)      = d.m;

          // Check updating conditions for infeasible points
          if (z.v > p.opt_err_tol && (ltred0_rho_mu(j) < p.update_con_1*lred0_0_mu(j) || qtred_rho_mu(j) < p.update_con_2*lred0_0_mu(j) || z.rho > z.kkt(0)*z.kkt(0))) {m_rho_mu(j) = INFTY;}

          // Check updating conditions for feasible points
          if (z.v <= p.opt_err_tol && qtred_rho_mu(j) < 0) {m_rho_mu(j) = INFTY;}
        }

        // Find minimum m for current rho
        Real m_min{m_rho_mu.minCoeff()};

        // Check for finite minimum
        if (m_min < INFTY)
        {
          // Loop through mu values
          for (Integer j{0}; j < p.mu_trials; ++j)
          {
            // Check condition
            if (m_rho_mu(j) <= p.update_con_3*m_min) {setMu(z, Mu(j));}
          }

          // Set condition check
          check = true;

          // Break loop
          break;
        }
      }

      // Check conditions
      if (check == false) {setRho(z, rho_curr); setMu(z, mu_curr);}

      // Evaluate merit
      evalMerit(z, i);
    }

    // Evaluate primal-dual right-hand side vector
    evalNewtonRhs(z, i);

    // Evaluate search direction
    evalNewtonStep(d, i, z);

    // Evaluate models
    evalModels(d, i, z);

    // Store last direction norm
    d.x_norm_ = d.x_norm;
  }

  // Evaluate and store trial step
  void evalTrialStep(Direction & d, Input & i, Direction & v)
  {
    // Set direction components
    v.x = d.x;
    if (i.nE > 0) {v.r1 = d.r1; v.r2 = d.r2;}
    if (i.nI > 0) {v.s1 = d.s1; v.s2 = d.s2;}
    if (i.nE > 0) {v.lE = d.lE;}
    if (i.nI > 0) {v.lI = d.lI;}
  }

  // Evaluate trial step cut by fraction-to-boundary rule
  void evalTrialStepCut(Direction & d, Input & i, Acceptance & a)
  {
    // Set direction components
    d.x = a.p*d.x ;
    if (i.nE > 0) {d.r1 = a.p*d.r1; d.r2 = a.p*d.r2;}
    if (i.nI > 0) {d.s1 = a.p*d.s1; d.s2 = a.p*d.s2;}
    if (i.nE > 0) {d.lE = a.d*d.lE;}
    if (i.nI > 0) {d.lI = a.d*d.lI;}
  }

  // Evaluate and store directions for parameter combinations
  void evalTrialSteps(Direction & d, Input & i, Iterate & z, Direction & d1,
    Direction & d2, Direction & d3)
  {
    // Store current penalty and interior-point parameters
    Real rho_curr{z.rho}, mu_curr{z.mu};

    // Evaluate direction for current penalty and interior-point parameters
    setRho(z, rho_curr);
    setMu(z, mu_curr);
    evalNewtonRhs(z, i);
    evalNewtonStep(d, i, z);
    evalTrialStep(d, i, d1);

    // Evaluate direction for zero interior-point parameter
    setRho(z, rho_curr);
    setMu(z, 0);
    evalNewtonRhs(z, i);
    evalNewtonStep(d, i, z);
    evalTrialStep(d, i, d2);

    // Evaluate direction for zero penalty parameter
    setRho(z, 0);
    setMu(z, mu_curr);
    evalNewtonRhs(z, i);
    evalNewtonStep(d, i, z);
    evalTrialStep(d, i, d3);
  }

  // Set direction
  void setDirection(Direction & d, Input & i, Vector const & dx, Vector const & dr1, Vector const & dr2,
    Vector const & ds1, Vector const & ds2, Vector const & dlE, Vector const & dlI, Real const dx_norm, Real const dl_norm)
  {
    // Set primal variables
    d.x = dx;
    if (i.nE > 0) {d.r1 = dr1; d.r2 = dr2; d.lE = dlE;}
    if (i.nI > 0) {d.s1 = ds1; d.s2 = ds2; d.lI = dlI;}
    d.x_norm = dx_norm;
    d.l_norm = dl_norm;
  }

} // namespace Pipal

#endif /* INCLUDE_PIPAL_DIRECTION_HH */
