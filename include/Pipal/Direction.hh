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

    Vector x;       /*!< Primal direction. */
    Real   x_norm;  /*!< Primal direction norm value. */
    Real   x_norm_; /*!< Primal direction norm last value. */
    Vector r1;      /*!< Equality constraint slack direction. */
    Vector r2;      /*!< Equality constraint slack direction. */
    Vector lE;      /*!< Equality constraint multiplier direction. */
    Vector s1;      /*!< Inequality constraint slack direction. */
    Vector s2;      /*!< Inequality constraint slack direction. */
    Vector lI;      /*!< Inequality constraint multiplier direction. */
    Real   l_norm;  /*!< Constraint multiplier direction norm. */
    Real   lred0;   /*!< Penalty-interior-point linear model value for zero penalty parameter. */
    Real   ltred0;  /*!< Penalty-interior-point linear model reduction value for zero penalty parameter. */
    Real   ltred;   /*!< Penalty-interior-point linear model reduction value. */
    Real   qtred;   /*!< Penalty-interior-point quadratic model reduction value. */
    Real   m;       /*!< Quality function value. */

  // Constructor
  Direction() {
    // Initialize last direction norm
    this->x_norm_ = INF;
  }

  // Evaluate linear combination of directions
  void evalLinearCombination(Input<Real> const & input, Direction<Real> const & d1, Direction<Real> const & d2,
    Direction<Real> const & d3, Vector const & a)
  {
    // Evaluate linear combinations
    this->x = a(1)*d1.x + a(2)*d2.x + a(3)*d3.x;
    if (input.nE > 0) {
      this->r1 = a(1)*d1.r1 + a(2)*d2.r1 + a(3)*d3.r1;
      this->r2 = a(1)*d1.r2 + a(2)*d2.r2 + a(3)*d3.r2;
    }
    if (input.nI > 0) {
      this->s1 = a(1)*d1.s1 + a(2)*d2.s1 + a(3)*d3.s1;
      this->s2 = a(1)*d1.s2 + a(2)*d2.s2 + a(3)*d3.s2;
    }
    if (input.nE > 0) {
      this->lE = a(1)*d1.lE + a(2)*d2.lE + a(3)*d3.lE;
    }
    if (input.nI > 0) {
      this->lI = a(1)*d1.lI + a(2)*d2.lI + a(3)*d3.lI;
    }

    // Evaluate primal direction norm
    this->x_norm = this->x.norm();

    // Evaluate dual direction norm
    this->l_norm = norm([this->lE; this->lI]);
  }

  // Evaluate model and model reductions
  void eval_models(i,z)
  {
    // Evaluate reduction in linear model of penalty-interior-point objective for zero penalty parameter
    this->lred0 = 0;
    if (input.nE > 0) {
      this->lred0 = this->lred0 - sum([1-z.mu./z.r1; 1-z.mu./z.r2].*[this->r1; this->r2]);
    }
    if (input.nI > 0) {
      this->lred0 = this->lred0 - sum([0-z.mu./z.s1; 1-z.mu./z.s2].*[this->s1; this->s2]);
    }

    // Evaluate remaining quantities only for nonzero penalty parameter
    if (z.rho > 0)
    {
      // Evaluate reduction in linear model of merit function for zero penalty parameter
      this->ltred0 = 0;
      if input.nE > 0, this->ltred0 = this->ltred0 - 0.5*full(sum(((1-z.mu./z.r1).*(-1+z.cE./(sqrt(z.cE.^2 +   z.mu^2)))+(1-z.mu./z.r2).*(1+z.cE./(sqrt(z.cE.^2 +   z.mu^2)))).*(z.JE*this->x))); end;
      if input.nI > 0, this->ltred0 = this->ltred0 - 0.5*full(sum(((0-z.mu./z.s1).*(-1+z.cI./(sqrt(z.cI.^2 + 4*z.mu^2)))+(1-z.mu./z.s2).*(1+z.cI./(sqrt(z.cI.^2 + 4*z.mu^2)))).*(z.JI*this->x))); end;

      // Evaluate reduction in linear model of merit function
      this->ltred = -z.rho*z.g.transpose() * this->x + this->ltred0;

      // Evaluate reduction in quadratic model of merit function
      this->qtred = this->ltred - 0.5*this->x.transpose()*z.H*this->x;
      if (input.nE > 0) {Jd = z.JE*this->x; Dinv = z.r1./(1+z.lE)+z.r2./(1-z.lE); this->qtred = this->qtred - 0.5*Jthis->transpose()*(Jthis->/Dinv);}
      if (input.nI > 0) {Jd = z.JI*this->x; Dinv = z.s1./(0+z.lI)+z.s2./(1-z.lI); this->qtred = this->qtred - 0.5*Jthis->transpose()*(Jthis->/Dinv);}

      // Initialize quality function vector
      vec.setZero(input.nV+2*input.nE+2*input.nI);

      // Set gradient of objective
      vec(1:input.nV) = z.rho*z.g;

      // Set gradient of Lagrangian for constraints
      if (input.nE > 0) {vec(Eigen::seq(1, input.nV)) += ((z.lE+this->lE).transpose() * z.JE).transpose();}
      if (input.nI > 0) {vec(Eigen::seq(1, input.nV)) += ((z.lI+this->lI).transpose() * z.JI).transpose();}

      // Set complementarity for constraint slacks
      if (input.nE > 0) {vec(1+input.nV       :input.nV+2*input.nE       ) = [(z.r1+this->r1).*(1 + (z.lE+this->lE)); (z.r2+this->r2).*(1 - (z.lE+this->lE))];}
      if (input.nI > 0) {vec(1+input.nV+2*input.nE:input.nV+2*input.nE+2*input.nI) = [(z.s1+this->s1).*(0 + (z.lI+this->lI)); (z.s2+this->s2).*(1 - (z.lI+this->lI))];}

      // Evaluate quality function
      this->m = norm(vec, inf);
    }
  }

  // Evaluate Newton step
  void eval_newton_step(i,z)
  {
    // Evaluate direction
    //dir = z.AS(:,z.AP)*(z.AL'\(z.AD\(z.AL\(z.AS(z.AP,:)*(-z.b)))));

    // Parse direction
    Integer idx_ini{0}, idx_end{input.nV-1};
    this->x = dir(Eigen::seq(0,input.nV-1));
    idx_ini += input.nV; idx_end += input.nE;
    if (input.nE > 0) {
      this->r1 = dir(Eigen::seq(idx_ini, idx_end));
      this->r2 = dir(Eigen::seq(idx_ini+input.nE, idx_end+input.nE));
    }
    idx_ini += 2*input.nE; idx_end += 2*input.nE+input.nI;
    if (input.nI > 0) {
      this->s1 = dir(Eigen::seq(idx_ini:input.nV+input.nE+input.nE+input.nI);
      this->s2 = dir(Eigen::seq(idx_ini+input.nI:input.nV+input.nE+input.nE+input.nI+input.nI);
    }
    if (input.nE > 0) {
      this->lE = dir(1+input.nV+input.nE+input.nE+input.nI+input.nI:input.nV+input.nE+input.nE+input.nI+input.nI+input.nE);
    }
    if (input.nI > 0) {
      this->lI = dir(1+input.nV+input.nE+input.nE+input.nI+input.nI+input.nE:input.nV+input.nE+input.nE+input.nI+input.nI+input.nE+input.nI);
    }

    // Evaluate primal direction norm
    this->x_norm = this->x.norm();

    // Evaluate dual direction norm
    this->l_norm = norm([this->lE;this->lI]);
  }

  // Evaluate search direction quantities
  void eval_step(p,i,c,z,a)
  {
    // Reset maximum exponent for interior-point parameter increases
    p.reset_muMaxExp;

    // Update penalty-interior-point parameters based on KKT errors
    z.update_parameters(p,i);

    // Evaluate matrices
    z.eval_matrices(p,i,c);

    // Set last penalty parameter
    z.set_rho_last(z.rho);

    // Check for aggressive algorithm
    if (p.algorithm == 1)
    {
      // Check KKT memory for potential mu increase limit
      if (z.kkt(2) > max(z.kkt_)) {p.mu_max_exp(0.0);}

      // Store current penalty and interior-point parameters
      rho_curr = z.rho; mu_curr = z.mu;

      // Evaluate trial steps
      this->eval_trial_steps(i,z,d1,d2,d3);

      // Set trial interior-point parameter values
      Mu = max(p.mu_min,min(p.mu_factor.^([p.mu_trials-1:-1:0]-p.mu_max_exp)*mu_curr,p.mu_max));

      // Initialize feasibility direction data
      lred0_0_mu = zeros(1,p.mu_trials);

      // Loop through interior-point parameter values
      for (j = 1:p.mu_trials)
      {
        // Set penalty and interior-point parameters
        z.set_rho(0.0); z.set_mu(Mu(j));

        // Evaluate direction
        this->evalLinearCombination(i,d1,d2,d3,[(z.rho/rho_curr+z.mu/mu_curr-1),(1-z.mu/mu_curr),(1-z.rho/rho_curr)]);

        // Cut length
        this->x = min(this->x_norm_/max(this->x_norm,1),1)*this->x;

        // Run fraction-to-boundary
        a.fraction_to_boundary(p,i,z,d);

        // Cut length
        this->eval_trials_step_cut(i,a);

        // Evaluate models
        this->eval_models(i,z);

        // Set feasibility direction data
        lred0_0_mu(j) = this->lred0;
      }

      // Initialize updating data
      ltred0_rho_mu.setZero(p.mu_trials);
      qtred_rho_mu.setZero(p.mu_trials);
      m_rho_mu.setZero(p.mu_trials);

      // Initialize check
      Integer check{0};

      // Loop through penalty parameter values
      for (k = 1:p.rho_trials)
      {
        // Set penalty parameter
        z.set_rho(max(p.rho_min,(p.rho_factor^(k-1))*rho_curr));

        // Set last penalty parameter
        if (rho_curr > z.kkt(1)*z.kkt(1)) {z.set_rho_last(z.rho);}

        // Loop through interior-point parameter values
        for (j = 1:p.mu_trials)
        {
          // Set interior-point parameter
          z.set_mu(Mu(j));

          // Evaluate direction
          this->evalLinearCombination(i,d1,d2,d3,[(z.rho/rho_curr+z.mu/mu_curr-1),(1-z.mu/mu_curr),(1-z.rho/rho_curr)]);

          // Run fraction-to-boundary
          a.fraction_to_boundary(p,i,z,d);

          // Cut steps
          this->eval_trials_step_cut(i,a);

          // Evaluate models
          this->eval_models(i,z);

          // Set updating data
          ltred0_rho_mu(j) = this->ltred0;
          qtred_rho_mu(j)  = this->qtred;
          m_rho_mu(j)      = this->m;

          // Check updating conditions for infeasible points
          if (z.v > p.opt_err_tol && (ltred0_rho_mu(j) < p.update_con_1*lred0_0_mu(j) ||
              qtred_rho_mu(j) < p.update_con_2*lred0_0_mu(j) || z.rho > z.kkt(1)^2)) {m_rho_mu(j) = INF;};

          // Check updating conditions for feasible points
          if (z.v <= p.opt_err_tol && qtred_rho_mu(j) < 0) {m_rho_mu(j) = INF;}

        }

        // Find minimum m for current rho
        m_min = min(m_rho_mu);

        // Check for finite minimum
        if (m_min < inf)
        {
          // Loop through mu values
          for (j = 1:p.mu_trials)
          {
            // Set mu
            mu = Mu(j);

            // Check condition
            if (m_rho_mu(j) <= p.update_con_3*m_min) {z.set_mu(mu);}
          }

          // Set condition check
          check = 1;

          // Break loop
          break;
        }
      }

      // Check conditions
      if (check == 0) {z.set_rho(rho_curr); z.set_mu(mu_curr);}

      // Evaluate merit
      z.evalMerit(i);
      }

      // Evaluate primal-dual right-hand side vector
      z.eval_newton_rhs(i);

      // Evaluate search direction
      this->eval_newton_step(i,z);

      // Evaluate models
      this->eval_models(i,z);

      // Store last direction norm
      this->x_norm_ = this->x_norm;
    }

    // Evaluate and store trial step
    void eval_trial_step(i,d)
    {
      // Set direction components
      d.x = this->x;
      if (i.nE > 0) {d.r1 = this->r1; d.r2 = this->r2;}
      if (i.nI > 0) {d.s1 = this->s1; d.s2 = this->s2;}
      if (i.nE > 0) {d.lE = this->lE;}
      if (i.nI > 0) {d.lI = this->lI;}
    }

    // Evaluate trial step cut by fraction-to-boundary rule
    void eval_trial_step_cut(i,a)
    {
      // Set direction components
      this->x = a.p*this->x;
      if (i.nE > 0) {this->r1 = a.p*this->r1; this->r2 = a.p*this->r2;}
      if (i.nI > 0) {this->s1 = a.p*this->s1; this->s2 = a.p*this->s2;}
      if (i.nE > 0) {this->lE = a.d*this->lE;}
      if (i.nI > 0) {this->lI = a.d*this->lI;}
    }

    // Evaluate and store directions for parameter combinations
    void eval_trial_steps(d,i,z,d1,d2,d3)

      // Store current penalty and interior-point parameters
      rho_curr = z.rho;
      mu_curr  = z.mu;

      // Evaluate direction for current penalty and interior-point parameters
      z.set_rho(rho_curr);
      z.set_mu(mu_curr);
      z.eval_newton_rhs(i);
      this->eval_newton_step(i, z);
      this->eval_trials_step(i, d1);

      // Evaluate direction for zero interior-point parameter
      z.set_rho(rho_curr);
      z.set_mu(0.0);
      z.eval_newton_rhs(i);
      this->eval_newton_step(i, z);
      this->eval_trials_step(i, d2);

      // Evaluate direction for zero penalty parameter
      z.set_rho(0.0);
      z.set_mu(mu_curr);
      z.eval_newton_rhs(i);
      this->eval_newton_step(i, z);
      this->eval_trials_step(i, d3);
    }

    // Set direction
    void setDirection(d,i,dx,dr1,dr2,ds1,ds2,dlE,dlI,dx_norm,dl_norm)
    {
      // Set primal variables
      this->x = dx;
      if (input.nE > 0) {this->r1 = dr1; this->r2 = dr2; this->lE = dlE;}
      if (input.nI > 0) {this->s1 = ds1; this->s2 = ds2; this->lI = dlI;}
      this->x_norm = dx_norm;
      this->l_norm = dl_norm;
    }

  }; // struct Direction

} // namespace Pipal

#endif /* INCLUDE_PIPAL_DIRECTION_HH */
