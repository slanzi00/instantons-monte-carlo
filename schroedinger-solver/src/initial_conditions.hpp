#ifndef SOLVER_INIT_CONDITIONS_HPP
#define SOLVER_INIT_CONDITIONS_HPP

#include <cmath>
#include <iostream>

#include <unsupported/Eigen/Polynomials>

#include "potentials.hpp"

struct x_limits
{
  double x_min;
  double x_max;
};

template <uint16_t polynomial_potential_degree>
class InitialConditions
{
  double m_max_energy;
  double m_step;
  x_limits m_boundaries;
  PolynomialPotential<polynomial_potential_degree> m_potential;

 public:
  InitialConditions(PolynomialPotential<polynomial_potential_degree> potential, double max_energy)
      :  m_max_energy{max_energy}, m_step{1. / (10. * std::sqrt(2. * m_max_energy))}, m_potential{std::move(potential)}
  {
    auto coeff = m_potential.coefficients;
    coeff(0) -= m_max_energy;
    Eigen::PolynomialSolver<double, polynomial_potential_degree> polynomial_solver(coeff);
    auto greatets_root = polynomial_solver.greatestRoot();
    m_boundaries = {-greatets_root.real(), greatets_root.real()};
  }

  double get_step()
  {
    return m_step;
  }

  uint16_t n_points()
  {
    return std::round(2. * (m_boundaries.x_max / m_step));
  }

  double x_min()
  {
    return m_boundaries.x_min;
  }

  double x_max()
  {
    return m_boundaries.x_max;
  }

  double kinetic()
  {
    return 1. / (2. * std::pow(m_step, 2));
  }

  double potential(double x)
  {
    return m_potential(x);
  }
};

#endif