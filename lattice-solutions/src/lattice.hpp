#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

#include "setting_values.hpp"

struct AnharmonicPotential
{
  double eta;
  double operator()(double x);
};

struct Lattice
{
  // potential choose
  AnharmonicPotential m_potential;
  
  // define lattice
  std::array<double, sv::n_lattice_points> positions;
  std::array<double, sv::n_lattice_points> euclidean_time;

  Lattice(AnharmonicPotential const& potential);

  // calculate action for metropolis sweeps: 
  // in this case we compute the action between neigh positions
  double calculate_action(size_t position_index);

  // there is also the possibility to compute the full action
  double calculate_complete_action();
};

#endif