#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

struct AnharmonicPotential
{
  double eta;
  double operator()(double x)
  {
    return std::pow(x * x - eta * eta, 2);
  }
};

template <size_t n_points>
struct Lattice
{
  double m_lattice_spacing;
  AnharmonicPotential m_potential;
  std::array<double, n_points> positions;
  std::array<double, n_points> euclidean_time;

  Lattice(double lattice_spacing, AnharmonicPotential const& potential)
      : m_lattice_spacing{lattice_spacing}, m_potential{potential}
  {
    for (size_t i{0}; i != n_points; ++i) {
#if defined COLD_START
      positions[i] = COLD_START;
#elif defined HOT_START
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> uniform(-HOT_START, HOT_START);
      positions[i] = uniform(gen);
#endif
      euclidean_time[i] = i * m_lattice_spacing;
    }
  }

  double calculate_action(size_t position_index)
  {
    auto x_pm = positions[position_index] - positions[position_index - 1];
    auto x_pp = positions[position_index + 1] - positions[position_index];
    return (1. / (4. * m_lattice_spacing)) * (x_pm * x_pm + x_pp * x_pp) +
           m_lattice_spacing * m_potential(positions[position_index]);
  }
};

#endif