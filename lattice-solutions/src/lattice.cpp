#include "lattice.hpp"

double AnharmonicPotential::operator()(double x)
{
  return std::pow(x * x - eta * eta, 2);
}

Lattice::Lattice(AnharmonicPotential const& potential) : m_potential{potential}
{
  for (size_t i = 0; i != sv::n_lattice_points; ++i) {
#if defined COLD_START
    positions[i] = COLD_START;
#elif defined HOT_START
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform(-HOT_START, HOT_START);
    positions[i] = uniform(gen);
#endif
    euclidean_time[i] = i * sv::lattice_spacing;
  }
}

double Lattice::calculate_action(size_t position_index)
{
  auto prev_index = (position_index + sv::n_lattice_points - 1) % (sv::n_lattice_points);
  auto next_index = (position_index + 1) % (sv::n_lattice_points);
  auto x_pm = positions[position_index] - positions[prev_index];
  auto x_pp = positions[next_index] - positions[position_index];
  return (1. / (4. * sv::lattice_spacing)) * (x_pm * x_pm + x_pp * x_pp) +
         sv::lattice_spacing * m_potential(positions[position_index]);
}

double Lattice::calculate_complete_action()
{
  auto result = 0.;
  for (auto position_idx = 0; position_idx != sv::n_lattice_points; ++position_idx) {
    auto next_index = (position_idx + 1) % (sv::n_lattice_points);
    // auto prev_index = (position_index )
    auto x_pp = positions[next_index] - positions[position_idx];
    result += (1. / (4. * sv::lattice_spacing)) * (x_pp * x_pp) +
              sv::lattice_spacing * m_potential(positions[position_idx]);
  }
  return result;
}