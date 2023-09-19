#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

#include <iostream>
#include <memory>
#include <vector>

#include <boost/histogram.hpp>

#include "correlators.hpp"
#include "lattice.hpp"

template <size_t lattice_points>
using lattice_ptr = std::shared_ptr<Lattice<lattice_points>>;

template <size_t correlators_points>
using correlators_ptr = std::shared_ptr<Correlators<correlators_points>>;

enum class EvolutionType { Normal, Cooling };

bool evolution_condition_verified(auto evolution_type,
                                  auto acceptance_probability,
                                  auto random_number)
{
  switch (evolution_type) {
    case EvolutionType::Normal:
      return acceptance_probability <= 1. && acceptance_probability <= random_number;
      break;
    case EvolutionType::Cooling:
      return acceptance_probability <= 1.;
      break;
    default:
      std::cerr << "Evolution condition not specified.\n";
      return false;
  }
}

struct InstantonDensity
{
  std::array<double, 1000> n_instantons;
  std::array<double, 1000> actions;
};

template <class Histogram, size_t lattice_points, size_t correlators_points>
class Metropolis
{
  lattice_ptr<lattice_points> m_lattice;
  correlators_ptr<correlators_points> m_correlators;
  correlators_ptr<correlators_points> m_correlators_cool;
  std::array<double, lattice_points> m_positions_storage;
  InstantonDensity m_instanton_density;
  std::mt19937 m_gen;
  Histogram m_histogram;

 public:
  Metropolis(lattice_ptr<lattice_points> lattice,
             correlators_ptr<correlators_points> correlators,
             correlators_ptr<correlators_points> correlators_cool,
             Histogram histogram)
      : m_lattice{std::move(lattice)}
      , m_correlators{std::move(correlators)}
      , m_correlators_cool{std::move(correlators_cool)}
      , m_positions_storage{}
      , m_gen{std::random_device{}()}
      , m_histogram{std::move(histogram)}
  {
  }

  void one_sweep_monte_carlo(auto evolution_t,
                             auto& gauss_step_dist,
                             auto& random_dist,
                             bool fill_histo_enabled = true)
  {
    for (int i = 0; i != lattice_points; ++i) {
      auto initial_action = m_lattice->calculate_action(i);
      auto dx = gauss_step_dist(m_gen);
      m_lattice->positions[i] += dx;
      auto final_action = m_lattice->calculate_action(i);
      auto acceptance_probability = std::exp(initial_action - final_action);
      if (evolution_condition_verified(evolution_t, acceptance_probability, random_dist(m_gen))) {
        m_lattice->positions[i] -= dx;
      }
    }
    if (fill_histo_enabled) {
      std::for_each(
          m_lattice->positions.begin(), m_lattice->positions.end(), std::ref(m_histogram));
    }
  }

  int calculate_number_instantons()
  {
    auto result = 0;
    for (size_t i = 0; i != m_lattice->positions.size(); ++i) {
      auto prev_index = (i + m_lattice->positions.size() - 1) % (m_lattice->positions.size());
      if (std::signbit(m_lattice->positions[prev_index]) != std::signbit(m_lattice->positions[i])) {
        ++result;
      }
    }
    return result;
  }

  template <uint64_t n_sweeps, uint32_t n_cooling_sweeps>
  void evolve_lattice()
  {
    std::normal_distribution gaussian_step{0., 0.5};
    std::uniform_real_distribution probability{0., 1.};

    for (uint16_t thermalization = 0; thermalization != 100; ++thermalization) {
      one_sweep_monte_carlo(EvolutionType::Normal, gaussian_step, probability, false);
    }

    for (uint64_t i_sweep = 0; i_sweep != n_sweeps; ++i_sweep) {
      if (i_sweep % 2 == 0) {
        m_correlators->fill_correlators(m_lattice->positions, m_gen, 400);
      }
      if (i_sweep % 200 == 0) {
        m_positions_storage = m_lattice->positions;
        for (uint32_t i_cool = 0; i_cool != n_cooling_sweeps; ++i_cool) {
          one_sweep_monte_carlo(EvolutionType::Cooling, gaussian_step, probability);
          m_instanton_density.n_instantons[i_cool] += calculate_number_instantons();
          m_instanton_density.actions[i_cool] += m_lattice->calculate_complete_action();
          if (i_cool % 2 == 0) {
            m_correlators_cool->fill_correlators(m_lattice->positions, m_gen, 100);
          }
        }
        m_lattice->positions = m_positions_storage;
      }
      one_sweep_monte_carlo(EvolutionType::Normal, gaussian_step, probability);
    }
    m_correlators->normalize_correlators();
    m_correlators->calculate_errors();
    m_correlators_cool->normalize_correlators();
    m_correlators_cool->calculate_errors();
  }

  decltype(auto) probability_histogram()
  {
    return m_histogram;
  }

  decltype(auto) instanton_density()
  {
    return m_instanton_density;
  }
};

#endif