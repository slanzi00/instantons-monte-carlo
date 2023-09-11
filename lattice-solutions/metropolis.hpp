#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

#include <iostream>
#include <memory>
#include <vector>

#include <boost/histogram.hpp>

#include "correlators.hpp"
#include "lattice.hpp"

template <class Histogram, size_t lattice_points, size_t correlators_points>
class Metropolis
{
  std::shared_ptr<Lattice<lattice_points>> m_lattice;
  std::shared_ptr<Correlators<correlators_points>> m_correlators;
  Histogram m_histogram;
  std::mt19937 m_gen{std::random_device{}()};

 public:
  Metropolis(std::shared_ptr<Lattice<lattice_points>> lattice,
             std::shared_ptr<Correlators<correlators_points>> correlators,
             Histogram histogram)
      : m_lattice{std::move(lattice)}
      , m_correlators{std::move(correlators)}
      , m_histogram{std::move(histogram)}
  {
  }

  void one_sweep_monte_carlo()
  {
  }

  template <uint32_t n_sweeps>
  void evolve_lattice()
  {
    std::normal_distribution gaussian_step{0., 0.5};
    std::uniform_real_distribution probability{0., 1.};
    for (uint32_t i_sweep = 0; i_sweep != 1000; ++i_sweep) {
      for (int i = 1; i != lattice_points - 1; ++i) {
        auto initial_action = m_lattice->calculate_action(i);
        auto dx = gaussian_step(m_gen);
        m_lattice->positions[i] += dx;
        auto final_action = m_lattice->calculate_action(i);
        auto acceptance_probability = std::exp(initial_action - final_action);
        if (acceptance_probability <= 1. && acceptance_probability <= probability(m_gen)) {
          m_lattice->positions[i] -= dx;
        }
      }
    }

    for (uint32_t i_sweep = 0; i_sweep != n_sweeps; ++i_sweep) {
      for (int i = 1; i != lattice_points - 1; ++i) {
        auto initial_action = m_lattice->calculate_action(i);
        auto dx = gaussian_step(m_gen);
        m_lattice->positions[i] += dx;
        auto final_action = m_lattice->calculate_action(i);
        auto acceptance_probability = std::exp(initial_action - final_action);
        if (acceptance_probability <= 1. && acceptance_probability <= probability(m_gen)) {
          m_lattice->positions[i] -= dx;
        }
      }
      std::for_each(
          m_lattice->positions.begin(), m_lattice->positions.end(), std::ref(m_histogram));
      m_correlators->fill_correlators(m_lattice->positions, m_gen, 20);
    }
  }

  decltype(auto) probability_histogram()
  {
    // auto h = boost::histogram::make_histogram(
    //     boost::histogram::axis::regular<>(num_bins, x_min, x_max, "x"));
    // std::for_each(m_lattice->positions.begin(), m_lattice->positions.end(), std::ref(h));
    return m_histogram;
  }
};

#endif