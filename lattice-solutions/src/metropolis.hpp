#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

#include <iostream>
#include <memory>

#include <boost/histogram.hpp>

#include "correlators.hpp"
#include "lattice.hpp"

// choose the evolution type
enum class EvolutionType { Normal, Cooling };

// verify evolution conditions
bool evolution_condition_verified(EvolutionType, double, double);

using histogram = boost::histogram::histogram<
    std::tuple<boost::histogram::axis::
                   regular<double, boost::use_default, boost::use_default, boost::use_default>>,
    boost::histogram::unlimited_storage<>>;

class Metropolis
{
  std::shared_ptr<Lattice> m_lattice;
  std::shared_ptr<Correlators> m_correlators;
  std::shared_ptr<Correlators> m_correlators_cool;
  std::array<double, sv::n_lattice_points> m_positions_storage;
  std::array<uint32_t, sv::n_sweeps_cool> m_n_instantons;
  std::array<double, sv::n_sweeps_cool> m_actions;
  std::mt19937 m_gen;
  histogram m_histogram;

 public:
  Metropolis(std::shared_ptr<Lattice> lattice,
             std::shared_ptr<Correlators> correlators,
             std::shared_ptr<Correlators> correlators_cool,
             histogram histogram);

  void one_sweep_monte_carlo(EvolutionType evolution_t,
                             std::normal_distribution<double>& gauss_step_dist,
                             std::uniform_real_distribution<double>& random_dist,
                             bool fill_histo_enabled = true);

  int calculate_number_instantons();

  void evolve_lattice(bool make_cooling = true);

  histogram probability_histogram();

  std::array<uint32_t, sv::n_sweeps_cool> number_instantons();

  std::array<double, sv::n_sweeps_cool> actions();
};

inline bool evolution_condition_verified(EvolutionType evolution_type,
                                         double acceptance_probability,
                                         double random_number)
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

#endif