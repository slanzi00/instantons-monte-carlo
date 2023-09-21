#include "metropolis.hpp"

Metropolis::Metropolis(std::shared_ptr<Lattice> lattice,
                       std::shared_ptr<Correlators> correlators,
                       std::shared_ptr<Correlators> correlators_cool,
                       histogram histogram)
    : m_lattice{std::move(lattice)}
    , m_correlators{std::move(correlators)}
    , m_correlators_cool{std::move(correlators_cool)}
    , m_positions_storage{}
    , m_n_instantons{}
    , m_actions{}
    , m_gen{std::random_device{}()}
    , m_histogram{std::move(histogram)}
{
}

void Metropolis::one_sweep_monte_carlo(EvolutionType evolution_t,
                                       std::normal_distribution<double>& gauss_step_dist,
                                       std::uniform_real_distribution<double>& random_dist,
                                       bool fill_histo_enabled)
{
  for (int i = 0; i != sv::n_lattice_points; ++i) {
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
    std::for_each(m_lattice->positions.begin(), m_lattice->positions.end(), std::ref(m_histogram));
  }
}

int Metropolis::calculate_number_instantons()
{
  int result = 0;
  for (size_t i = 0; i != sv::n_lattice_points; ++i) {
    auto next_index = (i + 1) % (sv::n_lattice_points);
    if (std::signbit(m_lattice->positions[i]) != std::signbit(m_lattice->positions[next_index])) {
      ++result;
    }
  }
  return result;
}

void Metropolis::evolve_lattice(bool make_cooling)
{
  std::normal_distribution gaussian_step{0., 0.45};
  std::uniform_real_distribution probability{0., 1.};

  for (uint16_t th = 0; th != sv::n_therm_sweeps; ++th) {
    one_sweep_monte_carlo(EvolutionType::Normal, gaussian_step, probability, false);
  }

  for (uint64_t i_sweep = 0; i_sweep != sv::n_sweeps; ++i_sweep) {
    if (i_sweep % sv::take_corr_meas_every == 0) {
      m_correlators->fill_correlators(m_lattice->positions, m_gen);
    }
    if ((i_sweep % sv::make_cool_every == 0) && make_cooling == true) {
      m_positions_storage = m_lattice->positions;
      for (uint32_t i_cool = 0; i_cool != sv::n_sweeps_cool; ++i_cool) {
        m_n_instantons[i_cool] += calculate_number_instantons();
        m_actions[i_cool] += m_lattice->calculate_complete_action();
        one_sweep_monte_carlo(EvolutionType::Cooling, gaussian_step, probability);
      }
      // correlation functions after cooling
      m_correlators_cool->fill_correlators(m_lattice->positions, m_gen);
      m_lattice->positions = m_positions_storage;
    }
    one_sweep_monte_carlo(EvolutionType::Normal, gaussian_step, probability);
  }
  m_correlators->normalize_correlators();
  m_correlators->calculate_errors();
  m_correlators_cool->normalize_correlators_cool();
  m_correlators_cool->calculate_errors_cool();
}

histogram Metropolis::probability_histogram()
{
  double norm = 0.,
         bin_width = (sv::x_max_histogram - sv::x_min_histogram) / (double)sv::n_histogram_bins;

  for (auto bin : m_histogram) {
    norm += bin * bin_width;
  }

  for (auto bin : m_histogram) {
    bin /= norm;
  }
  return m_histogram;
}

std::array<uint32_t, sv::n_sweeps_cool> Metropolis::number_instantons()
{
  return m_n_instantons;
}

std::array<double, sv::n_sweeps_cool> Metropolis::actions()
{
  return m_actions;
}