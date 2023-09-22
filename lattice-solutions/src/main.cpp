#define HOT_START 1.4

#include <format>
#include <fstream>
#include <iostream>

#include "metropolis.hpp"

void print_histogram_csv(histogram const& histogram)
{
  std::ofstream histo_f("data/probability_histogram.csv");
  histo_f << "probability\n";
  for (auto i : histogram) {
    histo_f << std::format("{}\n",(double)i);
  }
}

void print_correlators_csv(std::shared_ptr<Lattice> lattice,
                           std::shared_ptr<Correlators> correlators)
{
  std::ofstream corr_f("data/correlators_log_derivative.csv");
  corr_f << "time,corr_1,dcorr_1,corr_2,d_corr2,corr_3,dcorr_3,log_der_c1,dlog_der_c1,log_der_c2,"
            "dlog_der_c2,log_der_c3,dlog_der_c3\n";
  for (size_t i = 0; i != sv::n_correlator_points - 1; ++i) {
    corr_f << std::format(
        "{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},"
        "{:4.12f},{:4.12f},{:4.12f}\n",
        lattice->euclidean_time[i],
        correlators->correlators(0, i),
        correlators->correlators_errors(0, i),
        correlators->correlators(1, i),
        correlators->correlators_errors(1, i),
        correlators->correlators(2, i),
        correlators->correlators_errors(2, i),
        correlators->log_derivative(0, i),
        correlators->log_derivative_error(0, i),
        correlators->log_derivative(1, i),
        correlators->log_derivative_error(1, i),
        correlators->log_derivative(2, i),
        correlators->log_derivative_error(2, i));
  }
}

void print_correlators_cool_csv(std::shared_ptr<Lattice> lattice,
                                std::shared_ptr<Correlators> correlators_cool)
{
  std::ofstream corr_cool_f("data/correlators_log_derivative_cool.csv");
  corr_cool_f
      << "time,corr_1,dcorr_1,corr_2,d_corr2,corr_3,dcorr_3,log_der_c1,dlog_der_c1,log_der_c2,"
         "dlog_der_c2,log_der_c3,dlog_der_c3\n";
  for (size_t i = 0; i != sv::n_correlator_points - 1; ++i) {
    corr_cool_f << std::format(
        "{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},"
        "{:4.12f},{:4.12f},{:4.12f}\n",
        lattice->euclidean_time[i],
        correlators_cool->correlators(0, i),
        correlators_cool->correlators_errors(0, i),
        correlators_cool->correlators(1, i),
        correlators_cool->correlators_errors(1, i),
        correlators_cool->correlators(2, i),
        correlators_cool->correlators_errors(2, i),
        correlators_cool->log_derivative(0, i),
        correlators_cool->log_derivative_error(0, i),
        correlators_cool->log_derivative(1, i),
        correlators_cool->log_derivative_error(1, i),
        correlators_cool->log_derivative(2, i),
        correlators_cool->log_derivative_error(2, i));
  }
}

void print_instanton_density_csv(std::array<uint32_t, sv::n_sweeps_cool> const& instantons_density,
                                 std::array<double, sv::n_sweeps_cool> const& actions)
{
  std::ofstream id_f("data/instanton_density_16.csv");
  id_f << "n_cool,n_inst,action\n";
  for (size_t i = 0; i != sv::n_sweeps_cool; ++i) {
    id_f << std::format("{},{},{:6.1f}\n", i, instantons_density[i], actions[i]);
  }
}

int main()
{
  using namespace boost::histogram;
  AnharmonicPotential potential{HOT_START};
  auto lattice = std::make_shared<Lattice>(potential);
  auto correlators = std::make_shared<Correlators>();
  auto correlators_cool = std::make_shared<Correlators>();
  Metropolis metropolis_evolver{
      lattice,
      correlators,
      correlators_cool,
      make_histogram(
          axis::regular<>(sv::n_histogram_bins, sv::x_min_histogram, sv::x_max_histogram, "x"))};

  metropolis_evolver.evolve_lattice();
  print_histogram_csv(metropolis_evolver.probability_histogram());
  print_correlators_csv(lattice, correlators);
  print_correlators_cool_csv(lattice, correlators_cool);
  print_instanton_density_csv(metropolis_evolver.number_instantons(), metropolis_evolver.actions());
}
