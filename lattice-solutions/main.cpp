#define COLD_START 1.4

#include <iostream>

#include "lattice.hpp"
#include "metropolis.hpp"

void print_log_der(auto const& correlators)
{
  for (int i = 0; i != 29; ++i) {
    std::cout << i << ' ' << (correlators(0, i) - correlators(0, i + 1)) / correlators(0, i) / 0.05
              << '\n';
  }
  for (int i = 0; i != 29; ++i) {
    std::cout << i << ' ' << (correlators(2, i) - correlators(2, i + 1)) / correlators(2, i) / 0.05
              << '\n';
  }
}

int main()
{
  using namespace boost::histogram;
  AnharmonicPotential potential{COLD_START};
  auto lattice = std::make_shared<Lattice<800>>(0.05, potential);
  auto correlators = std::make_shared<Correlators<30>>();
  correlators->normalization = 1000000. * 200.;
  Metropolis metropolis_evolver{
      lattice, correlators, make_histogram(axis::regular<>(500, -3., 3., "x"))};

  metropolis_evolver.evolve_lattice<1000000, 200>();
  // correlators->normalize_correlators();
  // for (size_t j = 0; j != 29; ++j) {
  //   std::cout << lattice->euclidean_time[j] << '\t' << correlators->log_derivative(0, j, 0.05)
  //             << '\n';
  // }

  print_log_der(correlators->correlators);

  // auto h = metropolis_evolver.probability_histogram();
  // Correlators<30> pos_correlators;
  // pos_correlators(lattice->positions);

  // for (auto x : h) {
  //   std::cout << x << '\n';
  // }

  // for (int i{0}; i != 800; ++i) {
  //   std::cout << lattice->euclidean_time[i] << '\t' << lattice->positions[i] << '\n';
  // }
}
