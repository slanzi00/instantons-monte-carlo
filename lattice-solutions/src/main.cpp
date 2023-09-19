#define HOT_START 1.4

#include <iostream>

#include "lattice.hpp"
#include "metropolis.hpp"

int main()
{
  using namespace boost::histogram;
  AnharmonicPotential potential{HOT_START};
  auto lattice = std::make_shared<Lattice<800>>(0.05, potential);
  auto correlators = std::make_shared<Correlators<30>>();
  auto correlators_cool = std::make_shared<Correlators<30>>();
  correlators->normalization = 1000000. * 200.;
  correlators_cool->normalization = 100000. * 40;
  Metropolis metropolis_evolver{
      lattice, correlators, correlators_cool, make_histogram(axis::regular<>(500, -3., 3., "x"))};

  metropolis_evolver.evolve_lattice<100000, 200>();
  // for (size_t j = 0; j != 29; ++j) {
  //   std::cout << lattice->euclidean_time[j] << '\t' << correlators->log_derivative(1, j, 0.05)
  //             << '\t' << correlators->log_derivative_error(1, j, 0.05) << '\n';
  //   // std::cout << lattice->euclidean_time[j] << '\t' << correlators->correlators(1, j) << '\n';
  //   // std::cout << lattice->euclidean_time[j] << '\t' << correlators->correlators(2, j) << '\n';
  // }

  // print_log_der(correlators->correlators);

  // auto h = metropolis_evolver.probability_histogram();
  // Correlators<30> pos_correlators;
  // pos_correlators(lattice->positions);

  // for (auto x : h) {
  //   std::cout << x << '\n';
  // }
  
  auto id = metropolis_evolver.instanton_density();
  for (int i{0}; i != 200; ++i) {
    std::cout << i << '\t' << id.actions[i] / id.n_instantons[i] << '\n';
  }
}
