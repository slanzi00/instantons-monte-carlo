#define HOT_START 1.4

#include <iostream>

#include "lattice.hpp"
#include "metropolis.hpp"

int main()
{
  using namespace boost::histogram;
  AnharmonicPotential potential{HOT_START};
  auto lattice = std::make_shared<Lattice<1000>>(0.05, potential);
  auto correlators = std::make_shared<Correlators<30>>();
  Metropolis metropolis_evolver{
      lattice, correlators, make_histogram(axis::regular<>(50, -2., 2., "x"))};

  metropolis_evolver.evolve_lattice<1000000>();
  // for (size_t i{0}; i != correlators->correlators.size1(); ++i) {
  //   for (size_t j{0}; j != correlators->correlators.size2(); ++j) {
  //     std::cout << lattice->euclidean_time[j] << '\t' << correlators->correlators(i, j) / 20000000.
  //               << '\n';
  //   }
  // }

  auto h = metropolis_evolver.probability_histogram();
  // Correlators<30> pos_correlators;
  // pos_correlators(lattice->positions);

  // for (auto x : h) {
  //   std::cout << x << '\n';
  // }

  // for (int i{0}; i != 800; ++i) {
  //   std::cout << lattice->euclidean_time[i] << '\t' << lattice->positions[i] << '\n';
  // }
}
