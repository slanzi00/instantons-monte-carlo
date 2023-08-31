#include <iostream>

#include "initial_conditions.hpp"
#include "potentials.hpp"
#include "solver.hpp"

int main()
{
  // anharmonic potential : x^4 - 2 * eta^2 * x^2 + eta^4
  double eta = 1.4;
  PolynomialPotential<4> anharmonic_potential;
  anharmonic_potential.coefficients << std::pow(eta, 4), 0., -2. * std::pow(eta, 2), 0., 1.;

  InitialConditions ic{anharmonic_potential, 10.};

  SchroedingerSolver solver{ic};
  auto positions = solver.get_positions();
  auto energies = solver.energy_eigenvalues();
  auto wavefunctions = solver.wavefunctions();
  auto dipole = solver.position_matrix_elements<1>(wavefunctions);
  auto quadrupole = solver.position_matrix_elements<2>(wavefunctions);
  auto exapole = solver.position_matrix_elements<3>(wavefunctions);
  auto correlator_1 = solver.correlator(dipole, energies);
  auto correlator_2 = solver.correlator(quadrupole, energies);
  auto correlator_3 = solver.correlator(exapole, energies);

  // std::cout << "x limits: x_min = " << ic.x_min() << ", x_max = " << ic.x_max() << '\n'
  //           << "# points: " << ic.n_points() << '\n'
  //           << "step: " << ic.get_step() << '\n';

  // for (int i = 0; i != ic.n_points(); ++i) {
  //   std::cout << positions(i) << ' ' << wavefunctions(i, 0) << '\n';
  // }

  std::cout << correlator_1 << '\n';
}