#include <format>
#include <fstream>
#include <iostream>

#include "initial_conditions.hpp"
#include "potentials.hpp"
#include "solver.hpp"

void print_eigens_csv(auto const& positions, auto const& eigenvalues, auto const& wavefuncions)
{
  std::ofstream eigens_f("data/eigenvalues_eigenfunctions.csv");
  eigens_f << "positions,eigenvalues,eigenfuntion(0),eigenfuncion(1),eigenfunction(2)\n";
  for (auto i = 0; i != eigenvalues.size(); ++i) {
    eigens_f << std::format("{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f}\n",
                            positions[i],
                            eigenvalues[i],
                            -wavefuncions(i, 0),
                            wavefuncions(i, 1),
                            wavefuncions(i, 2));
  }
}

void print_correlators_csv(auto& solver,
                           auto const& eigenvalues,
                           auto const& dipole,
                           auto const& quadrupole,
                           auto const& exapole)
{
  std::ofstream corr_f("data/correlators_log_derivative.csv");
  corr_f << "time,<x(τ)x(0)>,<x(τ)²x(0)²>,<x(τ)³x(0)³>,dlog<x(τ)x(0)>/dτ,dlog<x(τ)²x(0)²>/"
            "dτ,dlog<x(τ)³x(0)³>/dτ\n";
  size_t n_tau = 200;
  double tau_max = 1.5;
  double tau = 0.;
  for (size_t i = 0; i != n_tau; ++i) {
    corr_f << std::format("{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f},{:4.12f}\n",
                          tau,
                          solver.correlator(dipole, eigenvalues, tau),
                          solver.correlator(quadrupole, eigenvalues, tau),
                          solver.correlator(exapole, eigenvalues, tau),
                          solver.log_derivative(dipole, eigenvalues, tau),
                          solver.log_derivative(quadrupole, eigenvalues, tau, true),
                          solver.log_derivative(exapole, eigenvalues, tau));
    tau += tau_max / static_cast<double>(n_tau);
  }
}

int main()
{
  // anharmonic potential : x^4 - 2 * eta^2 * x^2 + eta^4
  double eta = 1.4;
  PolynomialPotential<4> anharmonic_potential;
  anharmonic_potential.coefficients << std::pow(eta, 4), 0., -2. * std::pow(eta, 2), 0., 1.;

  InitialConditions ic{anharmonic_potential, 400.};

  SchroedingerSolver solver{ic};
  auto positions = solver.get_positions();
  auto energies = solver.energy_eigenvalues();
  auto wavefunctions = solver.wavefunctions();
  auto dipole = solver.position_matrix_elements<1>(wavefunctions);
  auto quadrupole = solver.position_matrix_elements<2>(wavefunctions);
  auto exapole = solver.position_matrix_elements<3>(wavefunctions);

  print_eigens_csv(positions, energies, wavefunctions);
  print_correlators_csv(solver, energies, dipole, quadrupole, exapole);
}