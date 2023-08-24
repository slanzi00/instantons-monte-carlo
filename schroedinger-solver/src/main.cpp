#include <format>
#include <iostream>

#include "solver_kernels.hpp"
#include "timer.hpp"

int main()
{
  solver::InitialConditions ic(100, -3., 3.);
  Eigen::MatrixXd hamiltonian = Eigen::MatrixXd::Zero(ic.n_points, ic.n_points);
  Eigen::VectorXd dipole_matrix_elements = Eigen::VectorXd::Zero(ic.n_points);
  Eigen::VectorXd quadrupole_matrix_elements = Eigen::VectorXd::Zero(ic.n_points);
  Eigen::VectorXd exapole_matrix_elements = Eigen::VectorXd::Zero(ic.n_points);

#ifdef ENABLE_TIME_MEASUREMENT
  {
    solver::Timer timer;
    solver::serial::fill_hamiltonian(hamiltonian, solver::AnharmonicPotential(1.4), ic);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);
    Eigen::VectorXd energy_eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd wavefunctions = solver.eigenvectors();
    std::cout << energy_eigenvalues << '\n';
    // std::cout << hamiltonian << '\n';
    // solver::serial::calculate_matrix_elements(dipole_matrix_elements, wavefunctions, 1, ic);
    // solver::serial::calculate_matrix_elements(quadrupole_matrix_elements, wavefunctions, 2, ic);
    // solver::serial::calculate_matrix_elements(exapole_matrix_elements, wavefunctions, 3, ic);
    // for (int i{0}; i != ic.n_points; ++i) {
    //   std::cout << dipole_matrix_elements(i) << '\t' << quadrupole_matrix_elements(i) << '\t'
    //             << exapole_matrix_elements(i) << '\n';
    // }
  }
#else
  solver::serial::fill_hamiltonian(hamiltonian, anharmonic_potential, ic);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);
  Eigen::VectorXd eigenvalues = solver.eigenvalues();
  Eigen::MatrixXd eigenvectors = solver.eigenvectors();
  solver::serial::calculate_matrix_elements(dipole_matrix_elements, wavefunctions, 1, ic);
  solver::serial::calculate_matrix_elements(quadrupole_matrix_elements, wavefunctions, 2, ic);
  solver::serial::calculate_matrix_elements(exapole_matrix_elements, wavefunctions, 3, ic);
#endif
}