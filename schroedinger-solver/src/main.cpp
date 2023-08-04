#include <format>
#include <iostream>

#include "solver_kernels.hpp"
#include "timer.hpp"

int main()
{
  solver::InitialConditions ic(1000, -4., 4.);
  solver::AnharmonicPotential anharmonic_potential(1.4);
  Eigen::MatrixXd hamiltonian = Eigen::MatrixXd::Zero(ic.n_points, ic.n_points);

#ifdef ENABLE_TIME_MEASUREMENT
  {
    solver::Timer timer;
    solver::serial::create_hamiltonian(hamiltonian, anharmonic_potential, ic);
  }
#else
  solver::serial::create_hamiltonian(hamiltonian, anharmonic_potential, ic);
#endif

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);
  Eigen::VectorXd eigenvalues = solver.eigenvalues();
  Eigen::MatrixXd eigenvectors = solver.eigenvectors();
}