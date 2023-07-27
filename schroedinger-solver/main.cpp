#include <cmath>
#include <iostream>

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <tbb/tbb.h>

double harmonic_potential(double x)
{
  return 0.5f * std::pow(x, 2);
}

int main()
{
  constexpr int n_points = 10;  // dimension
  constexpr double x_min = -5.f;
  constexpr double x_max = 5.f;
  constexpr double h = (x_max - x_min) / n_points;    // step
  constexpr double k = 1.f / (2.f * std::pow(h, 2));  // wave number

  Eigen::SparseMatrix<double> hamiltonian(n_points, n_points);

  // fill matrix in parallel on cpu
  tbb::parallel_for(0, n_points, [&](int i) {
    hamiltonian.insert(i, i) = 2 * k + harmonic_potential(x_min + (i + 1) * h);
    if (i > 0) {
      hamiltonian.insert(i - 1, i) = -k;
    }
    if (i < n_points - 1) {
      hamiltonian.insert(i + 1, i) = -k;
    }
  });

  // make sparse the hamiltonian
  hamiltonian.makeCompressed();

  std::cout << hamiltonian << '\n';

  // find eigenvalues (energy levels) and eigenvactors (wavefunction)
  Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<double>> solver(hamiltonian);
  Eigen::VectorXd energy_eigenvalues = solver.eigenvalues();
  Eigen::MatrixXd wavefunction = solver.eigenvectors();
  Eigen::SparseMatrix<double> wavefunction_sparse = wavefunction.sparseView();

  std::cout << energy_eigenvalues << '\n';

  std::cout << wavefunction_sparse << '\n';
}
