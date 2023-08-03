#include <Eigen/Dense>  // Assicurati di avere Eigen installato: https://eigen.tuxfamily.org
#include <cmath>
#include <iostream>
#include <vector>

using namespace Eigen;

const double h_bar = 1.0;  // Costante di Planck ridotta
const double m = 1.0;      // Massa della particella

// Potenziale armonico V(x) = 0.5 * x^2
double potential(double x)
{
  return 0.5 * std::pow(x, 2);
}

// Costruzione dell'operatore di derivata seconda (derivata centrale) con le differenze finite
MatrixXd second_derivative_operator(int N, double dx)
{
  MatrixXd operator_matrix = MatrixXd::Zero(N, N);
  for (int i = 1; i < N - 1; ++i) {
    operator_matrix(i, i) = -2.0;
    operator_matrix(i, i - 1) = operator_matrix(i, i + 1) = 1.0;
  }
  operator_matrix(0, 0) = operator_matrix(N - 1, N - 1) = -2.0;
  operator_matrix(0, 1) = operator_matrix(N - 1, N - 2) = 1.0;
  operator_matrix /= (dx * dx);
  return operator_matrix;
}

// Costruzione dell'operatore hamiltoniano H = -h_bar^2 / (2m) * d^2/dx^2 + V(x)
MatrixXd hamiltonian(const VectorXd& V, double h_bar, double m, int N, double dx)
{
  MatrixXd H = (-h_bar * h_bar / (2 * m)) * second_derivative_operator(N, dx);
  for (int i = 0; i < N; ++i) {
    H(i, i) += V(i);
  }
  return H;
}

// Risoluzione dell'equazione di Schrödinger stazionaria H(psi(x)) = E * psi(x)
VectorXd solve_stationary_schrodinger(const VectorXd& V, double h_bar, double m, int N, double dx)
{
  MatrixXd H = hamiltonian(V, h_bar, m, N, dx);

  SelfAdjointEigenSolver<MatrixXd> solver(H);
  VectorXd eigenvalues = solver.eigenvalues();
  MatrixXd eigenvectors = solver.eigenvectors();

  return eigenvectors.col(
      0);  // Restituisce il primo autovettore corrispondente all'energia più bassa
}

int main()
{
  double x_min = -1.;
  double x_max = 1;
  int N = 1000;
  double dx = (x_max - x_min) / (N - 1);

  VectorXd x(N);
  for (int i = 0; i < N; ++i) {
    x(i) = x_min + i * dx;
  }

  VectorXd V(N);
  for (int i = 0; i < N; ++i) {
    V(i) = potential(x(i));
  }

  VectorXd psi = solve_stationary_schrodinger(V, h_bar, m, N, dx);

  // Stampiamo gli autovettori su righe separate
  std::cout << "#Autovettori:" << std::endl;
  auto sum = 0.;
  for (int i = 0; i < N; ++i) {
    // std::cout << x_min + i * dx << ' ' << † << std::endl;
    sum += psi(i) * psi(i) * dx;
  }

  std::cout << sum << '\n';
}
