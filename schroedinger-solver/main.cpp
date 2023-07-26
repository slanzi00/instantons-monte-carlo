#include <iostream>
#include <cmath>
#include <Eigen/Dense>

auto main() -> int
{
  constexpr std::size_t N_POINTS = 3;
  Eigen::MatrixXd matrix(N_POINTS, N_POINTS);
  for (int i{}; i != N_POINTS; ++i) {
    for (auto j : {i - 1, i, i + 1}) {
      if (j != -1 && j != N_POINTS) {
        matrix(i, j) = 1.f;
      }
    }
  }
  std::cout << matrix << '\n';
  Eigen::EigenSolver<Eigen::MatrixXd> solver(matrix);
  Eigen::VectorXd eigenvalues = solver.eigenvalues().real();
  Eigen::MatrixXd eigenvectors = solver.eigenvectors().real();
  std::cout << "Eigenvalues:\n"
            << eigenvalues << std::endl;
  std::cout << "Eigenvectors:\n"
            << eigenvectors << std::endl;
}