#include "solver_kernels.hpp"

namespace solver {

namespace serial {

void fill_hamiltonian(Eigen::MatrixXd& hamiltonian,
                      Potential const& potential,
                      InitialConditions& ic)
{
  assert(hamiltonian.rows() == hamiltonian.cols());
  for (int i{0}; i != ic.n_points; ++i) {
    hamiltonian(i, i) = 2 * ic.get_kinetic() + potential(ic.x_min + i * ic.get_step());
    if (i > 0) {
      hamiltonian(i - 1, i) = -ic.get_kinetic();
    }
    if (i < ic.n_points - 1) {
      hamiltonian(i + 1, i) = -ic.get_kinetic();
    }
  }
}

void calculate_matrix_elements(Eigen::VectorXd& n_pole,
                               Eigen::MatrixXd const& wavefunctions,
                               uint8_t position_pow,
                               InitialConditions& ic)
{
  for (int n{0}; n != ic.n_points; ++n) {
    for (int i{0}; i != ic.n_points; ++i) {
      n_pole(n) += wavefunctions(i, 0) * std::pow(ic.x_min + i * ic.get_step(), position_pow) *
                   wavefunctions(i, n);
    }
  }
  n_pole = n_pole.array().pow(2);
}

};  // namespace serial

};  // namespace solver