#include "solver_kernels.hpp"

namespace solver {

namespace serial {

void fill_hamiltonian(Eigen::MatrixXd& hamiltonian,
                      Potential const& potential,
                      InitialConditions& ic)
{
  assert(hamiltonian.rows() == hamiltonian.cols());

  Eigen::MatrixXd kinetic_A = Eigen::MatrixXd::Zero(ic.n_points, ic.n_points);
  kinetic_A.diagonal(0).setConstant(2. / (2. * std::pow(ic.get_step(), 2)));
  kinetic_A.diagonal(1).setConstant(-1. / (2. * std::pow(ic.get_step(), 2)));
  kinetic_A.diagonal(-1).setConstant(-1. / (2. * std::pow(ic.get_step(), 2)));

  Eigen::MatrixXd kinetic_B = Eigen::MatrixXd::Zero(ic.n_points, ic.n_points);
  kinetic_B.diagonal(0).setConstant(-5. / 12.);
  kinetic_B.diagonal(1).setConstant(-1. / 24.);
  kinetic_B.diagonal(-1).setConstant(-1. / 24.);

  Eigen::MatrixXd potential_V = Eigen::MatrixXd::Zero(ic.n_points, ic.n_points);
  for (int i{}; i != ic.n_points; ++i) {
    potential_V(i, i) = potential(ic.x_min + i * ic.get_step());
  }

  hamiltonian = kinetic_B.inverse() * kinetic_A + potential_V;
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