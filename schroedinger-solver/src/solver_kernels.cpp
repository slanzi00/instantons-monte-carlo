#include "solver_kernels.hpp"

namespace solver {

namespace serial {

void create_hamiltonian(Eigen::MatrixXd& hamiltonian,
                        Potential const& potential,
                        InitialConditions& ic)
{
  assert(hamiltonian.rows() == hamiltonian.cols());
  for (int i{0}; i != ic.n_points; ++i) {
    hamiltonian(i, i) = 2 * ic.get_kinetic() + potential(ic.x_min + (i + 1) * ic.get_step());
    if (i > 0) {
      hamiltonian(i - 1, i) = -ic.get_kinetic();
    }
    if (i < ic.n_points - 1) {
      hamiltonian(i + 1, i) = -ic.get_kinetic();
    }
  }
}

};  // namespace serial

};  // namespace solver