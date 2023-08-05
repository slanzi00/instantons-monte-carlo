#ifndef SOLVER_KERNELS_TBB_HPP
#define SOLVER_KERNELS_TBB_HPP

#include <chrono>
#include <iostream>
#include <stdexcept>

#include <Eigen/Dense>

#include "potentials.hpp"
#include "solver_init_conditions.hpp"

namespace solver {

namespace serial {

void fill_hamiltonian(Eigen::MatrixXd& hamiltonian,
                      Potential const& potential,
                      InitialConditions& ic);

void calculate_matrix_elements(Eigen::VectorXd& n_pole,
                               Eigen::MatrixXd const& wavefunctions,
                               uint8_t position_pow,
                               InitialConditions& ic);

};  // namespace serial

};  // namespace solver

#endif