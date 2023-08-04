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

void create_hamiltonian(Eigen::MatrixXd& hamiltonian,
                        Potential const& potential,
                        InitialConditions& ic);

};  // namespace serial

};  // namespace solver

#endif