#include "solver_init_conditions.hpp"

namespace solver {

double InitialConditions::get_step()
{
  return (x_max - x_min) / (n_points - 1.);
}

double InitialConditions::get_kinetic()
{
  return 1. / (2. * this->get_step() * this->get_step());
}

};  // namespace solver