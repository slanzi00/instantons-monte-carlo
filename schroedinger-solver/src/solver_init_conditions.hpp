#ifndef SOLVER_INIT_CONDITIONS_HPP
#define SOLVER_INIT_CONDITIONS_HPP

namespace solver {

struct InitialConditions
{
  int n_points;
  double x_min;
  double x_max;
  double get_step();
  double get_kinetic();
};

};  // namespace solver
#endif