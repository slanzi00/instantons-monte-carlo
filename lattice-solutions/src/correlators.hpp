#ifndef CORRELATORS_HPP
#define CORRELATORS_HPP

#include <cmath>
#include <random>

#include <boost/numeric/ublas/matrix.hpp>

#include "setting_values.hpp"

using correlators_t = boost::numeric::ublas::matrix<double>;
using zeros = boost::numeric::ublas::zero_matrix<double>;

struct Correlators
{
  // matricies with as many rows as correlators pows:
  // <x(0)x(tau)>, <x2(0)x2(tau)>, <x3(0)x3(tau)>
  correlators_t correlators;

  // storage also the squared values in order to calculate errors
  correlators_t correlators_squared;

  // a matrix to strage errors
  correlators_t correlators_errors;

  // public member functions
  Correlators();
  void fill_correlators(std::array<double, sv::n_lattice_points> const& positions,
                        std::mt19937 random_generator);
  void normalize_correlators();
  void calculate_errors();
  double log_derivative(size_t correlator_row, size_t correlator_col);
  double log_derivative_error(size_t correlator_row, size_t correlator_col);
};

#endif