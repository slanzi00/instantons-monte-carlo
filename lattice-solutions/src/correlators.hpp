#ifndef CORRELATORS_HPP
#define CORRELATORS_HPP

#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>

using correlators_t = boost::numeric::ublas::matrix<double>;
using zeros = boost::numeric::ublas::zero_matrix<double>;

template <size_t n_points>
struct Correlators
{
  double normalization;
  correlators_t correlators;
  correlators_t correlators_squared;
  correlators_t correlators_errors;

  Correlators()
      : normalization{}
      , correlators{zeros(3, n_points)}
      , correlators_squared{zeros(3, n_points)}
      , correlators_errors{zeros(3, n_points)}
  {
  }

  void fill_correlators(auto const& positions, auto& random_generator, uint16_t n_measurement)
  {
    std::uniform_int_distribution<size_t> random_idx_generator{0, positions.size() - n_points - 1};
    for (uint16_t j = 0; j != n_measurement; ++j) {
      auto random_index = random_idx_generator(random_generator);
      for (size_t i = 0; i != n_points; ++i) {
        auto correlator_i = positions[random_index] * positions[random_index + i];
        correlators(0, i) += correlator_i;
        correlators_squared(0, i) += std::pow(correlator_i, 2);
        correlators(1, i) += std::pow(correlator_i, 2);
        correlators_squared(1, i) += std::pow(correlator_i, 4);
        correlators(2, i) += std::pow(correlator_i, 3);
        correlators_squared(2, i) += std::pow(correlator_i, 6);
      }
    }
  }

  void normalize_correlators()
  {
    for (size_t i = 0; i != n_points; ++i) {
      for (size_t j : {0, 1, 2}) {
        correlators(j, i) /= normalization;
        correlators_squared(j, i) /= normalization;
      }
    }
  }

  void calculate_errors()
  {
    for (size_t i = 0; i != n_points; ++i) {
      for (size_t j : {0, 1, 2}) {
        auto square_mean = correlators_squared(j, i);
        auto mean_square = std::pow(correlators(j, i), 2);
        correlators_errors(j, i) = std::sqrt((square_mean - mean_square) / normalization);
      }
    }
  }

  double log_derivative(auto correlator_row, auto correlator_col, auto lattice_spacing)
  {
    auto c = correlators(correlator_row, correlator_col);
    auto c_next = correlators(correlator_row, correlator_col + 1);
    auto c_2_last = correlators(1, n_points - 1);
    if (correlator_row != 1) {
      return (c - c_next) / (c * lattice_spacing);
    } else {
      return (c - c_next) / ((c - c_2_last) * lattice_spacing);
    }
  }

  double log_derivative_error(auto correlator_row, auto correlator_col, auto lattice_spacing)
  {
    double c = correlators(correlator_row, correlator_col);
    double c_next = correlators(correlator_row, correlator_col + 1);
    double c_2_sub = correlators(correlator_row, correlator_col) - correlators(1, n_points - 1);
    double c_2_sub_next =
        correlators(correlator_row, correlator_col + 1) - correlators(1, n_points - 1);

    double c_err = correlators_errors(correlator_row, correlator_col);
    double c_err_next = correlators_errors(correlator_row, correlator_col + 1);
    double c_2_err_last = correlators_errors(1, n_points - 1);
    double c_2_err_sub = std::sqrt(std::pow(c_err, 2) + std::pow(c_2_err_last, 2));
    double c_2_err_next_sub = std::sqrt(std::pow(c_err_next, 2) + std::pow(c_2_err_last, 2));

    if (correlator_row != 1) {
      double error_square = std::pow(c_err_next / c, 2) + std::pow((c_err * c_next) / (c * c), 2);
      return std::sqrt(error_square) / lattice_spacing;
    } else {
      double error_square = std::pow(c_2_err_next_sub / c_2_sub, 2) +
                            std::pow((c_2_err_sub * c_2_sub_next) / (c_2_sub * c_2_sub), 2);
      return std::sqrt(error_square) / lattice_spacing;
    }
  }
};

#endif