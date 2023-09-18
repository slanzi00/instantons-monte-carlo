#ifndef CORRELATORS_HPP
#define CORRELATORS_HPP

#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>

template <size_t n_points>
struct Correlators
{
  double normalization;
  boost::numeric::ublas::matrix<double> correlators;
  boost::numeric::ublas::matrix<double> correlators_squared;
  boost::numeric::ublas::matrix<double> correlators_errors;

  Correlators()
      : normalization{}
      , correlators{boost::numeric::ublas::zero_matrix<double>(3, n_points)}
      , correlators_squared{boost::numeric::ublas::zero_matrix<double>(3, n_points)}
      , correlators_errors{boost::numeric::ublas::zero_matrix<double>(3, n_points)}
  {
  }

  void fill_correlators(auto const& positions, auto& random_generator, uint16_t n_measurement)
  {
    std::uniform_int_distribution<size_t> random_idx_generator{0, positions.size() - n_points};
    for (size_t i = 0; i != n_points; ++i) {
      auto random_index = random_idx_generator(random_generator);
      for (uint16_t j = 0; j != n_measurement; ++j) {
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

  // void calculate_errors()
  // {
  //   for (size_t i = 0; i != n_points; ++i) {
  //     for (size_t j : {0, 1, 2}) {
  //       auto square_mean = correlators_squared(j, i);
  //       auto mean_square = std::pow(correlators(j, i), 2);
  //       correlators_errors(j, i) = std::sqrt((square_mean - mean_square) / normalization);
  //     }
  //   }
  // }

  // double log_derivative(auto correlator_row, auto correlator_col, auto lattice_spacing)
  // {
  //   return (correlators(correlator_row, correlator_col) -
  //           correlators(correlator_row, correlator_col + 1)) /
  //          correlators(correlator_row, correlator_col) / lattice_spacing;
  // }

  // double log_derivative_error()
  // {
  // }
};

#endif