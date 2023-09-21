#include "correlators.hpp"

// correlators normalization
double normalization = (sv::n_sweeps * static_cast<double>(sv::n_correlator_meas)) /
                       static_cast<double>(sv::take_corr_meas_every);

// initialize data member
Correlators::Correlators()
    : correlators{zeros(3, sv::n_correlator_points)}
    , correlators_squared{zeros(3, sv::n_correlator_points)}
    , correlators_errors{zeros(3, sv::n_correlator_points)}
{
}

void Correlators::fill_correlators(std::array<double, sv::n_lattice_points> const& positions,
                                   std::mt19937 random_generator)
{
  std::uniform_int_distribution random_idx_generator{
      0, sv::n_lattice_points - sv::n_correlator_points - 1};
  for (uint16_t j = 0; j != sv::n_correlator_meas; ++j) {
    // starting form a random index we make measurements to fill the correlator
    size_t random_index = random_idx_generator(random_generator);
    for (size_t i = 0; i != sv::n_correlator_points; ++i) {
      double correlator_i = positions[random_index] * positions[random_index + i];
      correlators(0, i) += correlator_i;
      correlators_squared(0, i) += std::pow(correlator_i, 2);
      correlators(1, i) += std::pow(correlator_i, 2);
      correlators_squared(1, i) += std::pow(correlator_i, 4);
      correlators(2, i) += std::pow(correlator_i, 3);
      correlators_squared(2, i) += std::pow(correlator_i, 6);
    }
  }
}

void Correlators::normalize_correlators()
{
  for (size_t i = 0; i != sv::n_correlator_points; ++i) {
    for (size_t j : {0, 1, 2}) {
      correlators(j, i) /= normalization;
      correlators_squared(j, i) /= normalization;
    }
  }
}

void Correlators::calculate_errors()
{
  for (size_t i = 0; i != sv::n_correlator_points; ++i) {
    for (size_t j : {0, 1, 2}) {
      double square_mean = correlators_squared(j, i);
      double mean_square = std::pow(correlators(j, i), 2);
      correlators_errors(j, i) = std::sqrt((square_mean - mean_square) / normalization);
    }
  }
}

double Correlators::log_derivative(size_t correlator_row, size_t correlator_col)
{
  double c = correlators(correlator_row, correlator_col);
  double c_next = correlators(correlator_row, correlator_col + 1);
  double c_2_last = correlators(1, sv::n_correlator_points - 1);
  if (correlator_row != 1) {
    return (c - c_next) / (c * sv::lattice_spacing);
  } else {
    return (c - c_next) / ((c - c_2_last) * sv::lattice_spacing);
  }
}

double Correlators::log_derivative_error(size_t correlator_row, size_t correlator_col)
{
  double c = correlators(correlator_row, correlator_col);
  double c_next = correlators(correlator_row, correlator_col + 1);
  double c_2_sub =
      correlators(correlator_row, correlator_col) - correlators(1, sv::n_correlator_points - 1);
  double c_2_sub_next =
      correlators(correlator_row, correlator_col + 1) - correlators(1, sv::n_correlator_points - 1);

  double c_err = correlators_errors(correlator_row, correlator_col);
  double c_err_next = correlators_errors(correlator_row, correlator_col + 1);
  double c_2_err_last = correlators_errors(1, sv::n_correlator_points - 1);
  double c_2_err_sub = std::sqrt(std::pow(c_err, 2) + std::pow(c_2_err_last, 2));
  double c_2_err_next_sub = std::sqrt(std::pow(c_err_next, 2) + std::pow(c_2_err_last, 2));

  if (correlator_row != 1) {
    double error_square = std::pow(c_err_next / c, 2) + std::pow((c_err * c_next) / (c * c), 2);
    return std::sqrt(error_square) / sv::lattice_spacing;
  } else {
    double error_square = std::pow(c_2_err_next_sub / c_2_sub, 2) +
                          std::pow((c_2_err_sub * c_2_sub_next) / (c_2_sub * c_2_sub), 2);
    return std::sqrt(error_square) / sv::lattice_spacing;
  }
}