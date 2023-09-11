#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>

template <size_t n_points>
struct Correlators
{
  boost::numeric::ublas::matrix<double> correlators;
  boost::numeric::ublas::matrix<double> correlators_squared;

  Correlators()
      : correlators{boost::numeric::ublas::zero_matrix<double>(3, n_points)}
      , correlators_squared{boost::numeric::ublas::zero_matrix<double>(3, n_points)}
  {
  }

  void fill_correlators(auto const& positions, auto& random_generator, uint8_t n_measurement)
  {
    std::uniform_int_distribution<size_t> random_idx_generator{0, positions.size() - n_points - 1};
    for (size_t i{0}; i != n_points; ++i) {
      auto random_index = random_idx_generator(random_generator);
      for (uint8_t j{0}; j != n_measurement; ++j) {
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
};