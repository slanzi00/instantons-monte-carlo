#ifndef POTENTIALS_HPP
#define POTENTIALS_HPP

#include <array>
#include <cmath>

#include <Eigen/Core>

template <uint16_t degree>
struct PolynomialPotential
{
  Eigen::Matrix<double, degree + 1, 1> coefficients;

  double operator()(double x)
  {
    double result = 0.;
    for (int pow = degree; pow >= 0; --pow) {
      result += coefficients(pow) * std::pow(x, pow);
    }
    return result;
  }
};

#endif