#include "potentials.hpp"

namespace solver {

double Potential::operator()(double x) const
{
  x = 0.;
  return x;
}

double HarmonicPotential::operator()(double x) const
{
  return 0.5 * std::pow(x, 2);
}

AnharmonicPotential::AnharmonicPotential(double eta) : m_eta{eta}
{
}

double AnharmonicPotential::operator()(double x) const
{
  return std::pow(std::pow(x, 2) - std::pow(m_eta, 2), 2);
}

};  // namespace solver