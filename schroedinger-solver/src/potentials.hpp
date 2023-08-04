#ifndef POTENTIALS_HPP
#define POTENTIALS_HPP

#include <cmath>

namespace solver {

class Potential
{
 public:
  virtual double operator()(double x) const;
};

class HarmonicPotential : public Potential
{
 public:
  double operator()(double x) const override;
};

class AnharmonicPotential : public Potential
{
  double m_eta;

 public:
  AnharmonicPotential(double eta);
  double operator()(double x) const override;
};

};  // namespace solver

#endif