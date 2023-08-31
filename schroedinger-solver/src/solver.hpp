#ifndef SOLVER_KERNELS_TBB_HPP
#define SOLVER_KERNELS_TBB_HPP

#include <iostream>

#include <Eigen/Eigenvalues>

#include "initial_conditions.hpp"

struct TridiagonalHamiltonian
{
  Eigen::VectorXd main_diagonal;
  Eigen::VectorXd sub_diagonal;
};

class SchroedingerSolver
{
  Eigen::VectorXd m_positions;
  TridiagonalHamiltonian m_hamiltonian;

 public:
  template <uint16_t polynomial_potential_degree>
  SchroedingerSolver(InitialConditions<polynomial_potential_degree>& ic)
  {
    m_positions = Eigen::VectorXd::Zero(ic.n_points());
    m_hamiltonian.main_diagonal = Eigen::VectorXd::Zero(ic.n_points());
    m_hamiltonian.sub_diagonal = Eigen::VectorXd::Zero(ic.n_points() - 1);

    for (int i{0}; i != ic.n_points(); ++i) {
      m_positions(i) = ic.x_min() + static_cast<double>(i) * ic.get_step();
      m_hamiltonian.main_diagonal(i) = 2. * ic.kinetic() + ic.potential(m_positions(i));
    }
    m_hamiltonian.sub_diagonal.setConstant(-ic.kinetic());
  }

  Eigen::VectorXd get_positions() const
  {
    return m_positions;
  }

  Eigen::VectorXd energy_eigenvalues()
  {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_calculator;
    eigen_calculator.computeFromTridiagonal(
        m_hamiltonian.main_diagonal, m_hamiltonian.sub_diagonal, Eigen::EigenvaluesOnly);
    return eigen_calculator.eigenvalues();
  }

  Eigen::MatrixXd wavefunctions()
  {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_calculator;
    eigen_calculator.computeFromTridiagonal(m_hamiltonian.main_diagonal,
                                            m_hamiltonian.sub_diagonal);
    return eigen_calculator.eigenvectors();
  }

  template <uint8_t position_pow>
  Eigen::VectorXd position_matrix_elements(Eigen::MatrixXd const& wavefunctions)
  {
    size_t num_elements = m_positions.size();
    Eigen::VectorXd n_pole = Eigen::VectorXd::Zero(num_elements);
    for (size_t n{0}; n != num_elements; ++n) {
      for (size_t i{0}; i != num_elements; ++i) {
        n_pole(n) +=
            wavefunctions(i, 0) * std::pow(m_positions(i), position_pow) * wavefunctions(i, n);
      }
    }
    return n_pole.array().pow(2);
  }

  Eigen::VectorXd correlator(Eigen::VectorXd const& multipole, Eigen::VectorXd const& eigenvalues)
  {
    int n_tau = 250;
    Eigen::VectorXd correlator = Eigen::VectorXd::Zero(n_tau);
    for (int tau{0}; tau != n_tau; ++tau) {
      double result = 0.;
      for (int n{0}; n != eigenvalues.size(); ++n) {
        result += multipole(n) * std::exp(-(eigenvalues(n) - eigenvalues(0)) * tau);
      }
      correlator(tau) = result;
    }
    return correlator;
  }
};

#endif