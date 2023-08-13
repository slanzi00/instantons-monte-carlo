#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <algorithm>

template <size_t array_size>
std::array<double, array_size> generate_random_array(std::mt19937& rng, double x_min, double x_max)
{
  static_assert(array_size > 0, "Array dimension must be > 0");
  std::uniform_real_distribution<double> dist(x_min, x_max);
  std::array<double, array_size> random_array{};
  for (size_t i = 0; i != array_size; ++i) {
    random_array[i] = dist(rng);
  }
  return random_array;
}

template <int n_updates, class Array, class Potential>
void evolve_using_metropolis(std::mt19937& rng,
                             Array& positions,
                             Potential const& potential,
                             double lattice_spacing,
                             int number_points,
                             int number_meas,
                             std::array<double,20>& x_cor,
                             std::array<double,20>& x_cor_square,
                             std::array<double,20>& x2_cor,
                             std::array<double,20>& x2_cor_square,
                             std::array<double,20>& x3_cor,
                             std::array<double,20>& x3_cor_square)                            
{
  auto calculate_action = [&]() {
    double action = 0.;
    for (size_t i = 0; i != positions.size(); ++i) {
      int next = (i + positions.size() - 1) % (positions.size()); //boundary conditions
      int prev = (i + 1) % (positions.size());
      double x_pm{(positions[i] - positions[prev]) / lattice_spacing};
      double x_pp{(positions[next] - positions[i]) / lattice_spacing};
      double kinetic{(1. / 4.) * (std::pow(x_pm, 2) + std::pow(x_pp, 2))};
      action += lattice_spacing * (kinetic + potential(positions[i]));
    }
    return action;
  };

  std::normal_distribution<double> gaussian_step(0., 0.5);
  
  for (int i = 0; i != n_updates; ++i) {
    //
    //simultaneus random generation
    auto randm{generate_random_array<(800)>(rng, 0., 1.)};
    //
    for (size_t j = 0; j != positions.size(); ++j) {
      double initial_action = calculate_action();
      double dx = gaussian_step(rng);
      positions[j] += dx;
      double final_action = calculate_action();
      if (std::exp(-(final_action - initial_action)) <= 1 &&
          std::exp(-(final_action - initial_action)) <= (randm[j])) {
        positions[j] -= dx;
      }
    }
     // correlation functions
    std::uniform_real_distribution<double> probability(0., 1.);
    for (int k = 0; k != number_meas; ++k) {
      //
      //specify site visiting order
      double ip0 = floor((800 - number_points) * probability(rng));
      double x0 = positions[ip0];
      //
      for (int l = 0; l != number_points; ++l) {    
        double x1 = positions[ip0 + l];
        double xcor = x0 * x1;
        x_cor[k] += xcor;
        x_cor_square[k] += std::pow(xcor,2);
        x2_cor[k] += std::pow(xcor,2);
        x2_cor_square[k] += std::pow(xcor,4);
        x3_cor[k] += std::pow(xcor,3);
        x3_cor_square[k] += std::pow(xcor,6);
      }
    }   
  } //end of MH
  //averages for correlation functions
  std::transform(x_cor.begin(), x_cor.end(), x_cor.begin(), [&](auto i) { return i / ( (n_updates + number_points) * 1.); });
  std::transform(x2_cor.begin(), x2_cor.end(), x2_cor.begin(), [&](auto i) { return i / ( (n_updates + number_points) * 1.); });
  std::transform(x3_cor.begin(), x3_cor.end(), x3_cor.begin(), [&](auto i) { return i / ( (n_updates + number_points) * 1.); });
}

template <int n_updates, class Array, class Potential>
void cooling(std::mt19937& rng,
                             Array& positions,
                             Potential const& potential,
                             double lattice_spacing)
{
  auto calculate_action = [&] {
    double action = 0.;
    for (size_t i = 1; i != positions.size() - 1; ++i) {
      double x_pm{(positions[i] - positions[i - 1]) / lattice_spacing};
      double x_pp{(positions[i + 1] - positions[i]) / lattice_spacing};
      double kinetic{(1. / 4.) * (std::pow(x_pm, 2) + std::pow(x_pp, 2))};
      action += lattice_spacing * (kinetic + potential(positions[i]));
    }
    return action;
  };

  std::normal_distribution<double> gaussian_step(0., 0.5);
  for (int i = 0; i != n_updates; ++i) {
    for (size_t j = 1; j != positions.size() - 1; ++j) {
      double initial_action = calculate_action();
      double dx = gaussian_step(rng);
      positions[j] += dx;
      double final_action = calculate_action();
      if (final_action >= initial_action) {
        positions[j] -= dx;
      }
    }
  }
}

auto main() -> int
{
  std::random_device rd;
  std::mt19937 rng(rd());
  //auto h{generate_random_array<800>(rng, -1.4, 1.4)}; // disordered (hot) start
  std::array<double,800> c{1.4}; // ordered (cold) start
  auto potential = [](double x) { return std::pow(std::pow(x, 2) - std::pow(1.4, 2), 2); };
  //
  //number of points in which cor functions are evaluated
  int constexpr np = 20;
  //
  //number measurements per sweep
  int constexpr nm = 5;
  //      
  std::array<double,np> xc{}; 
  std::array<double,np> xc_sq{}; 
  std::array<double,np> x2{}; 
  std::array<double,np> x2_sq{}; 
  std::array<double,np> x3{}; 
  std::array<double,np> x3_sq{};
  evolve_using_metropolis<10000>(rng, c, potential, 0.05, np, nm, xc, xc_sq, x2, x2_sq, x3, x3_sq);
  std::ofstream fdata;
  fdata.open("points.txt");
  for (auto const i : c) {
    fdata << i << '\n';
  }
  fdata.close();
  
  for (int i=0; i != np; ++i) {
    std::cout << i * 0.05 <<"  "<< xc[i] << '\n';
  }

  std::cout << '\n';
  std::cout << '\n';

  for (int i=0; i != np; ++i) {
    std::cout << i * 0.05 <<"  "<< x2[i] << '\n';
  }

  std::cout << '\n';
  std::cout << '\n';

  for (int i=0; i != np; ++i) {
    std::cout << i * 0.05 <<"  "<< x3[i] << '\n';
  }
}
