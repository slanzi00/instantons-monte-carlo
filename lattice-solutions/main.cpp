#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <algorithm>

template <size_t array_size>
std::array<float, array_size> generate_random_array(std::mt19937& rng, float x_min, float x_max)
{
  static_assert(array_size > 0, "Array dimension must be > 0");
  std::uniform_real_distribution<float> dist(x_min, x_max);
  std::array<float, array_size> random_array{};
  for (size_t i = 0; i != array_size; ++i) {
    random_array[i] = dist(rng);
  }
  return random_array;
}

template <int n_updates, class Array, class Potential>
void evolve_using_metropolis(std::mt19937& rng,
                             Array& positions,
                             Potential const& potential,
                             float lattice_spacing,
                             int number_points,
                             int number_meas,
                             std::array<float,30>& x_cor,
                             std::array<float,30>& x_cor_square,
                             std::array<float,30>& x2_cor,
                             std::array<float,30>& x2_cor_square,
                             std::array<float,30>& x3_cor,
                             std::array<float,30>& x3_cor_square)                            
{
  auto calculate_action = [&]() {
    float action = 0.f;
    for (size_t i = 0; i != positions.size(); ++i) {
      int next = (i + positions.size() - 1) % (positions.size()); //boundary conditions
      int prev = (i + 1) % (positions.size());
      float x_pm{(positions[i] - positions[prev]) / lattice_spacing};
      float x_pp{(positions[next] - positions[i]) / lattice_spacing};
      float kinetic{(1.f / 4.f) * static_cast<float>(std::pow(x_pm, 2) + std::pow(x_pp, 2))};
      action += lattice_spacing * (kinetic + potential(positions[i]));
    }
    return action;
  };

  std::normal_distribution<float> gaussian_step(0.f, 0.5f);
  
  for (int i = 0; i != n_updates; ++i) {
    //
    //simultaneus random generation
    auto randm{generate_random_array<(800 * 2)>(rng, 0.f, 1.f)};
    //
    //one sweep before actually taking observables
    for (size_t j = 0; j != positions.size(); ++j) {
      float initial_action = calculate_action();
      float dx = gaussian_step(rng);
      positions[j] += dx;
      float final_action = calculate_action();
      if (std::exp(-(final_action - initial_action)) <= 1 &&
          std::exp(-(final_action - initial_action)) <= (randm[j])) {
        positions[j] -= dx;
      }
    }
    //
    //actual configuration considered
    for (size_t j = 0; j != positions.size(); ++j) {
      float initial_action = calculate_action();
      float dx = gaussian_step(rng);
      positions[j] += dx;
      float final_action = calculate_action();
      if (std::exp(-(final_action - initial_action)) <= 1 &&
          std::exp(-(final_action - initial_action)) <= (randm[800 + j])) {
        positions[j] -= dx;
      }
    }
     // correlation functions
    std::uniform_real_distribution<float> probability(0.f, 1.f);
    for (int k = 0; k != number_meas; ++k) {
      //
      //specify site visiting order
      int ip0 = int ((800 - number_points) * probability(rng));
      float x0 = positions[ip0];
      //
      for (int l = 0; l != number_points; ++l) {    
        float x1 = positions[ip0 + l];
        float xcor = x0 * x1;
        x_cor[l] += xcor;
        x_cor_square[l] += static_cast<float>(std::pow(xcor,2));
        x2_cor[l] += static_cast<float>(std::pow(xcor,2));
        x2_cor_square[l] += static_cast<float>(std::pow(xcor,4));
        x3_cor[l] += static_cast<float>(std::pow(xcor,3));
        x3_cor_square[l] += static_cast<float>(std::pow(xcor,6));
      }
    }   
  } //end of MH
  //averages for correlation functions 
  std::transform(x_cor.begin(), x_cor.end(), x_cor.begin(), [&](auto i) { return i / ( (n_updates * number_meas) * 1.f); });
  std::transform(x2_cor.begin(), x2_cor.end(), x2_cor.begin(), [&](auto i) { return i / ( (n_updates * number_meas) * 1.f); });
  std::transform(x3_cor.begin(), x3_cor.end(), x3_cor.begin(), [&](auto i) { return i / ( (n_updates * number_meas) * 1.f); });
}

template <int n_updates, class Array, class Potential>
void cooling(std::mt19937& rng,
                             Array& positions,
                             Potential const& potential,
                             float lattice_spacing)
{
  auto calculate_action = [&] {
    float action = 0.;
    for (size_t i = 1; i != positions.size() - 1; ++i) {
      float x_pm{(positions[i] - positions[i - 1]) / lattice_spacing};
      float x_pp{(positions[i + 1] - positions[i]) / lattice_spacing};
      float kinetic{(1. / 4.) * (std::pow(x_pm, 2) + std::pow(x_pp, 2))};
      action += lattice_spacing * (kinetic + potential(positions[i]));
    }
    return action;
  };

  std::normal_distribution<float> gaussian_step(0., 0.5);
  for (int i = 0; i != n_updates; ++i) {
    for (size_t j = 1; j != positions.size() - 1; ++j) {
      float initial_action = calculate_action();
      float dx = gaussian_step(rng);
      positions[j] += dx;
      float final_action = calculate_action();
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
  std::array<float,800> c{1.4}; // ordered (cold) start
  auto potential = [](float x) { return std::pow(std::pow(x, 2) - std::pow(1.4, 2), 2); };
  //
  //number of points in which cor functions are evaluated
  int constexpr np = 30;
  //
  //number measurements per sweep
  int constexpr nm = 10;
  //      
  std::array<float,np> xc{}; 
  std::array<float,np> xc_sq{}; 
  std::array<float,np> x2{}; 
  std::array<float,np> x2_sq{}; 
  std::array<float,np> x3{}; 
  std::array<float,np> x3_sq{};
  evolve_using_metropolis<10000>(rng, c, potential, 0.05, np, nm, xc, xc_sq, x2, x2_sq, x3, x3_sq);
  std::ofstream fdata;
  fdata.open("points.txt");
  for (auto const i : c) {
    fdata << i << '\n';
  }
  fdata.close();
  std::ofstream cordata;
  cordata.open("correlations.txt");
  for (int i=0; i != np; ++i) {
    cordata << i * 0.05 <<"  "<< xc[i] << '\n';
  }

  cordata << '\n';
  cordata << '\n';

  for (int i=0; i != np; ++i) {
    cordata << i * 0.05 <<"  "<< x2[i] << '\n';
  }

  cordata << '\n';
  cordata << '\n';

  for (int i=0; i != np; ++i) {
    cordata << i * 0.05 <<"  "<< x3[i] << '\n';
  }

}
