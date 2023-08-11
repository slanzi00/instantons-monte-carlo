#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>

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
  int a = 0;
  int r = 0;
  for (int i = 0; i != n_updates; ++i) {
    for (size_t j = 1; j != positions.size() - 1; ++j) {
      double initial_action = calculate_action();
      double dx = gaussian_step(rng);
      positions[j] += dx;
      double final_action = calculate_action();
      std::uniform_real_distribution<double> probability(0., 1.);
      if (std::exp(-(final_action - initial_action)) <= 1 &&
          std::exp(-(final_action - initial_action)) <= probability(rng)) {
        positions[j] -= dx;
      }
    }
  }
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
  int a = 0;
  int r = 0;
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
  auto h{generate_random_array<800>(rng, -1.4, 1.4)}; // disordered (hot) start
  std::array<double,800> c{1.4}; // ordered (cold) start
  auto potential = [](double x) { return std::pow(std::pow(x, 2) - std::pow(1.4, 2), 2); };
  std::ofstream fdata;
  fdata.open("points.txt");
  evolve_using_metropolis<1000>(rng, c, potential, 0.05);
  for (int i=0; i != c.size(); ++i) {
    fdata << i * 0.05 <<"  "<< c[i] << '\n';
  }

  fdata << '\n';
  fdata << '\n';

  cooling<100000>(rng, c, potential, 0.05);
  for (int i=0; i != c.size(); ++i) {
    fdata << i * 0.05 <<"  "<< c[i] << '\n';
  }
  
  fdata << '\n';
  fdata << '\n';

  evolve_using_metropolis<1000>(rng, h, potential, 0.05);
  for (int i=0; i != h.size(); ++i) {
    fdata << i * 0.05 <<"  "<< h[i] << '\n';
  }

  fdata << '\n';
  fdata << '\n';

  cooling<100000>(rng, h, potential, 0.05);
  for (int i=0; i != h.size(); ++i) {
    fdata << i * 0.05 <<"  "<< h[i] << '\n';
  }
  fdata.close();
}
