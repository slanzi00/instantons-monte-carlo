#include <array>
#include <cmath>
#include <iostream>
#include <random>

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
std::array<float, n_updates> evolve_using_metropolis(std::mt19937& rng,
                                                     Array& positions,
                                                     Potential const& potential,
                                                     float lattice_spacing)
{
  auto calculate_action = [&] {
    float action = 0.;
    for (size_t i = 1; i != positions.size(); ++i) {
      action += 1. / (4. * lattice_spacing) * std::pow(positions[i] - positions[i - 1], 2) +
                lattice_spacing * potential(positions[i]);
    }
    return action;
  };

  std::normal_distribution<float> gaussian_step(0., .5);
  std::array<float, n_updates> x_vals{};

  for (int i = 0; i != n_updates; ++i) {
    float initial_action = calculate_action();
    float dx = gaussian_step(rng);
    size_t random_index = std::uniform_int_distribution<size_t>(0, positions.size() - 1)(rng);
    positions[random_index] += dx;
    float final_action = calculate_action();
    if (std::exp(-(final_action - initial_action)) > 1) {
      positions[random_index] -= dx;
    }
    auto x_val = (1. / 800.) * std::accumulate(positions.begin(), positions.end(), 0.);
    x_vals[i] = x_val;
  }
  return x_vals;
}

auto main() -> int
{
  std::random_device rd;
  std::mt19937 rng(rd());
  auto p{generate_random_array<800>(rng, -1.4, 1.4)};
  auto potential = [](float x) { return std::pow(std::pow(x, 2) - std::pow(1.4, 2), 2); };
  auto a = evolve_using_metropolis<10000>(rng, p, potential, 0.05);
  for (auto const i : a) {
    std::cout << i << '\n';
  }
}