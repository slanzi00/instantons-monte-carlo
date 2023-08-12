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
  std::uniform_real_distribution<double> probability(0., 1.);
  int a = 0;
  int r = 0; 
  for (int i = 0; i != n_updates; ++i) {
    for (size_t j = 1; j != positions.size() - 1; ++j) {
      double initial_action = calculate_action();
      double dx = gaussian_step(rng);
      positions[j] += dx;
      double final_action = calculate_action();
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

  // cooling routines /////////////////////////////////////////////////////////////////
  std::ofstream fdata;
  fdata.open("points.txt");
  evolve_using_metropolis<10000>(rng, c, potential, 0.05);
  for (int i=0; i != c.size(); ++i) {
    fdata << i * 0.05 <<"  "<< c[i] << '\n';
  }

  //fdata << '\n';
  //fdata << '\n';
  fdata.close();
  /*
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
  fdata.close();*/
  ///////////////////////////////////////////////////////////////////////////////////
  std::uniform_real_distribution<double> probability(0., 1.);
  int constexpr number_points = 20;
  int constexpr number_measurements = 5;
  std::array<double,number_points> x_cor{}; 
  std::array<double,number_points> x_cor_square{}; 
  std::array<double,number_points> x2_cor{}; 
  std::array<double,number_points> x2_cor_square{}; 
  std::array<double,number_points> x3_cor{}; 
  std::array<double,number_points> x3_cor_square{};
  auto square = [&](double sum_so_far, double a) { return sum_so_far + std::pow(a,2); };
  auto third = [&](double sum_so_far, double a) { return sum_so_far + std::pow(a,3); };
  auto fourth = [&](double sum_so_far, double a) { return sum_so_far + std::pow(a,4); };
  auto sixth = [&](double sum_so_far, double a) { return sum_so_far + std::pow(a,6); };
  // correlation functions
  for (int k = 0; k != number_points; ++k) {
    std::array<double,number_measurements>meas_x{};
    for (int l = 0; l != number_measurements; ++l){
    int p0 = static_cast<int>( (800 - number_points) * probability(rng) ); // WHAT IS HAPPENING HERE?   
      meas_x[l] = ( c[p0] * c[p0 + k] );
    }
    x_cor[k] = (1. / (number_measurements * 1.) ) * ( std::accumulate(meas_x.begin(), meas_x.end(), 0.) );    // mean values
    x_cor_square[k] = (1. / (number_measurements * 1.) ) * ( std::accumulate(meas_x.begin(), meas_x.end(), 0., square) );
    x2_cor[k] = (1. / (number_measurements * 1.) ) * ( std::accumulate(meas_x.begin(), meas_x.end(), 0., square) );
    x2_cor_square[k] = (1. / (number_measurements * 1.)) * ( std::accumulate(meas_x.begin(), meas_x.end(), 0., fourth) );
    x3_cor[k] = (1. / (number_measurements * 1.) ) * ( std::accumulate(meas_x.begin(), meas_x.end(), 0., third) );
    x3_cor_square[k] = (1. / (number_measurements * 1.) ) * ( std::accumulate(meas_x.begin(), meas_x.end(), 0., sixth) );
  }   
  std::ofstream cor_data;
  cor_data.open("correlations.txt");
  for (int i=0; i != x_cor.size(); ++i) {
    cor_data << i * 0.05 <<"  "<< x_cor[i] << '\n';
  }

  cor_data << '\n';
  cor_data << '\n';

  for (int i=0; i != x_cor.size(); ++i) {
    cor_data << i * 0.05 <<"  "<< x2_cor[i] << '\n';
  }

  cor_data << '\n';
  cor_data << '\n';

  for (int i=0; i != x_cor.size(); ++i) {
    cor_data << i * 0.05 <<"  "<< x3_cor[i] << '\n';
  }
  cor_data.close();
}
