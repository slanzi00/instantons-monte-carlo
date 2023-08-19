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
                             std::array<double,30>& x_cor,
                             std::array<double,30>& x_cor_square,
                             std::array<double,30>& x2_cor,
                             std::array<double,30>& x2_cor_square,
                             std::array<double,30>& x3_cor,
                             std::array<double,30>& x3_cor_square,
                             std::array<double,50> hist 
                             )                            
{
  auto calculate_action = [&]() {
    double action = 0.f;
    for (size_t i = 0; i != positions.size(); ++i) {
      int next = (i + positions.size() - 1) % (positions.size()); //boundary conditions
      int prev = (i + 1) % (positions.size());
      double x_pm{(positions[i] - positions[prev]) / lattice_spacing};
      double x_pp{(positions[next] - positions[i]) / lattice_spacing};
      double kinetic{(1.f / 4.f) * (std::pow(x_pm, 2) + std::pow(x_pp, 2))};
      action += lattice_spacing * (kinetic + potential(positions[i]));
    }
    return action;
  };

  std::normal_distribution<double> gaussian_step(0., 0.5);
  std::uniform_real_distribution<double> probability(0., 1.);

  //toss 100 sweeps before actually taking observables
    for (int h = 0; h != 100; ++h) {
      for (size_t j = 0; j != positions.size(); ++j) {
        double initial_action = calculate_action();
        double dx = gaussian_step(rng);
        positions[j] += dx;
        double final_action = calculate_action();
          if (std::exp(-(final_action - initial_action)) <= 1 &&
            std::exp(-(final_action - initial_action)) <= (probability(rng))) {
            positions[j] -= dx;
          }
        }
      }  
  // histogram parameters
  
  //double constexpr minimum_value = -1.5 * 1.4;
  //double constexpr bin_width = 3 * 1.4 / hist.size();

  for (int i = 0; i != n_updates; ++i) {
    //one sweep before taking observables
    for (size_t j = 0; j != positions.size(); ++j) {
        double initial_action = calculate_action();
        double dx = gaussian_step(rng);
        positions[j] += dx;
        double final_action = calculate_action();
          if (std::exp(-(final_action - initial_action)) <= 1 &&
            std::exp(-(final_action - initial_action)) <= (probability(rng))) {
            positions[j] -= dx;
          }
        }
    //actual configuration considered for observables
    //
    for (size_t j = 0; j != positions.size(); ++j) {
      double initial_action = calculate_action();
      double dx = gaussian_step(rng);
      positions[j] += dx;
      double final_action = calculate_action();
      if (std::exp(-(final_action - initial_action)) <= 1 &&
          std::exp(-(final_action - initial_action)) <= (probability(rng))) {
        positions[j] -= dx;
      }
    }
    // filling histogram array
    //
    /*for (size_t j = 0; j != positions.size(); ++j) {
      double candidate_histo = (positions[j] - minimum_value) / bin_width + 1.;
      if (candidate_histo < 0.)candidate_histo = 0.;
      if (candidate_histo > hist.size())candidate_histo = hist.size();
      hist[candidate_histo] += 1;
    }*/
    
    // correlation functions
    for (int k = 0; k != number_meas; ++k) {
     //
     //specify site visiting order
     int ip0 = int ((800 - number_points) * probability(rng));
     double x0 = positions[ip0];
     //
     for (int l = 0; l != number_points; ++l) {    
        double x1 = positions[ip0 + l];
        double xcor = x0 * x1;
        x_cor[l] += xcor;
        x_cor_square[l] += std::pow(xcor,2);
        x2_cor[l] += std::pow(xcor,2);
        x2_cor_square[l] += std::pow(xcor,4);
        x3_cor[l] += std::pow(xcor,3);
        x3_cor_square[l] += std::pow(xcor,6);
      }
    }   
  } //end of MH
}

template <int ncorr, class Array>
void mean_and_errors (                    
                      Array& x_cor,
                      Array& x_cor_square,
                      Array& x_cor_err
                      ) 
{
  for (size_t i = 0; i != x_cor_err.size(); ++i){
    x_cor[i] = x_cor[i] / ncorr;
    double del2 = (x_cor_square[i] / std::pow(ncorr,2)) - (std::pow(x_cor[i],2) / ncorr);
    if (del2 < 0.f)del2 = 0.f;
    x_cor_err[i] = std::sqrt(del2);
  }
}   


template <int n_updates, class Array, class Potential>
void cooling(std::mt19937& rng,
                             Array& positions,
                             Potential const& potential,
                             double lattice_spacing)
{
  auto calculate_action = [&]() {
    double action = 0.f;
    for (size_t i = 0; i != positions.size(); ++i) {
      int next = (i + positions.size() - 1) % (positions.size()); //boundary conditions
      int prev = (i + 1) % (positions.size());
      double x_pm{(positions[i] - positions[prev]) / lattice_spacing};
      double x_pp{(positions[next] - positions[i]) / lattice_spacing};
      double kinetic{(1.f / 4.f) * (std::pow(x_pm, 2) + std::pow(x_pp, 2))};
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
  //auto c{generate_random_array<800>(rng, -1.4, 1.4)}; // disordered (hot) start
  std::array<double,800> c{1.4}; // ordered (cold) start
  auto potential = [](double x) { return std::pow(std::pow(x, 2) - std::pow(1.4, 2), 2); };
  //
  //number of points in which cor functions are evaluated
  int constexpr np = 30;
  //
  //number measurements per sweep
  int constexpr nm = 6;
  //      
  //histogram array
  int constexpr number_bins = 50;
  std::array<double,number_bins> hist{};
  //
  std::array<double,np> xc{};
  std::array<double,np> x2{};
  std::array<double,np> x3{};
  std::array<double,np> xc_sq{};
  std::array<double,np> x2_sq{};
  std::array<double,np> x3_sq{};
  std::array<double,np> xc_er{};
  std::array<double,np> x2_er{};
  std::array<double,np> x3_er{};

  evolve_using_metropolis<10000>(rng, c, potential, 0.05, np, nm, 
                                  xc, xc_sq, x2, x2_sq, x3, x3_sq, 
                                  hist);
  mean_and_errors<10000*nm>(xc, xc_sq, xc_er);
  mean_and_errors<10000*nm>(x2, x2_sq, x2_er);
  mean_and_errors<10000*nm>(x3, x3_sq, x3_er);
  std::ofstream fdata;
  fdata.open("points.txt");
  for (auto const i : c) {
    fdata << i << '\n';
  }
  fdata.close();
  std::ofstream cordata;
  cordata.open("correlations.txt");
  for (int i=0; i != np; ++i) {
    cordata << i * 0.05 <<"  "<< xc[i] <<"  "<< xc_er[i] << '\n';
  }

  cordata << '\n';
  cordata << '\n';

  for (int i=0; i != np; ++i) {
    cordata << i * 0.05 <<"  "<< x2[i] <<"  "<< x2_er[i] << '\n';
  }

  cordata << '\n';
  cordata << '\n';

  for (int i=0; i != np; ++i) {
    cordata << i * 0.05 <<"  "<< x3[i] <<"  "<< x3_er[i] << '\n';
  }
  cordata.close();
  std::ofstream histogram;
  histogram.open("probability.txt");
  for (auto const i : hist) {
    histogram << i << '\n';
  }
  histogram.close();
}
