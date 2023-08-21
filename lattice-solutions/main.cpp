#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <algorithm>
#include <assert.h>

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

template <int n_updates, class Array, class Potential, class Histo, class Correlator>
void evolve_using_metropolis(std::mt19937& rng,
                             Array& positions,
                             Potential const& potential,
                             double lattice_spacing,
                             int number_points,
                             int number_meas,
                             int number_bins,
                             double max_value,
                             Histo& histo_array,
                             Correlator& x_cor,
                             Correlator& x2_cor,
                             Correlator& x3_cor,
                             Correlator& x_cor_err,
                             Correlator& x2_cor_err,
                             Correlator& x3_cor_err
                             )                            
{
  auto calculate_action = [&]() {
    double action = 0.;
    for (size_t i = 0; i != positions.size(); ++i) {
      double next = (i + positions.size() - 1) % (positions.size()); //boundary conditions
      double prev = (i + 1) % (positions.size());
      double x_pm{(positions[i] - positions[prev]) / lattice_spacing};
      double x_pp{(positions[next] - positions[i]) / lattice_spacing};
      double kinetic{(1. / 4.) * (std::pow(x_pm, 2) + std::pow(x_pp, 2))};
      action += lattice_spacing * (kinetic + potential(positions[i]));
    }
    return action;
  };

  auto filling_histo = [&]() {
    double bin_w = (2. * max_value) / (number_bins * 1.);
    for(int i=0; i != positions.size(); ++i) {
        int pos = floor((positions[i] + max_value)/bin_w);
        if(pos < 0)pos = 0;
        if(pos > number_bins)pos = number_bins;
        histo_array[pos] += 1.;
    }
  };

  std::normal_distribution<double> gaussian_step(0., 0.5);
  std::uniform_real_distribution<double> probability(0., 1.);

  //toss 200 sweeps before actually taking observables  
    for (int h = 0; h != 200; ++h) {
      for (size_t j = 0; j != positions.size(); ++j) {
        double initial_action = calculate_action();
        double dx = gaussian_step(rng);
        positions[j] += dx;
        double final_action = calculate_action();
          if (std::exp(-(final_action - initial_action)) <= 1. &&
            std::exp(-(final_action - initial_action)) <= (probability(rng))) {
            positions[j] -= dx;
          }
        }
      }  

  Correlator x_cor_square{};
  Correlator x2_cor_square{};
  Correlator x3_cor_square{};

  for (int i = 0; i != n_updates; ++i) {
    //one sweep before taking observables
    for (size_t j = 0; j != positions.size(); ++j) {
        double initial_action = calculate_action();
        double dx = gaussian_step(rng);
        positions[j] += dx;
        double final_action = calculate_action();
          if (std::exp(-(final_action - initial_action)) <= 1. &&
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
      if (std::exp(-(final_action - initial_action)) <= 1. &&
          std::exp(-(final_action - initial_action)) <= (probability(rng))) {
        positions[j] -= dx;
      }
    } 
    //filling histogram
    filling_histo();   

    // correlation functions
    for (size_t k = 0; k != number_meas; ++k) {
     //
     //specify site visiting order
     int ip0 = floor((800. - number_points) * probability(rng));
     double x0 = positions[ip0];
     //
     for (size_t l = 0; l != number_points; ++l) {    
        double xcor = x0*positions[ip0 + l];
        double x2cor = std::pow(x0*positions[ip0 + l],2);
        x_cor[l] += xcor;
        x_cor_square[l] += std::pow(xcor,2);
        x2_cor[l] += x2cor;
        x2_cor_square[l] += std::pow(x2cor,2);
        x3_cor[l] += std::pow(xcor,3);
        x3_cor_square[l] += std::pow(xcor,6);
      }
    }
       
  } //end of MH

  double ncor = number_meas * n_updates;
  // errors of corr functions
  for (size_t i = 0; i != x_cor_err.size(); ++i){
    x_cor[i] = x_cor[i]/ncor;
    double del2 = (x_cor_square[i] / std::pow(ncor,2)) - (std::pow(x_cor[i],2) / ncor);
    if (del2 < 0.)del2 = 0.;
    x_cor_err[i] = std::sqrt(del2);
  }
  for (size_t i = 0; i != x2_cor_err.size(); ++i){
    x2_cor[i] = x2_cor[i]/ncor;
    double del2 = (x2_cor_square[i] / std::pow(ncor,2)) - (std::pow(x2_cor[i],2) / ncor);
    if (del2 < 0.)del2 = 0.;
    x2_cor_err[i] = std::sqrt(del2);
  }
  for (size_t i = 0; i != x3_cor_err.size(); ++i){
    x3_cor[i] = x3_cor[i]/ncor;
    double del2 = (x3_cor_square[i] / std::pow(ncor,2)) - (std::pow(x3_cor[i],2) / ncor);
    if (del2 < 0.)del2 = 0.;
    x3_cor_err[i] = std::sqrt(del2);
  }
}

auto main() -> int
{
  std::random_device rd;
  std::mt19937 rng(rd());
  double constexpr spacing = 0.05;
  auto c{generate_random_array<800>(rng, -1.4, 1.4)}; // disordered (hot) start
  //std::array<double,800> c{1.4f}; // ordered (cold) start
  auto potential = [](double x) { return std::pow(std::pow(x,2) - std::pow(1.4,2),2); };
  //
  int constexpr monte_carlo_sweeps = 10000;
  //number of points in which cor functions are evaluated
  int constexpr np = 22;
  //
  //number measurements per sweep
  int constexpr nm = 10;
  //
  std::array<double,np> xc{};
  std::array<double,np> x2{};
  std::array<double,np> x3{};
  std::array<double,np> xc_er{};
  std::array<double,np> x2_er{};
  std::array<double,np> x3_er{};

  //histogram
  //
  //number of bins
  int constexpr nbins = 60;
  //
  //max value for histogram
  double constexpr max_val = 2.5;
  //
  //bin width
  double constexpr bin_w = (2.*max_val) / (nbins*1.);
  std::array<double,nbins> h{};

  evolve_using_metropolis<monte_carlo_sweeps>(rng, c, potential, spacing, np, nm, nbins, max_val, h,
                                xc, x2, x3, xc_er, x2_er, x3_er);
  
  std::ofstream fdata;
  fdata.open("points.txt");
  for (auto const i : c) {
    fdata << i << '\n';
  }
  fdata.close();

  std::ofstream cordata;
  cordata.open("correlations.txt");
  
  //log derivative xc
  for (size_t i = 0; i != np-1; ++i) {
    double dx = (xc[i] - xc[i+1]) / xc[i] / spacing;
    double dxe2 = std::pow(xc_er[i+1]/xc[i],2) + std::pow((xc_er[i]*xc[i+1]) / std::pow(xc[i],2),2);
    double dxe = std::sqrt(dxe2) / spacing;
    cordata << i * spacing <<"  "<< xc[i] <<"  "<< xc_er[i] <<"  "<< dx <<"  "<< dxe << '\n';
  }

  cordata << '\n';
  cordata << '\n';

  //log derivative x2
  for (size_t i = 0; i != np-1; ++i) {
    double x2_sub = x2[i] - x2[np-1];
    double x2_sub_next = x2[i+1] - x2[np-1];
    //std::cout << x2[i] <<" - "<< x2[np-1] <<"  =  "<< x2_sub <<'\n';
    double error_sub_next = std::sqrt(std::pow(x2_er[i+1],2) + std::pow(x2_er[np-1],2));
    double error_sub = std::sqrt(std::pow(x2_er[i],2) + std::pow(x2_er[np-1],2));

    double dx = (x2_sub - x2_sub_next) / x2_sub / spacing;
    double dxe2 = std::pow(error_sub_next / x2_sub,2) + std::pow((error_sub*x2_sub_next) / std::pow(x2_sub,2),2);
    double dxe = std::sqrt(dxe2)/spacing;

    cordata << i * spacing <<"  "<< x2[i] <<"  "<< x2_er[i] <<"  "<< dx <<"  "<<  dxe << '\n';
  }

  cordata << '\n';
  cordata << '\n';
  
  //log derivative x3
  for (size_t i = 0; i != np-1; ++i) {
    double dx = (x3[i] - x3[i+1]) / x3[i] / spacing;
    double dxe2 = std::pow(x3_er[i+1] / x3[i],2) + std::pow((x3_er[i]*x3[i+1]) / std::pow(x3[i],2),2);
    double dxe = std::sqrt(dxe2) / spacing;
    cordata << i * spacing <<"  "<< x3[i] <<"  "<< x3_er[i] <<"  "<< dx <<"  "<< dxe << '\n';
  }

  cordata.close();
  std::ofstream probability;
  probability.open("ground_state_probability.txt");

  double xnorm;
  for(int i=0; i != nbins; ++i) {
    xnorm += (h[i] * bin_w);
  }
  for(int i=0; i != nbins; ++i) {
    double xx = (-max_val) + (double)i * bin_w;
    probability << xx <<"  "<< h[i]/xnorm <<'\n';
    }
  probability.close();
}
