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

template <int n_updates, int unconsidered_configurations, int thermalization_steps, class Array, class Potential, class Histo, class Correlator>
void evolve_using_metropolis(double lattice_spacing,
                             std::mt19937& rng,
                             Array& positions,
                             Array& position_cooled,
                             Potential const& potential,
                             int number_bins, // histogram
                             double max_value,
                             Histo& histo_array,
                             int number_points, // correlation functions
                             int number_meas,
                             Correlator& x_cor,
                             Correlator& x2_cor,
                             Correlator& x3_cor,
                             Correlator& x_cor_err,
                             Correlator& x2_cor_err,
                             Correlator& x3_cor_err,
                             Correlator& x_cooled,
                             Correlator& x2_cooled,
                             Correlator& x3_cooled,
                             Correlator& x_cooled_err,
                             Correlator& x2_cooled_err,
                             Correlator& x3_cooled_err
                             )                            
{
  auto calculate_action (Array& path, double lattice_spacing, Potential const& potential) {
    return [&]() {
        double action = 0.;
          for (size_t i = 0; i != path.size(); ++i) {
            double next = (i + path.size() - 1) % (path.size()); //boundary conditions
            double prev = (i + 1) % (path.size());
            double x_pm{(path[i] - path[prev]) / lattice_spacing};
            double x_pp{(path[next] - path[i]) / lattice_spacing};
            double kinetic{(1. / 4.) * (std::pow(x_pm, 2) + std::pow(x_pp, 2))};
            action += lattice_spacing * (kinetic + potential(path[i]));
            }
            return action;
          };
  }//end calculate_action

  auto calculate_correlations ( std::mt19937& rng,
                                int number_meas, 
                                int number_points,
                                Array& path, 
                                Correlator& x, 
                                Correlator& x2, 
                                Correlator& x3, 
                                Correlator& x_square,
                                Correlator& x2_square, 
                                Correlator& x3_square) { return [&]() {
      for (size_t k = 0; k != number_meas; ++k) {  
        ncor += 1;
        int ip0 = floor((800. - number_points) * probability(rng));   //specify site visiting order
        for (size_t l = 0; l != number_points; ++l) {    
          double xcor = path[ip0] * path[ip0 + l];
          x[l] += xcor;
          x_square[l] += std::pow(xcor,2);
          x2[l] += std::pow(xcor,2);
          x2_square[l] += std::pow(xcor,4);
          x3[l] += std::pow(xcor,3);
          x3_square[l] += std::pow(xcor,6);
        }
      }
    };
  }//end calculate_correlations

  auto means_errors_correlations (int ncor,
                                  Correlator& x, 
                                  Correlator& x2,  
                                  Correlator& x3, 
                                  Correlator& x_square,
                                  Correlator& x2_square, 
                                  Correlator& x3_square,
                                  Correlator& x_err,
                                  Correlator& x2_err,
                                  Correlator& x3_err) {
    std::cout <<"ncor ="<< (ncor*1.) << '\n';
    for (size_t i = 0; i != x_cor_err.size(); ++i){ // means and errors of corr functions
      x_cor[i] = x_cor[i]/(ncor*1.);
      double del2 = (x_cor_square[i] / std::pow((ncor*1.),2)) - (std::pow(x_cor[i],2) / (ncor*1.));
      if (del2 < 0.)del2 = 0.;
      x_cor_err[i] = std::sqrt(del2);
    }
    for (size_t i = 0; i != x2_cor_err.size(); ++i){
      x2_cor[i] = x2_cor[i]/(ncor*1.);
      double del2 = (x2_cor_square[i] / std::pow((ncor*1.),2)) - (std::pow(x2_cor[i],2) / (ncor*1.));
      if (del2 < 0.)del2 = 0.;
      x2_cor_err[i] = std::sqrt(del2);
    }
    for (size_t i = 0; i != x3_cor_err.size(); ++i){
      x3_cor[i] = x3_cor[i]/(ncor*1.);
      double del2 = (x3_cor_square[i] / std::pow((ncor*1.),2)) - (std::pow(x3_cor[i],2) / (ncor*1.));
      if (del2 < 0.)del2 = 0.;
      x3_cor_err[i] = std::sqrt(del2);
    }
  }

  auto filling_histo = [&]() {
    double bin_w = (2. * max_value) / (number_bins * 1.);
    for(int i=0; i != positions.size(); ++i) {
        int pos = floor((positions[i] + max_value)/bin_w);
        if(pos < 0)pos = 0;
        if(pos > number_bins)pos = number_bins;
        histo_array[pos] += 1.;
    }
  };//end filling histo

  std::normal_distribution<double> gaussian_step(0., 0.5);
  std::uniform_real_distribution<double> probability(0., 1.);
  int ncor;
  Correlator x_cor_square{};
  Correlator x2_cor_square{};
  Correlator x3_cor_square{};

  Correlator x_square_cooled{};
  Correlator x2_square_cooled{};
  Correlator x3_square_cooled{};
  
  for (int i = 0; i != thermalization_steps + n_updates; ++i) { //Monte Carlo
    for (size_t j = 0; j != positions.size(); ++j) {
      double initial_action = calculate_action(positions, lattice_spacing, potential);
      double initial_action_cooled = calculate_action(position_cooled, lattice_spacing, potential);
      double dx = gaussian_step(rng);
      positions[j] += dx;
      position_cooled[j] += dx;
      double final_action = calculate_action(positions, lattice_spacing, potential);
      double final_action_cooled = calculate_action(position_cooled, lattice_spacing, potential);
      if (std::exp(-(final_action - initial_action)) <= 1. &&
          std::exp(-(final_action - initial_action)) <= (probability(rng))) {
        positions[j] -= dx;
      }
      if (initial_action >= final_action) {
        position_cooled[j] -= dx;
      }
    } 
    if (i > thermalization_steps && std::fmod(i,unconsidered_configurations) == 0) //consider only certain configurations after certain number of steps
    {
      filling_histo();   //filling histogram
      calculate_correlations(rng,number_meas, number_points, positions, 
                            x_cor, x2_cor, x3_cor, x_cor_square, x2_cor_square, x3_cor_square);  // correlation functions
      calculate_correlations(rng,number_meas, number_points, position_cooled, 
                            x_cooled, x2_cooled, x3_cooled, x_square_cooled, x2_square_cooled, x3_square_cooled);  // cooled correlation functions
    } //end if
  } //end Monte Carlo
  means_errors_correlations(ncor, x_cor, x2_cor, x3_cor, x_cor_square, x2_cor_square, x3_cor_square, x_cor_err, x2_cor_err, x3_cor_err);
  means_errors_correlations(ncor, x_cooled, x2_cooled, x3_cooled, x_square_cooled, x2_square_cooled, x3_square_cooled, x_cooled_err, x2_cooled_err, x3_cooled_err);
  
}


auto main() -> int
{
  std::random_device rd;
  std::mt19937 rng(rd());
  double constexpr spacing = 0.05;
  //
  int constexpr total_sweeps = 20000; //total number of monte carlo sweeps
  //
  int constexpr take_every = 2; //configurations to discard between every measurement
  //
  int constexpr unconsider = 200; //first configurations to discard for equilibration purpose, before start
  //
  int constexpr np = 30; //number of points in which cor functions are evaluated
  //
  int constexpr nm = 400; //number measurements per sweep
  //
  auto potential = [](double x) { return std::pow(std::pow(x,2) - std::pow(1.4,2),2); };
  auto c{generate_random_array<800>(rng, -1.4, 1.4)}; // disordered (hot) start
  auto cool{generate_random_array<800>(rng, -1.4, 1.4)}; //cooling array
  std::array<double,np> xc{};
  std::array<double,np> x2{};
  std::array<double,np> x3{};
  std::array<double,np> xc_er{};
  std::array<double,np> x2_er{};
  std::array<double,np> x3_er{};

  std::array<double,np> xc_c{};
  std::array<double,np> x2_c{};
  std::array<double,np> x3_c{};
  std::array<double,np> xc_er_c{};
  std::array<double,np> x2_er_c{};
  std::array<double,np> x3_er_c{};

  int constexpr nbins = 60; //number of bins of histogram
  //
  double constexpr max_val = 2.5; //max x-range value for histogram
  //
  double constexpr bin_w = (2.*max_val) / (nbins*1.); //bin width
  //
  std::array<double,nbins> h{};

  evolve_using_metropolis<total_sweeps, take_every, unconsider>(spacing, rng, c, cool, potential, nbins, max_val, h, 
                            np, nm, xc, x2, x3, xc_er, x2_er, x3_er, xc_c, x2_c, x3_c, xc_er_c, x2_er_c, x3_er_c);
  
  std::ofstream fdata;
  fdata.open("points.txt");
  for (auto const i : c) {
    fdata << i << '\n';
  }
  fdata << '\n';
  fdata << '\n';
  for (auto const i : cool) {
    fdata << i << '\n';
  }
  fdata.close();
}
