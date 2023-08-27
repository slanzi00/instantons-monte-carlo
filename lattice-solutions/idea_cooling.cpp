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

template <int n_updates, int unconsidered_configurations, int thermalization_steps,
          int c_updates, int unconsidered_cool, class Array, class Potential, class Density_storer>
void evolve_using_metropolis(std::mt19937& rng,
                             Array& positions,
                             Array& positions_cool,
                             Potential const& potential,
                             double lattice_spacing,
                             Density_storer& n_inst,
                             Density_storer& n_inst_err,
                             Density_storer& action_cooled,
                             Density_storer& action_cooled_err
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

  auto calculate_action_cool = [&]() {
    double action = 0.;
    for (size_t i = 0; i != positions_cool.size(); ++i) {
      double next = (i + positions_cool.size() - 1) % (positions_cool.size()); //boundary conditions
      double prev = (i + 1) % (positions_cool.size());
      double x_pm{(positions_cool[i] - positions_cool[prev]) / lattice_spacing};
      double x_pp{(positions_cool[next] - positions_cool[i]) / lattice_spacing};
      double kinetic{(1. / 4.) * (std::pow(x_pm, 2) + std::pow(x_pp, 2))};
      action += lattice_spacing * (kinetic + potential(positions_cool[i]));
    }
    return action;
  };
  
  std::normal_distribution<double> gaussian_step(0., 0.5);
  std::uniform_real_distribution<double> probability(0., 1.);
  Density_storer n_inst_square{};
  Density_storer action_square{};
  int number_cool_config;

  for (int i = 0; i <= thermalization_steps + n_updates; ++i) { //Monte Carlo
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
    if (i > thermalization_steps && i % (unconsidered_configurations*unconsidered_cool) == 0 )
    { 
      number_cool_config += 1;
      double nin = positions_cool.size()-1;
      for (size_t i=0; i != positions_cool.size(); ++i){
        positions_cool[i] = positions[i];
        double prev = (i + 1) % (positions_cool.size());
        if (std::signbit(positions_cool[i]) == std::signbit(positions_cool[prev])) {
              nin -= 1.;
            }
      }
      double sdiscr = calculate_action_cool();
      n_inst[0] += nin;
      n_inst_square[0] += std::pow(nin,2);
      action_cooled[0] += sdiscr;
      action_square[0] += std::pow(sdiscr,2);
      
      for (int l = 1; l != c_updates; ++l) { //cooling and modify nin
        for (size_t j = 0; j != positions_cool.size(); ++j) {
          double initial_action = calculate_action_cool();
          double dx = gaussian_step(rng);
          positions_cool[j] += dx;
          double final_action = calculate_action_cool();
          if (final_action >= initial_action) {
              positions_cool[j] -= dx;
          }
          int nin = positions_cool.size()-1;
          double prev = (j + 1) % (positions_cool.size());
          if (std::signbit(positions_cool[j]) == std::signbit(positions_cool[prev])) {
              nin -= 1.;
            }
          }
        sdiscr = calculate_action_cool();
        n_inst[l] += nin;
        n_inst_square[l] += std::pow(nin,2);
        action_cooled[l] += sdiscr;
        action_square[l] += std::pow(sdiscr,2);
      }//end cooling procedure
    }//end if
  } //end Monte Carlo
  //instanton density, cooled action
  //std::cout << (n_updates*1.) / ((unconsidered_configurations*unconsidered_cool)*1.) <<"   vs   "<< number_cool_config <<'\n';
  for (size_t i = 0; i != n_inst.size(); ++i){
    n_inst[i] = n_inst[i]/(number_cool_config*1.);
    double del2 = (n_inst_square[i] / std::pow((number_cool_config*1.),2)) - (std::pow(n_inst[i],2) / (number_cool_config*1.));
    if (del2 < 0.)del2 = 0.;
    n_inst_err[i] = std::sqrt(del2);
  }
  for (size_t i = 0; i != action_cooled_err.size(); ++i){
    action_cooled[i] = action_cooled[i]/(number_cool_config*1.);
    double del2 = (action_square[i] / std::pow((number_cool_config*1.),2)) - (std::pow(action_cooled[i],2) / (number_cool_config*1.));
    if (del2 < 0.)del2 = 0.;
    action_cooled_err[i] = std::sqrt(del2);
  }
}

auto main() -> int
{
  std::random_device rd;
  std::mt19937 rng(rd());
  double constexpr spacing = 0.05;
  //
  int constexpr total_sweeps = 20000; //total number of monte carlo sweeps
  //
  int constexpr take_every = 2; //sweep to discard between every configuration used for measurements
  //
  int constexpr equilibration = 200; //first configurations to discard for equilibration purpose, before start
  //
  int constexpr cooling_sweeps = 400;
  //
  int constexpr take_every_cool = 20; //configurations to discard between every cooleing procedure
  //
  auto potential = [](double x) { return std::pow(std::pow(x,2) - std::pow(1.4,2),2); };
  auto c{generate_random_array<800>(rng, -1.4, 1.4)}; // disordered (hot) start
  std::array<double,800>cool{};
  std::array<double,cooling_sweeps> nin{};
  std::array<double,cooling_sweeps> nin_er{};
  std::array<double,cooling_sweeps> action_cool{};
  std::array<double,cooling_sweeps> action_cool_er{};

  evolve_using_metropolis<total_sweeps, take_every, equilibration, cooling_sweeps, take_every_cool>(rng, c, cool, potential, spacing, 
                          nin, nin_er, action_cool, action_cool_er);                          
  
  //instanton density
  double constexpr s0 = (4./3.) * std::pow(1.4,3);
  double constexpr instanton_density_oneloop = 8.* std::sqrt(2.0/M_PI) * std::pow(1.4,2.5) * std::exp(-s0);
  double constexpr instanton_density_twoloops = instanton_density_oneloop * (1. - 71./(72.*s0));
  
  double si_av;
  double del2;
  double si_er;
  for(int i=0; i != nin.size(); ++i) {
    si_av = action_cool[i]/nin[i];
    del2 = std::pow(action_cool_er[i]/action_cool[i],2) + std::pow(nin_er[i]/nin[i],2);
    si_er = si_av * std::sqrt(del2);
    std::cout << i <<" "<< si_av <<" "<< si_er <<" "<< s0 <<'\n';
    }

}
