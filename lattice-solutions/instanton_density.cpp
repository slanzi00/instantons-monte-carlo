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
  for (auto& value : random_array) {
        value = dist(rng);
    }
  return random_array;
}

template <typename Array>
constexpr double calculate_kinetic(const Array& positions, size_t i, const size_t positions_size, const double& current_position, double lattice_spacing) {
    size_t prev = (i + positions_size - 1) % (positions_size);
    size_t next = (i + 1) % (positions_size);
    double x_pm = (current_position - positions[prev]) / lattice_spacing;
    double x_pp = (positions[next] - current_position) / lattice_spacing;
    return 0.25 * (std::pow(x_pm, 2) + std::pow(x_pp, 2));
}

template <typename Array, typename Potential>
double calculate_action(const Array& positions, Potential& potential, const size_t& positions_size, double lattice_spacing) {
    double action = 0.;
    for (size_t i = 0; i < positions_size; ++i) {
        const double& current_position = positions[i];
        double kinetic = calculate_kinetic(positions, i, positions_size, current_position, lattice_spacing);
        action += lattice_spacing * (kinetic + potential(current_position));
    }
    return action;
}

template <typename Array>
void evolve_positions(std::mt19937& rng, Array& positions, double lattice_spacing,
                      double probability_threshold) {
    std::normal_distribution<double> gaussian_step(0.0, 0.5);
    for (double& position : positions) {
        double initial_action = calculate_action(positions, lattice_spacing);
        double dx = gaussian_step(rng);
        position += dx;
        double final_action = calculate_action(positions, lattice_spacing);
        if (std::exp(final_action - initial_action) <= probability_threshold) {
            position -= dx;
        }
    }
}

template <typename T>
void calculate_statistics(T& dataset, T& square_dataset, T& errors_dataset, const double& number_correlator) {
    for (size_t i = 0; i != dataset.size(); ++i){     
    double& dataset_current = dataset[i];
    double& square_current = square_dataset[i];
    dataset_current /= number_correlator;
    square_current = (square_current / (number_correlator * number_correlator)) - ((dataset_current * dataset_current) / number_correlator);
    errors_dataset[i] = std::sqrt(std::max(square_current, 0.));
  }
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
  std::normal_distribution<double> gaussian_step(0., 0.5);
  std::uniform_real_distribution<double> probability(0., 1.);
  Density_storer n_inst_square{};
  Density_storer action_square{};
  const double configurations_cooled = (1.*n_updates)/(1.*(unconsidered_configurations* unconsidered_cool));
  const size_t positions_size = positions_cool.size();

  for (int i = 0; i <= thermalization_steps + n_updates; ++i) { 
    //Monte Carlo
    for (auto& pos : positions) {
        double initial_action = calculate_action(positions, potential, positions_size, lattice_spacing);
        double dx = gaussian_step(rng);
        pos += dx;
        double final_action = calculate_action(positions, potential, positions_size, lattice_spacing);
        double acceptance_prob = std::exp(initial_action - final_action);
        if (acceptance_prob <= 1.0 && acceptance_prob <= probability(rng)) {
            pos -= dx;
        }
    } 
    if (i > thermalization_steps && i % (unconsidered_configurations*unconsidered_cool) == 0 )
    { 
      double nin = 0;
      for (size_t i=0; i != positions_size; ++i){
        positions_cool[i] = positions[i];
        size_t prev = (i + positions_size - 1) % (positions_size);
        size_t next = (i + 1) % (positions_size);
        if (std::signbit(positions_cool[next]) != std::signbit(positions_cool[prev])) {
              nin += 1.;
            }
      }
      double sdiscr = calculate_action(positions, potential, positions_size, lattice_spacing);
      n_inst[0] += nin;
      n_inst_square[0] += std::pow(nin,2);
      action_cooled[0] += sdiscr;
      action_square[0] += std::pow(sdiscr,2);
      
      for (int l = 1; l != c_updates; ++l) { //cooling and modify nin
        double number_in = 0;
        for (size_t j = 0; j != positions_size; ++j) {
          double initial_action = calculate_action(positions_cool, potential, positions_size, lattice_spacing);
          double dx = gaussian_step(rng);
          double& position = positions_cool[j];
          position += dx;
          double final_action = calculate_action(positions_cool, potential, positions_size, lattice_spacing);
          if (final_action >= initial_action) {
              position -= dx;
          }
          size_t prev = (j + positions_size - 1) % (positions_size);
          size_t next = (j + 1) % (positions_size);
          if (std::signbit(positions_cool[next]) != std::signbit(positions_cool[prev])) {
              number_in += 1.;
            }
          }
        sdiscr = calculate_action(positions_cool, potential, positions_size, lattice_spacing);
        n_inst[l] += number_in;
        n_inst_square[l] += std::pow(number_in,2);
        action_cooled[l] += sdiscr;
        action_square[l] += std::pow(sdiscr,2);
      }//end cooling procedure
    }//end if
  } //end Monte Carlo
  //instanton density, cooled action

  calculate_statistics(n_inst, n_inst_square, n_inst_err, configurations_cooled);
  calculate_statistics(action_cooled, action_square, action_cooled_err, configurations_cooled);
}

auto main() -> int
{
  std::random_device rd;
  std::mt19937 rng(rd());
  double constexpr eta = 1.4;
  //
  double constexpr spacing = 0.05;
  //
  int constexpr total_sweeps = 8e3; //total number of monte carlo sweeps 
  //                          2e5
  int constexpr take_every = 2; //sweep to discard between every configuration used for measurements
  //
  int constexpr equilibration = 200; //first configurations to discard for equilibration purpose, before start
  //
  int constexpr cooling_sweeps = 660;
  //                            660
  int constexpr take_every_cool = 15; //configurations to discard between every cooleing procedure
  //                             250   
  auto potential = [](const double& x) { return std::pow(std::pow(x,2) - std::pow(eta,2),2); };
  auto c{generate_random_array<800>(rng, -eta, eta)}; // disordered (hot) start
  std::array<double,800>cool{};
  std::array<double,cooling_sweeps> nin{};
  std::array<double,cooling_sweeps> nin_er{};
  std::array<double,cooling_sweeps> action_cool{};
  std::array<double,cooling_sweeps> action_cool_er{};

  evolve_using_metropolis<total_sweeps, take_every, equilibration, cooling_sweeps, take_every_cool>(rng, c, cool, potential, spacing, 
                          nin, nin_er, action_cool, action_cool_er);                          
  
  //instanton density
  double constexpr s0 = (4./3.) * std::pow(eta,3);
  double constexpr instanton_density_oneloop = 8.* std::sqrt(2.0/M_PI) * std::pow(eta,2.5) * std::exp(-s0);
  double constexpr instanton_density_twoloops = instanton_density_oneloop * (1. - 71./(72.*s0));
  
  for(size_t i=0; i != nin.size(); ++i) {
    std::cout << i <<" "<< nin[i] / (spacing * 800.) <<" "<< nin_er[i] / (spacing * 800.) <<" "<< 
                instanton_density_oneloop  <<" "<< 
                instanton_density_twoloops  <<'\n';
  }
  std::cout << '\n' << '\n';
  
  for(size_t i=0; i != nin.size(); ++i) {
    
    std::cout << i <<" "<< action_cool[i] <<" "<<
             action_cool_er[i] <<" "<< nin[i] * s0 <<'\n';
  }
  std::cout << '\n' << '\n';

  double si_av;
  double del2;
  double si_er;
  for(size_t i=0; i != nin.size(); ++i) {
    si_av = action_cool[i] / nin[i];
    del2 = std::pow(action_cool_er[i]/action_cool[i],2) + std::pow(nin_er[i]/nin[i],2);
    si_er = si_av * std::sqrt(del2);
    std::cout << i <<" "<< si_av <<" "<< si_er <<" "<< s0 <<'\n';
    }

}
