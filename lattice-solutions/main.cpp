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

template <int n_updates, int unconsidered_configurations, int thermalization_steps, int c_updates, int unconsidered_cool, class Array, class Potential, class Histo, class Correlator>
void evolve_using_metropolis(std::mt19937& rng,
                             Array& positions,
                             Array& positions_cool,
                             Potential const& potential,
                             double lattice_spacing,
                             int const number_points,
                             int const number_meas,
                             int const number_bins,
                             double max_value,
                             Histo& histo_array,
                             Histo& histo_array_cool,
                             Correlator& x_cor,
                             Correlator& x2_cor,
                             Correlator& x3_cor,
                             Correlator& x_cor_err,
                             Correlator& x2_cor_err,
                             Correlator& x3_cor_err,
                             Correlator& x_cool,
                             Correlator& x2_cool,
                             Correlator& x3_cool,
                             Correlator& x_err_cool,
                             Correlator& x2_err_cool,
                             Correlator& x3_err_cool
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

  auto filling_histo = [&]() {
    double bin_w = (2. * max_value) / (number_bins * 1.);
    for(int i=0; i != positions.size(); ++i) {
        int pos = floor((positions[i] + max_value)/bin_w);
        if(pos < 0)pos = 0;
        if(pos > number_bins-1)pos = number_bins-1;
        histo_array[pos] += 1.;
    }
  };

  auto filling_histo_cool = [&]() {
    double bin_w = (2. * max_value) / (number_bins * 1.);
    for(int i=0; i != positions_cool.size(); ++i) {
        int pos = floor((positions_cool[i] + max_value)/bin_w);
        if(pos < 0)pos = 0;
        if(pos > number_bins-1)pos = number_bins-1;
        histo_array_cool[pos] += 1.;
    }
  };

  std::normal_distribution<double> gaussian_step(0., 0.5);
  std::uniform_real_distribution<double> probability(0., 1.);
  int const ncor = (n_updates / unconsidered_configurations) * number_meas;
  int const ncor_cool = (n_updates / (unconsidered_configurations * unconsidered_cool)) * number_meas;
  Correlator x_cor_square{};
  Correlator x2_cor_square{};
  Correlator x3_cor_square{};
  Correlator x_square_cool{};
  Correlator x2_square_cool{};
  Correlator x3_square_cool{};
  
  for (int i = 0; i != thermalization_steps + n_updates + 1; ++i) { //Monte Carlo
    for (size_t j = 0; j != positions.size(); ++j) {
      double initial_action = calculate_action();
      double dx = gaussian_step(rng);
      positions[j] += dx;
      double final_action = calculate_action();
      if (std::exp(-(final_action - initial_action)) <= 1. &&
          std::exp(-(final_action - initial_action)) <= (probability(rng))) {
        positions[j] -= dx;
      }
      filling_histo(); 
    }
    if (i > thermalization_steps && (i % unconsidered_configurations) == 0) //consider only certain configurations after certain number of steps
    {   
      for (size_t k = 0; k != number_meas; ++k) {   // correlation functions
        int ip0 = floor((800. - number_points) * probability(rng));   //specify site visiting order
        for (size_t l = 0; l != number_points; ++l) {    
          double xcor = positions[ip0] * positions[ip0 + l];
          x_cor[l] += xcor;
          x_cor_square[l] += std::pow(xcor,2);
          x2_cor[l] += std::pow(xcor,2);
          x2_cor_square[l] += std::pow(xcor,4);
          x3_cor[l] += std::pow(xcor,3);
          x3_cor_square[l] += std::pow(xcor,6);
        }
      }
    }//end if
    if (i > thermalization_steps && (i % (unconsidered_configurations*unconsidered_cool)) == 0)
    { 
      std::copy(positions.begin(), positions.end(), positions_cool.begin()); 
      for (int l = 0; l != c_updates; ++l) { //cooling
        for (size_t j = 0; j != positions_cool.size(); ++j) {
          double initial_action = calculate_action_cool();
          double dx = gaussian_step(rng);
          positions_cool[j] += dx;
          double final_action = calculate_action_cool();
          if (final_action >= initial_action) {
          positions_cool[j] -= dx;
          }
        } 
      }
      for (size_t k = 0; k != number_meas; ++k) {   // correlation functions
        int ip0 = floor((800. - number_points) * probability(rng));   //specify site visiting order
        for (size_t l = 0; l != number_points; ++l) {    
          double xcor = positions_cool[ip0] * positions_cool[ip0 + l];
          x_cool[l] += xcor;
          x_square_cool[l] += std::pow(xcor,2);
          x2_cool[l] += std::pow(xcor,2);
          x2_square_cool[l] += std::pow(xcor,4);
          x3_cool[l] += std::pow(xcor,3);
          x3_square_cool[l] += std::pow(xcor,6);
        }
      }
    }//end if
  } //end Monte Carlo
 
  std::cout <<"ncor ="<< (ncor*1.) << '\n';
  std::cout <<"ncor_cool ="<< (ncor_cool*1.) << '\n';
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

  for (size_t i = 0; i != x_err_cool.size(); ++i){ // means and errors of cooled corr functions
    x_cool[i] = x_cool[i]/(ncor_cool*1.);
    double del2 = (x_square_cool[i] / std::pow((ncor_cool*1.),2)) - (std::pow(x_cool[i],2) / (ncor_cool*1.));
    if (del2 < 0.)del2 = 0.;
    x_err_cool[i] = std::sqrt(del2);
  }
  for (size_t i = 0; i != x2_err_cool.size(); ++i){
    x2_cool[i] = x2_cool[i]/(ncor_cool*1.);
    double del2 = (x2_square_cool[i] / std::pow((ncor_cool*1.),2)) - (std::pow(x2_cool[i],2) / (ncor_cool*1.));
    if (del2 < 0.)del2 = 0.;
    x2_err_cool[i] = std::sqrt(del2);
  }
  for (size_t i = 0; i != x3_err_cool.size(); ++i){
    x3_cool[i] = x3_cool[i]/(ncor_cool*1.);
    double del2 = (x3_square_cool[i] / std::pow((ncor_cool*1.),2)) - (std::pow(x3_cool[i],2) / (ncor_cool*1.));
    if (del2 < 0.)del2 = 0.;
    x3_err_cool[i] = std::sqrt(del2);
  }
}

template <class Correlator>
void print_corr_and_log (std::ofstream& ofile, Correlator& x, Correlator& x_err, double spacing) {
  for (size_t i = 0; i != x.size()-1; ++i) {
    double dx = (x[i] - x[i+1]) / x[i] / spacing;
    double dxe2 = std::pow(x_err[i+1]/x[i],2) + std::pow((x_err[i]*x[i+1]) / std::pow(x[i],2),2);
    double dxe = std::sqrt(dxe2) / spacing;
    ofile << i * spacing <<"  "<< x[i] <<"  "<< x_err[i] <<"  "<< dx <<"  "<< dxe << '\n';
  }
}

template <class Correlator>
void print_corr_and_log_subtracted (std::ofstream& ofile, Correlator& x, Correlator& x_err, double spacing) {
  for (size_t i = 0; i != x.size()-1; ++i) {
    double x_sub = x[i] - x.back();
    double x_sub_next = x[i+1] - x.back();
    double error_sub_next = std::sqrt(std::pow(x_err[i+1],2) + std::pow(x_err.back(),2));
    double error_sub = std::sqrt(std::pow(x_err[i],2) + std::pow(x_err.back(),2));
    double dx = (x_sub - x_sub_next) / x_sub / spacing;
    double dxe2 = std::pow(error_sub_next / x_sub,2) + std::pow((error_sub*x_sub_next) / std::pow(x_sub,2),2);
    double dxe = std::sqrt(dxe2)/spacing;
    ofile << i * spacing <<"  "<< x[i] <<"  "<< x_err[i] <<"  "<< dx <<"  "<<  dxe << '\n';
  }

}

template <class Array>
void print_probability_density (std::ofstream& ofile, Array& h, int nbins, double max_val, double bin_w) {
  double xnorm;
  for(int i=0; i != nbins; ++i) {
    xnorm += (h[i] * bin_w);
  }
  for(int i=0; i != nbins; ++i) {
    double xx = (-max_val) + (double)i * bin_w;
    ofile << xx <<"  "<< h[i]/xnorm <<'\n';
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
  int constexpr np = 30; //number of points in which cor functions are evaluated
  //
  int constexpr nm = 400; //number measurements per sweep
  //
  auto potential = [](double x) { return std::pow(std::pow(x,2) - std::pow(1.4,2),2); };
  auto c{generate_random_array<800>(rng, -1.4, 1.4)}; // disordered (hot) start
  std::array<double,800>cool{};
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
  std::array<double,nbins> h_c{};

  evolve_using_metropolis<total_sweeps, take_every, equilibration, cooling_sweeps, take_every_cool>(rng, c, cool, potential, spacing, 
                          np, nm, nbins, max_val, h, h_c, xc, x2, x3, xc_er, x2_er, x3_er, xc_c, x2_c, x3_c, xc_er_c, x2_er_c, x3_er_c);                          
  
  std::ofstream fdata;
  fdata.open("points.txt");
  for (auto const i : c)fdata << i << '\n';
  fdata.close();
  std::ofstream cooldata;
  cooldata.open("points_cooled.txt");
  for (auto const i : cool)cooldata << i << '\n';
  cooldata.close();

  std::ofstream cordata;
  cordata.open("correlations.txt");
  print_corr_and_log<std::array<double,np>>(cordata, xc, xc_er, spacing); //log derivative xc
  cordata << '\n';
  cordata << '\n';
  print_corr_and_log_subtracted<std::array<double,np>>(cordata, x2, x2_er, spacing);//log derivative x2
  cordata << '\n';
  cordata << '\n';
  print_corr_and_log<std::array<double,np>>(cordata, x3, x3_er, spacing); //log derivative x3
  cordata.close();
  std::ofstream corcoldata;
  corcoldata.open("correlations_cooled.txt");
  print_corr_and_log<std::array<double,np>>(corcoldata, xc_c, xc_er_c, spacing); //log derivative xc
  corcoldata << '\n';
  corcoldata << '\n';
  print_corr_and_log_subtracted<std::array<double,np>>(corcoldata, x2_c, x2_er_c, spacing);//log derivative x2
  corcoldata << '\n';
  corcoldata << '\n';
  print_corr_and_log<std::array<double,np>>(corcoldata, x3_c, x3_er_c, spacing); //log derivative x3
  corcoldata.close();

  std::ofstream probability;
  probability.open("ground_state_probability.txt");
  print_probability_density<std::array<double,nbins>>(probability, h, nbins, max_val, bin_w);
  probability.close();

  std::ofstream probabilitycol;
  probabilitycol.open("ground_state_probability_cooled.txt");
  print_probability_density<std::array<double,nbins>>(probabilitycol, h_c, nbins, max_val, bin_w);
  probabilitycol.close();
}
