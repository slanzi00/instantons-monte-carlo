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
    for (double& value : random_array) {
        value = dist(rng);
    }
    return random_array;
}

template <typename Array>
constexpr double calculate_kinetic(const Array& positions, size_t i, const size_t positions_size, const double& current, double lattice_spacing) {
    size_t prev = (i + positions_size - 1) % (positions_size);
    size_t next = (i + 1) % (positions_size);
    double x_pm = (current - positions[prev]) / lattice_spacing;
    double x_pp = (positions[next] - current) / lattice_spacing;
    return 0.25 * (std::pow(x_pm, 2) + std::pow(x_pp, 2));
}

template <typename Array, typename Potential>
double calculate_action(const Array& positions, Potential& potential, const size_t& positions_size, const double& lattice_spacing) {
    double action = 0.;
    for (size_t i=0; i != positions_size; ++i) {
        const double& current = positions[i];
        double kinetic = calculate_kinetic(positions, i, positions_size, current, lattice_spacing);
        action += lattice_spacing * (kinetic + potential(current));
    }
    return action;
}

/*//could be an idea
template <typename Array>
void evolve_positions(std::mt19937& rng, Array& positions, Potential& potential, const size_t& positions_size, double lattice_spacing,
                      double probability_threshold) {
    std::normal_distribution<double> gaussian_step(0.0, 0.5);
    std::uniform_real_distribution<double> probability(0., 1.);
    for (double& position : positions) {
        double initial_action = calculate_action(positions, potential, positions_size, lattice_spacing);
        double dx = gaussian_step(rng);
        position += dx;
        double final_action = calculate_action(positions, potential, positions_size, lattice_spacing);
        double acceptance_prob = std::exp(initial_action - final_action);
        if (acceptance_prob <= 1.0 && acceptance_prob <= probability(rng)) {
            pos -= dx;
        }
    } 
}
*/
template <typename Histo>
void filling_histo(Histo& histo_array, double& current_position, const double& bin_width, const double& max_value) {
        int pos = std::max(floor((current_position + max_value)/bin_width), 0.);
        histo_array[pos] += 1.0;
}

template <typename T>
void calculate_statistics(T& dataset, T& square_dataset, T& errors_dataset, const double& number_correlator) {
    for (size_t i = 0; i != dataset.size(); ++i){ // means and errors of corr functions
    dataset[i] /= number_correlator;
    square_dataset[i] = (square_dataset[i] / std::pow(number_correlator,2)) - (std::pow(dataset[i], 2) / number_correlator);
    errors_dataset[i] = std::sqrt(std::max(square_dataset[i], 0.));
  }
}

template <int n_updates, int unconsidered_configurations, int thermalization_steps, int c_updates, 
          int unconsidered_cool, class Array, class Potential, class Histo, class Correlator>
void evolve_using_metropolis(std::mt19937& rng,  
                             Array& positions,
                             Array& positions_cool,
                             Potential const& potential,
                             double lattice_spacing,
                             double const& number_points,
                             double const& number_meas,
                             double const& bin_width,
                             double const& max_value,
                             Histo& histo_array,
                             Histo& histo_cool,
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
  std::normal_distribution<double> gaussian_step(0., 0.5);
  std::uniform_real_distribution<double> probability(0., 1.);
  const double ncor = (n_updates / unconsidered_configurations) * number_meas;
  const double number_cool_config = n_updates / (unconsidered_configurations * unconsidered_cool);
  const double ncor_cool = number_cool_config * number_meas;
  const size_t positions_size = positions.size();
  Correlator x_cor_square{};
  Correlator x2_cor_square{};
  Correlator x3_cor_square{};
  Correlator x_square_cool{};
  Correlator x2_square_cool{};
  Correlator x3_square_cool{};
  
  for (int m = 0; m != thermalization_steps + n_updates + 1; ++ m) { 
    //Monte Carlo
    for (auto& current_position : positions) {
        double initial_action = calculate_action(positions, potential, positions_size, lattice_spacing);
        double dx = gaussian_step(rng);
        current_position += dx;
        double final_action = calculate_action(positions, potential, positions_size, lattice_spacing);
        double acceptance_prob = std::exp(initial_action - final_action);
        if (acceptance_prob <= 1.0 && acceptance_prob <= probability(rng)) {
            current_position -= dx;
        }
        filling_histo(histo_array, current_position, bin_width, max_value);
    } 
       
     if (m > thermalization_steps && (m % unconsidered_configurations) == 0) 
     {

        for (size_t k = 0; k != number_meas; ++k) {   // correlation functions
            size_t start = std::floor((positions_size - number_points) * probability(rng));
            for (size_t l = 0; l != number_points; ++l) {
                size_t next = (start + l + 1) % (positions_size);
                double xcor = positions[start] * positions[next];
                x_cor[l] += xcor;
                x_cor_square[l] += std::pow(xcor,2);
                x2_cor[l] += std::pow(xcor,2);
                x2_cor_square[l] += std::pow(xcor,4);
                x3_cor[l] += std::pow(xcor,3);
                x3_cor_square[l] += std::pow(xcor,6);
            }
        }
    }//end if

    if (m > thermalization_steps && (m % (unconsidered_configurations * unconsidered_cool)) == 0)
    { 
      
      std::copy(positions.begin(), positions.end(), positions_cool.begin());

      for (int l = 0; l != c_updates; ++l) { 
        //cooling
        for (auto& current_cool : positions_cool) {
          double initial_action = calculate_action(positions_cool, potential, positions_size, lattice_spacing);
          double dx = gaussian_step(rng);
          current_cool += dx;
          double final_action = calculate_action(positions_cool, potential, positions_size, lattice_spacing);
          if (final_action >= initial_action) {
              current_cool -= dx;
          }
          filling_histo(histo_cool, current_cool, bin_width, max_value);
        } 
      }
      
      for (size_t k = 0; k != number_meas; ++k) {   
        // correlation functions after cooling 
          int start = std::floor((positions_size - number_points) * probability(rng));   
          for (size_t l = 0; l != number_points; ++l) {  
            size_t next = (start + l + 1) % (positions_size);
            double xcor = positions_cool[start] * positions_cool[next];
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
  
  calculate_statistics(x_cor, x_cor_square, x_cor_err, ncor);
  calculate_statistics(x2_cor, x2_cor_square, x2_cor_err, ncor);
  calculate_statistics(x3_cor, x3_cor_square, x3_cor_err, ncor);

  calculate_statistics(x_cool, x_square_cool, x_err_cool, ncor_cool);
  calculate_statistics(x2_cool, x2_square_cool, x2_err_cool, ncor_cool);
  calculate_statistics(x3_cool, x3_square_cool, x3_err_cool, ncor_cool);
}

template <class Correlator>
void print_corr_and_log (std::ofstream& ofile, Correlator& x, Correlator& x_err, double spacing) {
  for (size_t i = 0; i != x.size()-1; ++i) {
    double& position = x[i];
    double dx = (position - x[i+1]) / position / spacing;
    double dxe2 = std::pow(x_err[i+1] / position,2) + std::pow((x_err[i]*x[i+1]) / (position * position),2);
    double dxe = std::sqrt(dxe2) / spacing;
    ofile << i * spacing <<"  "<< position <<"  "<< x_err[i] <<"  "<< dx <<"  "<< dxe << '\n';
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
  double xnorm = 0;
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
  int constexpr total_sweeps = 2e4; //total number of monte carlo sweeps
  //
  int constexpr take_every = 2; //sweep to discard between every configuration used for measurements
  //
  int constexpr equilibration = 50; //first configurations to discard for equilibration purpose, before start
  //
  int constexpr cooling_sweeps = 660;
  //
  int constexpr take_every_cool = 30; //configurations to discard between every cooleing procedure
  //
  int constexpr np = 30; //number of points in which cor functions are evaluated
  //
  int constexpr nm = 400; //number measurements per sweep
  //
  auto potential = [](const double& x) { return std::pow(std::pow(x,2) - std::pow(1.4,2),2); };
  auto c{generate_random_array<1480>(rng, -1.4, 1.4)}; // disordered (hot) start
  std::array<double,1480>cool{};
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
  double constexpr bin_w = (2.* max_val) / static_cast<double>(nbins); //bin width
  //
  std::array<double,nbins> h{};
  std::array<double,nbins> h_c{};

  evolve_using_metropolis<total_sweeps, take_every, equilibration, cooling_sweeps, take_every_cool>(rng, c, cool, potential, spacing, 
                          np, nm, bin_w, max_val, h, h_c, xc, x2, x3, xc_er, x2_er, x3_er, xc_c, x2_c, x3_c, xc_er_c, x2_er_c, x3_er_c);                          
  
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
