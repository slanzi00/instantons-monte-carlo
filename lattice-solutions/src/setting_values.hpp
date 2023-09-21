#ifndef SETTING_VALUES_HPP
#define SETTING_VALUES_HPP

namespace sv {

constexpr auto n_lattice_points = 800;
constexpr auto lattice_spacing = 0.05;
constexpr auto n_correlator_points = 30;
constexpr auto n_correlator_meas = 400;
constexpr auto take_corr_meas_every = 2;
constexpr auto n_histogram_bins = 500;
constexpr auto n_therm_sweeps = 100;
constexpr auto n_sweeps = 1e4;
constexpr auto n_sweeps_cool = 5000;
constexpr auto make_cool_every = 1000;
constexpr auto x_min_histogram = -3.;
constexpr auto x_max_histogram = 3.;

};  // namespace sv

#endif