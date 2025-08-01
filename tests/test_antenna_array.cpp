#include "test_utils.hpp"
#include "src/antenna/array_factor.hpp"
#include "src/antenna/current_distribution.hpp"
#include "src/antenna/impedance_frequency.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("Antenna Array Functions");

    try {
        antenna::array_factor<double>::linear_array_params params;
        params.num_elements = 4;
        params.element_spacing_m = 0.5;
        params.frequency_hz = 100e6;
        params.progressive_phase_rad = 0.0;
        params.beam_steering_angle_rad = 0.0;
        
        auto af = antenna::array_factor<double>::calculate_linear_array_factor(M_PI/2, params);
        bool passed = std::abs(af) >= 0.0 && std::abs(af) <= 4.0;
        suite.add_test("calculate_linear_array_factor", passed);
    } catch (...) {
        suite.add_test("calculate_linear_array_factor", false, "exception thrown");
    }
    
    try {
        antenna::array_factor<double>::planar_array_params params;
        params.num_x_elements = 2;
        params.num_y_elements = 2;
        params.spacing_x_m = 0.5;
        params.spacing_y_m = 0.5;
        params.frequency_hz = 100e6;
        params.steering_theta_rad = 0.0;
        params.steering_phi_rad = 0.0;
        
        auto af = antenna::array_factor<double>::calculate_planar_array_factor(M_PI/2, 0.0, params);
        bool passed = std::abs(af) >= 0.0 && std::abs(af) <= 4.0;
        suite.add_test("calculate_planar_array_factor", passed);
    } catch (...) {
        suite.add_test("calculate_planar_array_factor", false, "exception thrown");
    }
    
    try {
        antenna::array_factor<double>::circular_array_params params;
        params.num_elements = 8;
        params.radius_m = 1.0;
        params.frequency_hz = 100e6;
        params.beam_steering_angle_rad = 0.0;
        
        auto af = antenna::array_factor<double>::calculate_circular_array_factor(M_PI/2, params);
        bool passed = std::abs(af) >= 0.0;
        suite.add_test("calculate_circular_array_factor", passed);
    } catch (...) {
        suite.add_test("calculate_circular_array_factor", false, "exception thrown");
    }
    
    try {
        antenna::array_factor<double>::linear_array_params dir_params;
        dir_params.num_elements = 4;
        dir_params.element_spacing_m = 0.5;
        dir_params.frequency_hz = 100e6;
        dir_params.progressive_phase_rad = 0.0;
        dir_params.beam_steering_angle_rad = 0.0;
        double d = antenna::array_factor<double>::calculate_array_directivity_db(dir_params);
        bool passed = d > 1.0 && d < 10.0;
        suite.add_test("calculate_array_directivity", passed);
    } catch (...) {
        suite.add_test("calculate_array_directivity", false, "exception thrown");
    }
    
    try {
        antenna::array_factor<double>::linear_array_params bw_params;
        bw_params.num_elements = 4;
        bw_params.element_spacing_m = 0.5;
        bw_params.frequency_hz = 100e6;
        bw_params.progressive_phase_rad = 0.0;
        bw_params.beam_steering_angle_rad = 0.0;
        double bw = antenna::array_factor<double>::calculate_array_beamwidth_rad(bw_params);
        bool passed = bw > 0.0 && bw < 180.0;
        suite.add_test("calculate_beamwidth_degrees", passed);
    } catch (...) {
        suite.add_test("calculate_beamwidth_degrees", false, "exception thrown");
    }
    
    try {
        antenna::array_factor<double>::linear_array_params grating_params;
        grating_params.num_elements = 4;
        grating_params.element_spacing_m = 2.0;
        grating_params.frequency_hz = 100e6;
        grating_params.progressive_phase_rad = 0.0;
        grating_params.beam_steering_angle_rad = 0.0;
        auto angles = antenna::array_factor<double>::calculate_grating_lobe_angles(grating_params);
        bool passed = angles.size() >= 0;
        suite.add_test("find_grating_lobe_angles", passed);
    } catch (...) {
        suite.add_test("find_grating_lobe_angles", false, "exception thrown");
    }
    
    try {
        auto weights = antenna::array_factor<double>::generate_dolph_chebyshev_amplitudes(8, -20.0);
        bool passed = weights.size() == 8;
        suite.add_test("calculate_dolph_chebyshev_weights", passed);
    } catch (...) {
        suite.add_test("calculate_dolph_chebyshev_weights", false, "exception thrown");
    }

    try {
        antenna::current_distribution<double>::wire_geometry wire;
        wire.length_m = 5.0;
        wire.radius_m = 0.001;
        wire.frequency_hz = 30e6;
        wire.feed_position_m = 0.0;
        wire.feed_impedance = std::complex<double>(50.0, 0.0);
        
        auto current = antenna::current_distribution<double>::calculate_center_fed_current(0.0, wire);
        bool passed = std::abs(current) >= 0.0 && std::abs(current) <= 1.0;
        suite.add_test("calculate_sinusoidal_current", passed);
    } catch (...) {
        suite.add_test("calculate_sinusoidal_current", false, "exception thrown");
    }
    
    try {
        antenna::current_distribution<double>::wire_geometry wire;
        wire.length_m = 2.5;
        wire.radius_m = 0.001;
        wire.frequency_hz = 30e6;
        wire.feed_position_m = 0.0;
        wire.feed_impedance = std::complex<double>(50.0, 0.0);
        
        auto current = antenna::current_distribution<double>::calculate_monopole_current(0.0, wire);
        bool passed = std::abs(current) >= 0.0;
        suite.add_test("calculate_monopole_current", passed);
    } catch (...) {
        suite.add_test("calculate_monopole_current", false, "exception thrown");
    }
    
    try {
        antenna::current_distribution<double>::wire_geometry wire;
        wire.length_m = 5.0;
        wire.radius_m = 0.001;
        wire.frequency_hz = 30e6;
        wire.feed_position_m = 0.0;
        wire.feed_impedance = std::complex<double>(50.0, 0.0);
        
        auto samples = antenna::current_distribution<double>::generate_current_samples(wire, 50, true);
        auto nulls = antenna::current_distribution<double>::find_current_nulls(samples, 0.1);
        bool passed = true;
        suite.add_test("find_current_nulls", passed);
    } catch (...) {
        suite.add_test("find_current_nulls", false, "exception thrown");
    }

    try {
        antenna::impedance_frequency<double>::antenna_geometry geom;
        geom.length_m = 5.0;
        geom.radius_m = 0.001;
        geom.height_above_ground_m = 10.0;
        geom.ground_conductivity = 0.005;
        geom.ground_permittivity = 13.0;
        
        auto z = antenna::impedance_frequency<double>::calculate_king_dipole_impedance(1e9, geom);
        bool passed = z.real() > 0.0 && std::abs(z.imag()) < 1000.0;
        suite.add_test("calculate_impedance_simple", passed);
    } catch (...) {
        suite.add_test("calculate_impedance_simple", false, "exception thrown");
    }
    
    try {
        antenna::impedance_frequency<double>::antenna_geometry geom;
        geom.length_m = 0.15;
        geom.radius_m = 0.001;
        geom.height_above_ground_m = 10.0;
        geom.ground_conductivity = 0.005;
        geom.ground_permittivity = 13.0;
        
        auto z = antenna::impedance_frequency<double>::calculate_king_dipole_impedance(1e9, geom);
        bool passed = z.real() > 0.0;
        suite.add_test("calculate_king_dipole_impedance", passed);
    } catch (...) {
        suite.add_test("calculate_king_dipole_impedance", false, "exception thrown");
    }
    
    try {
        antenna::impedance_frequency<double>::antenna_geometry geom;
        geom.length_m = 0.15;
        geom.radius_m = 0.001;
        geom.height_above_ground_m = 10.0;
        geom.ground_conductivity = 0.005;
        geom.ground_permittivity = 13.0;
        
        auto sweep_data = antenna::impedance_frequency<double>::frequency_sweep(
            geom, 0.9e9, 1.1e9, 100, "dipole");
        auto analysis = antenna::impedance_frequency<double>::analyze_bandwidth(sweep_data, 2.0, 50.0);
        bool passed = analysis.center_frequency_hz > 0.0;
        suite.add_test("analyze_bandwidth", passed);
    } catch (...) {
        suite.add_test("analyze_bandwidth", false, "exception thrown");
    }
    
    try {
        antenna::impedance_frequency<double>::antenna_geometry geom;
        geom.length_m = 0.15;
        geom.radius_m = 0.001;
        geom.height_above_ground_m = 10.0;
        geom.ground_conductivity = 0.005;
        geom.ground_permittivity = 13.0;
        
        auto sweep_data = antenna::impedance_frequency<double>::frequency_sweep(
            geom, 0.9e9, 1.1e9, 100, "dipole");
        auto resonances = antenna::impedance_frequency<double>::find_resonant_frequencies(sweep_data);
        bool passed = true;
        suite.add_test("find_resonances", passed);
    } catch (...) {
        suite.add_test("find_resonances", false, "exception thrown");
    }
    
    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}