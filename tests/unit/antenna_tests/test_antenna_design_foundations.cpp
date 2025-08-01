#include "src/antenna/current_distribution.hpp"
#include "src/antenna/array_factor.hpp"
#include "src/antenna/impedance_frequency.hpp"
#include "src/feedline/smith_chart.hpp"
#include "src/antenna/radiation_efficiency.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>

using namespace rvl::antenna;
using namespace rvl::feedline;

void test_current_distribution() {
    std::cout << "Testing Current Distribution Calculations..." << std::endl;

    current_distribution_d::wire_geometry wire;
    wire.length_m = 0.15;        
    wire.radius_m = 0.001;       
    wire.frequency_hz = 1e9;     
    wire.feed_position_m = 0.0;  

    auto center_current = current_distribution_d::calculate_center_fed_current(0.0, wire);
    std::cout << "  Center current magnitude: " << std::abs(center_current) << std::endl;

    auto end_current = current_distribution_d::calculate_center_fed_current(wire.length_m/2.0 - 0.001, wire);
    std::cout << "  End current magnitude: " << std::abs(end_current) << std::endl;

    auto samples = current_distribution_d::generate_current_samples(wire, 50, true);
    std::cout << "  Generated " << samples.size() << " current samples" << std::endl;

    auto nulls = current_distribution_d::find_current_nulls(samples, 0.1);
    std::cout << "  Found " << nulls.size() << " current nulls" << std::endl;

    auto monopole_current = current_distribution_d::calculate_monopole_current(wire.length_m/4.0, wire);
    std::cout << "  Monopole current at λ/8: " << std::abs(monopole_current) << std::endl;

    assert(std::abs(center_current) > std::abs(end_current));
    assert(samples.size() == 50);
    assert(std::abs(monopole_current) > 0.0);

    std::cout << "  Current distribution tests passed" << std::endl;
}

void test_array_factor() {
    std::cout << "Testing Array Factor Calculations..." << std::endl;

    array_factor_d::linear_array_params params;
    params.num_elements = 4;
    params.element_spacing_m = 0.15;  
    params.frequency_hz = 1e9;
    params.progressive_phase_rad = 0.0;
    params.beam_steering_angle_rad = 0.0;

    auto af_broadside = array_factor_d::calculate_linear_array_factor(0.0, params);
    std::cout << "  Array factor magnitude (broadside): " << std::abs(af_broadside) << std::endl;

    auto af_db = array_factor_d::calculate_array_factor_db(0.0, params);
    std::cout << "  Array factor (broadside): " << af_db << " dB" << std::endl;

    params.element_spacing_m = 0.35; 
    auto grating_lobes = array_factor_d::calculate_grating_lobe_angles(params);
    std::cout << "  Found " << grating_lobes.size() << " grating lobe angles" << std::endl;

    auto dolph_amplitudes = array_factor_d::generate_dolph_chebyshev_amplitudes(8, -25.0);
    std::cout << "  Generated Dolph-Chebyshev amplitudes for " << dolph_amplitudes.size() << " elements" << std::endl;

    array_factor_d::planar_array_params planar_params;
    planar_params.num_x_elements = 4;
    planar_params.num_y_elements = 4;
    planar_params.spacing_x_m = 0.15;
    planar_params.spacing_y_m = 0.15;
    planar_params.frequency_hz = 1e9;
    planar_params.steering_theta_rad = 0.0;
    planar_params.steering_phi_rad = 0.0;

    auto planar_af = array_factor_d::calculate_planar_array_factor(0.0, 0.0, planar_params);
    std::cout << "  Planar array factor magnitude: " << std::abs(planar_af) << std::endl;

    params.element_spacing_m = 0.15; 
    auto directivity = array_factor_d::calculate_array_directivity_db(params);
    std::cout << "  Array directivity: " << directivity << " dBi" << std::endl;

    assert(std::abs(af_broadside) > 0.0);
    assert(af_db > -100.0 && af_db < 100.0);
    assert(dolph_amplitudes.size() == 8);
    assert(directivity > 0.0);

    std::cout << "  Array factor tests passed" << std::endl;
}

void test_impedance_frequency() {
    std::cout << "Testing Impedance vs Frequency Analysis..." << std::endl;

    impedance_frequency_d::antenna_geometry geometry;
    geometry.length_m = 0.150;              
    geometry.radius_m = 0.001;              
    geometry.height_above_ground_m = 10.0;  
    geometry.ground_conductivity = 0.005;   
    geometry.ground_permittivity = 13.0;

    auto impedance_1ghz = impedance_frequency_d::calculate_arbitrary_length_dipole_impedance(1e9, geometry);
    std::cout << "  Dipole impedance at 1 GHz: " << impedance_1ghz.real() << " + j" 
              << impedance_1ghz.imag() << " Ω" << std::endl;

    auto king_impedance = impedance_frequency_d::calculate_king_dipole_impedance(1e9, geometry);
    std::cout << "  King's impedance at 1 GHz: " << king_impedance.real() << " + j" 
              << king_impedance.imag() << " Ω" << std::endl;

    geometry.length_m = 0.1; 
    auto loop_impedance = impedance_frequency_d::calculate_loop_impedance(1e9, geometry);
    std::cout << "  Small loop impedance: " << loop_impedance.real() << " + j" 
              << loop_impedance.imag() << " Ω" << std::endl;

    geometry.length_m = 0.150; 
    auto sweep_data = impedance_frequency_d::frequency_sweep(
        geometry, 0.8e9, 1.2e9, 100, "dipole");

    std::cout << "  Frequency sweep: " << sweep_data.size() << " points" << std::endl;
    std::cout << "  First point - Freq: " << sweep_data[0].frequency_hz/1e6 << " MHz, "
              << "VSWR: " << std::fixed << std::setprecision(2) << sweep_data[0].vswr << std::endl;

    auto bandwidth = impedance_frequency_d::analyze_bandwidth(sweep_data, 2.0, 50.0);
    std::cout << "  Center frequency: " << bandwidth.center_frequency_hz/1e6 << " MHz" << std::endl;
    std::cout << "  Bandwidth: " << bandwidth.bandwidth_hz/1e6 << " MHz" << std::endl;
    std::cout << "  Fractional BW: " << bandwidth.fractional_bandwidth * 100.0 << "%" << std::endl;
    std::cout << "  Q factor: " << bandwidth.q_factor << std::endl;
    std::cout << "  Resonances found: " << bandwidth.resonant_frequencies.size() << std::endl;

    assert(impedance_1ghz.real() > 0.0);
    assert(king_impedance.real() > 0.0);
    assert(sweep_data.size() == 100);

    if (bandwidth.bandwidth_hz <= 0.0) {
        std::cout << "  Warning: No frequencies found with VSWR < 2.0" << std::endl;
    }

    std::cout << "  Impedance frequency tests passed" << std::endl;
}

void test_smith_chart() {
    std::cout << "Testing Smith Chart Calculations..." << std::endl;

    std::complex<double> test_impedance(75.0, 25.0);
    auto gamma = smith_chart_d::impedance_to_reflection_coefficient(test_impedance, 50.0);
    std::cout << "  Reflection coefficient: " << gamma.real() << " + j" << gamma.imag() << std::endl;

    auto impedance_back = smith_chart_d::reflection_coefficient_to_impedance(gamma, 50.0);
    std::cout << "  Impedance back: " << impedance_back.real() << " + j" 
              << impedance_back.imag() << " Ω" << std::endl;

    auto smith_point = smith_chart_d::calculate_smith_point(test_impedance, 50.0);
    std::cout << "  VSWR: " << smith_point.vswr << std::endl;
    std::cout << "  Return loss: " << smith_point.return_loss_db << " dB" << std::endl;

    smith_chart_d::transmission_line_params line;
    line.characteristic_impedance = 50.0;
    line.electrical_length_degrees = 90.0; 
    line.velocity_factor = 0.66;
    line.loss_db_per_100m = 1.0;
    line.frequency_hz = 1e9;

    auto transformed_z = smith_chart_d::transform_impedance_through_line(test_impedance, line);
    std::cout << "  Transformed impedance (90°): " << transformed_z.real() << " + j" 
              << transformed_z.imag() << " Ω" << std::endl;

    std::complex<double> source_z(50.0, 0.0);
    std::complex<double> load_z(25.0, 10.0);
    auto l_network = smith_chart_d::design_l_network(source_z, load_z, 1e9);

    std::cout << "  L-network design:" << std::endl;
    std::cout << "    Q factor: " << l_network.q_factor << std::endl;
    std::cout << "    Bandwidth: " << l_network.bandwidth_hz/1e6 << " MHz" << std::endl;
    std::cout << "    Components: " << l_network.component_types.size() << std::endl;

    auto stub_network = smith_chart_d::design_stub_match(load_z, 50.0, 1e9);
    std::cout << "  Stub matching:" << std::endl;
    std::cout << "    Line length: " << stub_network.component_values[0].real() << "°" << std::endl;
    std::cout << "    Stub length: " << stub_network.component_values[1].real() << "°" << std::endl;

    auto vswr_circle = smith_chart_d::calculate_vswr_circle(2.0, 36);
    std::cout << "  VSWR circle points: " << vswr_circle.size() << std::endl;

    assert(std::abs(impedance_back.real() - test_impedance.real()) < 0.001);
    assert(std::abs(impedance_back.imag() - test_impedance.imag()) < 0.001);
    assert(smith_point.vswr >= 1.0);
    assert(l_network.component_types.size() >= 2);
    assert(vswr_circle.size() == 36);

    std::cout << "  Smith chart tests passed" << std::endl;
}

void test_radiation_efficiency() {
    std::cout << "Testing Radiation Efficiency Calculations..." << std::endl;

    radiation_efficiency_d::antenna_parameters antenna;
    antenna.frequency_hz = 1e9;
    antenna.physical_length_m = 0.15;
    antenna.effective_length_m = 0.15;
    antenna.radiation_resistance = 73.1;
    antenna.input_impedance_magnitude = 75.0;
    antenna.input_impedance = std::complex<double>(73.1, 42.5);

    auto copper = radiation_efficiency_d::create_conductor_material("copper");
    copper.diameter_m = 0.002; 
    std::cout << "  Copper conductivity: " << copper.conductivity_s_per_m/1e6 << " MS/m" << std::endl;

    auto air = radiation_efficiency_d::create_dielectric_material("air");
    std::cout << "  Air loss tangent: " << air.loss_tangent << std::endl;

    auto conductor_loss = radiation_efficiency_d::calculate_conductor_loss_resistance(
        antenna.frequency_hz, antenna.physical_length_m, copper);
    std::cout << "  Conductor loss resistance: " << conductor_loss << " Ω" << std::endl;

    auto dielectric_loss = radiation_efficiency_d::calculate_dielectric_loss_resistance(
        antenna.frequency_hz, air, antenna);
    std::cout << "  Dielectric loss resistance: " << dielectric_loss << " Ω" << std::endl;

    auto efficiency = radiation_efficiency_d::analyze_antenna_efficiency(
        antenna, copper, air, "dipole", 50.0);

    std::cout << "  Efficiency Analysis:" << std::endl;
    std::cout << "    Radiation efficiency: " << efficiency.radiation_efficiency * 100.0 << "%" << std::endl;
    std::cout << "    Mismatch efficiency: " << efficiency.mismatch_efficiency * 100.0 << "%" << std::endl;  
    std::cout << "    Total efficiency: " << efficiency.total_efficiency * 100.0 << "%" << std::endl;
    std::cout << "    Ohmic loss: " << efficiency.ohmic_loss_db << " dB" << std::endl;
    std::cout << "    Mismatch loss: " << efficiency.mismatch_loss_db << " dB" << std::endl;
    std::cout << "    Total gain reduction: " << efficiency.gain_reduction_db << " dB" << std::endl;

    auto aluminum = radiation_efficiency_d::create_conductor_material("aluminum");
    auto steel = radiation_efficiency_d::create_conductor_material("steel");

    std::vector<radiation_efficiency_d::conductor_properties> materials = {copper, aluminum, steel};
    auto material_comparison = radiation_efficiency_d::compare_conductor_materials(
        antenna, materials, air, "dipole");

    std::cout << "  Material comparison:" << std::endl;
    std::vector<std::string> material_names = {"Copper", "Aluminum", "Steel"};
    for (size_t i = 0; i < material_comparison.size(); ++i) {
        std::cout << "    " << material_names[i] << ": " 
                  << material_comparison[i].radiation_efficiency * 100.0 << "%" << std::endl;
    }

    auto aperture_eff = radiation_efficiency_d::calculate_aperture_efficiency(
        1.0, 0.8, 0.55, 0.9, 0.95); 
    std::cout << "  Aperture efficiency: " << aperture_eff * 100.0 << "%" << std::endl;

    assert(conductor_loss > 0.0);
    assert(efficiency.radiation_efficiency > 0.0 && efficiency.radiation_efficiency <= 1.0);
    assert(efficiency.total_efficiency > 0.0 && efficiency.total_efficiency <= 1.0);
    assert(material_comparison.size() == 3);
    assert(aperture_eff > 0.0 && aperture_eff <= 1.0);

    std::cout << "  Radiation efficiency tests passed" << std::endl;
}

void test_ham_antenna_design_examples() {
    std::cout << "Testing Ham Antenna Design Examples..." << std::endl;

    std::cout << "  20m Band Dipole Analysis:" << std::endl;

    impedance_frequency_d::antenna_geometry dipole_20m;
    dipole_20m.length_m = 10.0;              
    dipole_20m.radius_m = 0.001;             
    dipole_20m.height_above_ground_m = 15.0; 
    dipole_20m.ground_conductivity = 0.005;
    dipole_20m.ground_permittivity = 13.0;

    auto dipole_z = impedance_frequency_d::calculate_arbitrary_length_dipole_impedance(14.2e6, dipole_20m);
    std::cout << "    Impedance at 14.2 MHz: " << dipole_z.real() << " + j" << dipole_z.imag() << " Ω" << std::endl;

    auto smith_point = smith_chart_d::calculate_smith_point(dipole_z, 50.0);
    std::cout << "    VSWR (50Ω): " << smith_point.vswr << std::endl;

    std::cout << "  2m Yagi Array Analysis:" << std::endl;

    array_factor_d::linear_array_params yagi;
    yagi.num_elements = 5;              
    yagi.element_spacing_m = 0.35;      
    yagi.frequency_hz = 144.2e6;        
    yagi.progressive_phase_rad = 0.0;
    yagi.beam_steering_angle_rad = 0.0;

    auto yagi_directivity = array_factor_d::calculate_array_directivity_db(yagi);
    auto yagi_beamwidth = array_factor_d::calculate_array_beamwidth_rad(yagi);

    std::cout << "    Directivity: " << yagi_directivity << " dBi" << std::endl;
    std::cout << "    Beamwidth: " << yagi_beamwidth * 180.0 / M_PI << "°" << std::endl;

    std::cout << "  Wire Size Efficiency Comparison (40m band):" << std::endl;

    radiation_efficiency_d::antenna_parameters antenna_40m;
    antenna_40m.frequency_hz = 7.1e6;
    antenna_40m.physical_length_m = 20.0;     
    antenna_40m.radiation_resistance = 73.1;
    antenna_40m.input_impedance = std::complex<double>(73.1, 42.5);
    antenna_40m.input_impedance_magnitude = 75.0;

    std::vector<double> wire_diameters = {0.001, 0.002, 0.003, 0.005}; 
    std::vector<std::string> wire_descriptions = {"AWG 12", "AWG 10", "AWG 8", "3/16 inch"};

    auto air_dielectric = radiation_efficiency_d::create_dielectric_material("air");

    for (size_t i = 0; i < wire_diameters.size(); ++i) {
        auto conductor = radiation_efficiency_d::create_conductor_material("copper");
        conductor.diameter_m = wire_diameters[i];

        auto eff_analysis = radiation_efficiency_d::analyze_antenna_efficiency(
            antenna_40m, conductor, air_dielectric, "dipole");

        std::cout << "    " << wire_descriptions[i] << " (" << wire_diameters[i]*1000.0 << "mm): " 
                  << eff_analysis.radiation_efficiency * 100.0 << "%" << std::endl;
    }

    std::cout << "  OCF Dipole Matching Network:" << std::endl;

    std::complex<double> ocf_impedance(200.0, 100.0); 
    std::complex<double> feedline_z(50.0, 0.0);

    auto matching_network = smith_chart_d::design_l_network(feedline_z, ocf_impedance, 14.2e6);
    std::cout << "    Network Q: " << matching_network.q_factor << std::endl;
    std::cout << "    Bandwidth: " << matching_network.bandwidth_hz / 1000.0 << " kHz" << std::endl;

    assert(dipole_z.real() > 0.0);
    assert(smith_point.vswr >= 1.0);
    assert(yagi_directivity > 0.0);
    assert(matching_network.q_factor > 0.0);

    std::cout << "  Ham antenna design examples passed" << std::endl;
}

int main() {
    std::cout << "=== RadioVectorLib Antenna Design Foundations Test ===" << std::endl;
    std::cout << std::endl;

    try {
        test_current_distribution();
        std::cout << std::endl;

        test_array_factor();
        std::cout << std::endl;

        test_impedance_frequency();
        std::cout << std::endl;

        test_smith_chart();
        std::cout << std::endl;

        test_radiation_efficiency();
        std::cout << std::endl;

        test_ham_antenna_design_examples();
        std::cout << std::endl;

        std::cout << "=== ALL ANTENNA DESIGN FOUNDATION TESTS PASSED! ===" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}