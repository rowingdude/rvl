#include "src/propagation/tropospheric_scatter.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

using namespace rvl::propagation;

void test_atmospheric_conditions() {
    std::cout << "Testing Atmospheric Conditions Creation..." << std::endl;
    auto standard_atm = tropospheric_scatter_d::create_standard_atmosphere();
    std::cout << "  Standard Atmosphere (1 km height):" << std::endl;
    std::cout << "    Temperature: " << standard_atm.temperature_k << " K" << std::endl;
    std::cout << "    Pressure: " << standard_atm.pressure_hpa << " hPa" << std::endl;
    std::cout << "    Humidity: " << standard_atm.humidity_percent << "%" << std::endl;
    std::cout << "    Cn²: " << standard_atm.turbulence_strength_cn2 << " m⁻²/³" << std::endl;
    std::cout << "    Outer scale: " << standard_atm.outer_scale_length_m << " m" << std::endl;
    auto summer_atm = tropospheric_scatter_d::create_standard_atmosphere(2.0, 1.5);
    auto winter_atm = tropospheric_scatter_d::create_standard_atmosphere(2.0, 0.7);
    std::cout << "  Summer vs Winter turbulence:" << std::endl;
    std::cout << "    Summer Cn²: " << summer_atm.turbulence_strength_cn2 << std::endl;
    std::cout << "    Winter Cn²: " << winter_atm.turbulence_strength_cn2 << std::endl;
    assert(standard_atm.temperature_k > 250.0 && standard_atm.temperature_k < 300.0);
    assert(standard_atm.pressure_hpa > 800.0 && standard_atm.pressure_hpa < 1100.0);
    assert(summer_atm.turbulence_strength_cn2 > winter_atm.turbulence_strength_cn2);
    std::cout << "  Atmospheric conditions tests passed" << std::endl;
}

void test_path_geometry() {
    std::cout << "Testing Path Geometry Creation..." << std::endl;
    double tx_lat = 40.7128;
    double tx_lon = -74.0060;
    double rx_lat = 34.0522;
    double rx_lon = -118.2437;
    auto geometry = tropospheric_scatter_d::create_path_geometry(
        tx_lat, tx_lon, 30.0,
        rx_lat, rx_lon, 30.0
    );
    std::cout << "  NYC to LA Path:" << std::endl;
    std::cout << "    Distance: " << geometry.path_distance_km << " km" << std::endl;
    std::cout << "    Scatter angle: " << geometry.scattering_angle_rad * 180.0 / M_PI << "°" << std::endl;
    std::cout << "    Common volume height: " << geometry.common_volume_height_m << " m" << std::endl;
    auto short_path = tropospheric_scatter_d::create_path_geometry(
        40.0, -75.0, 50.0,
        41.0, -74.0, 50.0
    );
    std::cout << "  Short Path (~150 km):" << std::endl;
    std::cout << "    Distance: " << short_path.path_distance_km << " km" << std::endl;
    std::cout << "    Beyond horizon: " << std::boolalpha 
              << tropospheric_scatter_d::is_beyond_horizon(50.0, 50.0, short_path.path_distance_km) << std::endl;
    assert(geometry.path_distance_km > 3500.0 && geometry.path_distance_km < 4500.0);
    assert(short_path.path_distance_km > 100.0 && short_path.path_distance_km < 200.0);
    assert(geometry.scattering_angle_rad > 0.0);
    std::cout << "  Path geometry tests passed" << std::endl;
}

void test_enhanced_scatter_analysis() {
    std::cout << "Testing Enhanced Scatter Analysis..." << std::endl;
    auto geometry = tropospheric_scatter_d::create_path_geometry(
        42.0, -71.0, 100.0,
        41.0, -73.0, 100.0
    );
    auto atmosphere = tropospheric_scatter_d::create_standard_atmosphere(1.5, 1.2);
    double frequency_vhf = 144e6;
    double tx_gain = 12.0;
    double rx_gain = 10.0;
    auto analysis = tropospheric_scatter_d::calculate_enhanced_scatter_analysis(
        frequency_vhf, geometry, atmosphere, tx_gain, rx_gain);
    std::cout << "  VHF Troposcatter Analysis (144 MHz, " << geometry.path_distance_km << " km):" << std::endl;
    std::cout << "    Total path loss: " << std::setprecision(1) << std::fixed 
              << analysis.path_loss_db << " dB" << std::endl;
    std::cout << "    Basic transmission loss: " << analysis.basic_transmission_loss_db << " dB" << std::endl;
    std::cout << "    Scattering loss: " << analysis.scattering_loss_db << " dB" << std::endl;
    std::cout << "    Atmospheric absorption: " << analysis.atmospheric_absorption_db << " dB" << std::endl;
    std::cout << "    Fade margin required: " << analysis.fade_margin_db << " dB" << std::endl;
    std::cout << "    Path availability: " << analysis.availability_percent << "%" << std::endl;
    std::cout << "    Coherence bandwidth: " << analysis.coherence_bandwidth_hz / 1e3 << " kHz" << std::endl;
    std::cout << "    Coherence time: " << analysis.coherence_time_ms << " ms" << std::endl;
    std::cout << "    Signal std deviation: " << analysis.signal_std_deviation_db << " dB" << std::endl;
    double frequency_uhf = 432e6;
    auto analysis_uhf = tropospheric_scatter_d::calculate_enhanced_scatter_analysis(
        frequency_uhf, geometry, atmosphere, tx_gain, rx_gain);
    std::cout << "  UHF Comparison (432 MHz):" << std::endl;
    std::cout << "    Path loss difference: " << std::showpos 
              << (analysis_uhf.path_loss_db - analysis.path_loss_db) << " dB" << std::noshowpos << std::endl;
    std::cout << "    Coherence bandwidth: " << analysis_uhf.coherence_bandwidth_hz / 1e3 << " kHz" << std::endl;
    assert(analysis.path_loss_db > 50.0);
    assert(analysis.fade_margin_db > 10.0);
    assert(analysis.coherence_bandwidth_hz > 1000.0);
    assert(analysis_uhf.path_loss_db > analysis.path_loss_db);
    std::cout << "  Enhanced scatter analysis tests passed" << std::endl;
}

void test_availability_statistics() {
    std::cout << "Testing Availability Statistics..." << std::endl;
    auto geometry = tropospheric_scatter_d::create_path_geometry(
        40.0, -80.0, 50.0,
        42.0, -83.0, 50.0
    );
    auto atmosphere = tropospheric_scatter_d::create_standard_atmosphere(1.0, 1.0);
    auto analysis = tropospheric_scatter_d::calculate_enhanced_scatter_analysis(
        1296e6, geometry, atmosphere, 15.0, 15.0);
    auto stats = tropospheric_scatter_d::calculate_availability_statistics(analysis, 20);
    std::cout << "  Availability Statistics (1296 MHz, " << geometry.path_distance_km << " km):" << std::endl;
    std::cout << "    50% availability: " << stats.availability_50_percent_db << " dB" << std::endl;
    std::cout << "    90% availability: " << stats.availability_90_percent_db << " dB" << std::endl;
    std::cout << "    99% availability: " << stats.availability_99_percent_db << " dB" << std::endl;
    std::cout << "  Time vs Signal Level Percentiles:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), stats.time_percentiles.size()); ++i) {
        std::cout << "    " << std::setw(3) << (int)stats.time_percentiles[i] 
                  << "% time: " << std::setprecision(1) << stats.signal_levels_db[i] << " dB" << std::endl;
    }
    assert(stats.availability_99_percent_db > stats.availability_90_percent_db);
    assert(stats.availability_90_percent_db > stats.availability_50_percent_db);
    assert(stats.time_percentiles.size() == 20);
    std::cout << "  Availability statistics tests passed" << std::endl;
}

void test_diversity_calculations() {
    std::cout << "Testing Diversity Calculations..." << std::endl;
    auto geometry = tropospheric_scatter_d::create_path_geometry(
        45.0, -75.0, 100.0,
        43.0, -79.0, 100.0
    );
    auto atmosphere = tropospheric_scatter_d::create_standard_atmosphere();
    double freq1 = 1280e6;
    double freq2 = 1296e6;
    double freq3 = 2320e6;
    auto freq_diversity_close = tropospheric_scatter_d::calculate_frequency_diversity_improvement(
        freq1, freq2, geometry, atmosphere);
    auto freq_diversity_wide = tropospheric_scatter_d::calculate_frequency_diversity_improvement(
        freq1, freq3, geometry, atmosphere);
    std::cout << "  Frequency Diversity:" << std::endl;
    std::cout << "    Close spacing (16 MHz): " << freq_diversity_close << " dB improvement" << std::endl;
    std::cout << "    Wide spacing (1040 MHz): " << freq_diversity_wide << " dB improvement" << std::endl;
    double antenna_separation_close = 10.0;
    double antenna_separation_wide = 100.0;
    auto space_diversity_close = tropospheric_scatter_d::calculate_space_diversity_improvement(
        antenna_separation_close, freq1, geometry);
    auto space_diversity_wide = tropospheric_scatter_d::calculate_space_diversity_improvement(
        antenna_separation_wide, freq1, geometry);
    std::cout << "  Space Diversity (1280 MHz):" << std::endl;
    std::cout << "    10m separation: " << space_diversity_close << " dB improvement" << std::endl;
    std::cout << "    100m separation: " << space_diversity_wide << " dB improvement" << std::endl;
    assert(freq_diversity_wide > freq_diversity_close);
    assert(space_diversity_wide > space_diversity_close);
    assert(freq_diversity_wide <= 6.0);
    assert(space_diversity_wide <= 5.0);
    std::cout << "  Diversity calculation tests passed" << std::endl;
}

void test_frequency_optimization() {
    std::cout << "Testing Frequency Optimization..." << std::endl;
    auto geometry = tropospheric_scatter_d::create_path_geometry(
        39.0, -77.0, 30.0,
        41.5, -81.0, 30.0
    );
    auto atmosphere = tropospheric_scatter_d::create_standard_atmosphere(1.0, 0.8);
    std::vector<double> test_frequencies = {
        144e6,
        432e6,
        1296e6,
        2320e6,
        5760e6,
        10368e6
    };
    tropospheric_scatter_d::vector_type freq_vector(test_frequencies.begin(), test_frequencies.end());
    auto optimal_freq = tropospheric_scatter_d::find_optimal_frequency(
        freq_vector, geometry, atmosphere, 90.0);
    std::cout << "  Frequency Optimization (" << geometry.path_distance_km << " km, winter):" << std::endl;
    std::cout << "    Optimal frequency for 90% availability: " << optimal_freq / 1e6 << " MHz" << std::endl;
    auto freq_sweep = tropospheric_scatter_d::analyze_frequency_sweep(
        freq_vector, geometry, atmosphere);
    std::cout << "  Frequency Sweep Results:" << std::endl;
    for (size_t i = 0; i < test_frequencies.size(); ++i) {
        std::cout << "    " << std::setw(5) << (int)(test_frequencies[i] / 1e6) << " MHz: "
                  << std::setprecision(1) << freq_sweep[i].path_loss_db << " dB, "
                  << std::setprecision(0) << freq_sweep[i].availability_percent << "%" << std::endl;
    }
    assert(optimal_freq >= test_frequencies[0] && optimal_freq <= test_frequencies.back());
    assert(freq_sweep.size() == test_frequencies.size());
    std::cout << "  Frequency optimization tests passed" << std::endl;
}

void test_ham_radio_scenarios() {
    std::cout << "Testing Ham Radio Scenarios..." << std::endl;
    std::cout << "  Scenario 1: VHF/UHF Weak Signal Communication" << std::endl;
    auto moderate_path = tropospheric_scatter_d::create_path_geometry(
        42.36, -71.06, 100.0,
        40.71, -74.01, 100.0
    );
    auto good_conditions = tropospheric_scatter_d::create_standard_atmosphere(1.2, 1.3);
    auto vhf_analysis = tropospheric_scatter_d::calculate_enhanced_scatter_analysis(
        144e6, moderate_path, good_conditions, 17.0, 17.0);
    std::cout << "    144 MHz, " << moderate_path.path_distance_km << " km:" << std::endl;
    std::cout << "      Path loss: " << vhf_analysis.path_loss_db << " dB" << std::endl;
    std::cout << "      Required fade margin: " << vhf_analysis.fade_margin_db << " dB" << std::endl;
    std::cout << "      Link margin (1kW, -130dBm RX): " << 
                 (90.0 - vhf_analysis.path_loss_db - vhf_analysis.fade_margin_db) << " dB" << std::endl;
    std::cout << "  Scenario 2: Microwave DX (10 GHz)" << std::endl;
    auto long_path = tropospheric_scatter_d::create_path_geometry(
        47.61, -122.33, 200.0,
        45.52, -122.68, 200.0
    );
    auto microwave_analysis = tropospheric_scatter_d::calculate_enhanced_scatter_analysis(
        10368e6, long_path, good_conditions, 30.0, 30.0);
    std::cout << "    10368 MHz, " << long_path.path_distance_km << " km:" << std::endl;
    std::cout << "      Path loss: " << microwave_analysis.path_loss_db << " dB" << std::endl;
    std::cout << "      Coherence bandwidth: " << microwave_analysis.coherence_bandwidth_hz / 1e6 << " MHz" << std::endl;
    std::cout << "      Signal stability: " << microwave_analysis.signal_std_deviation_db << " dB std dev" << std::endl;
    std::cout << "  Scenario 3: Contest Operation with Diversity" << std::endl;
    auto contest_path = tropospheric_scatter_d::create_path_geometry(
        41.7, -91.5, 50.0,
        44.9, -93.3, 50.0
    );
    auto single_antenna = tropospheric_scatter_d::calculate_enhanced_scatter_analysis(
        432e6, contest_path, good_conditions, 12.0, 12.0);
    auto frequency_diversity = tropospheric_scatter_d::calculate_frequency_diversity_improvement(
        432e6, 435e6, contest_path, good_conditions);
    auto space_diversity = tropospheric_scatter_d::calculate_space_diversity_improvement(
        50.0, 432e6, contest_path);
    std::cout << "    432 MHz contest setup:" << std::endl;
    std::cout << "      Base path loss: " << single_antenna.path_loss_db << " dB" << std::endl;
    std::cout << "      Frequency diversity gain: " << frequency_diversity << " dB" << std::endl;
    std::cout << "      Space diversity gain: " << space_diversity << " dB" << std::endl;
    std::cout << "      Combined improvement: " << (frequency_diversity + space_diversity) << " dB" << std::endl;
    assert(vhf_analysis.path_loss_db > 70.0);
    assert(microwave_analysis.path_loss_db > vhf_analysis.path_loss_db);
    assert(frequency_diversity > 0.0);
    std::cout << "  Ham radio scenario tests passed" << std::endl;
}

int main() {
    std::cout << "=== Enhanced Tropospheric Scatter Test Suite ===" << std::endl;
    std::cout << std::endl;
    try {
        test_atmospheric_conditions();
        std::cout << std::endl;
        test_path_geometry();
        std::cout << std::endl;
        test_enhanced_scatter_analysis();
        std::cout << std::endl;
        test_availability_statistics();
        std::cout << std::endl;
        test_diversity_calculations();
        std::cout << std::endl;
        test_frequency_optimization();
        std::cout << std::endl;
        test_ham_radio_scenarios();
        std::cout << std::endl;
        std::cout << "=== ALL TROPOSPHERIC SCATTER TESTS PASSED! ===" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}