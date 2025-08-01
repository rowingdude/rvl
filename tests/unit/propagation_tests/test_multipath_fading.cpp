#include "src/propagation/two_ray_model.hpp"
#include "src/propagation/rayleigh_fading.hpp"
#include "src/propagation/rician_fading.hpp"
#include "src/propagation/deterministic_ray_tracing.hpp"
#include "src/rf_systems/link_budget.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <random>

using namespace rvl::propagation;
using namespace rvl::rf_systems;

void test_two_ray_model() {
    std::cout << "Testing Two-Ray Ground Reflection Model..." << std::endl;

    const double frequency = 900e6; 
    const double distance = 1000.0; 
    const double tx_height = 30.0; 
    const double rx_height = 1.5;  

    two_ray_model_d::two_ray_geometry geometry;
    geometry.distance_m = distance;
    geometry.operating_frequency_hz = frequency;
    geometry.transmitter_height_m = tx_height;
    geometry.receiver_height_m = rx_height;

    two_ray_model_d::ground_reflection_parameters ground;
    ground.relative_permittivity = 15.0;
    ground.conductivity_s_per_m = 0.01;
    ground.surface_roughness_factor = 0.02;

    auto path_loss_db = two_ray_model_d::calculate_two_ray_path_loss_db(
        geometry, ground, true);

    std::cout << "  Two-ray path loss: " << std::fixed << std::setprecision(2) 
              << path_loss_db << " dB" << std::endl;

    const double grazing_angle = two_ray_model_d::calculate_grazing_angle_rad(geometry);

    auto refl_coeff_h = two_ray_model_d::calculate_ground_reflection_coefficient(
        grazing_angle, frequency, ground, true);
    auto refl_coeff_v = two_ray_model_d::calculate_ground_reflection_coefficient(
        grazing_angle, frequency, ground, false);

    std::cout << "  Horizontal reflection coeff: " << std::abs(refl_coeff_h) 
              << " ∠" << std::arg(refl_coeff_h) * 180.0 / M_PI << "°" << std::endl;
    std::cout << "  Vertical reflection coeff: " << std::abs(refl_coeff_v) 
              << " ∠" << std::arg(refl_coeff_v) * 180.0 / M_PI << "°" << std::endl;

    assert(path_loss_db > 0.0 && path_loss_db < 200.0);
    assert(std::abs(refl_coeff_h) <= 1.0);
    assert(std::abs(refl_coeff_v) <= 1.0);

    std::cout << "  Two-ray model tests passed" << std::endl;
}

void test_rayleigh_fading() {
    std::cout << "Testing Rayleigh Fading Model..." << std::endl;

    const double sigma = 1.0;
    const double envelope = 1.5;

    auto pdf = rayleigh_fading_d::calculate_rayleigh_pdf(envelope, sigma);
    auto cdf = rayleigh_fading_d::calculate_rayleigh_cdf(envelope, sigma);

    std::cout << "  PDF(1.5, σ=1.0): " << std::fixed << std::setprecision(4) 
              << pdf << std::endl;
    std::cout << "  CDF(1.5, σ=1.0): " << std::fixed << std::setprecision(4) 
              << cdf << std::endl;

    auto mean = rayleigh_fading_d::calculate_rayleigh_mean(sigma);
    auto variance = rayleigh_fading_d::calculate_rayleigh_variance(sigma);

    std::cout << "  Mean: " << mean << ", Variance: " << variance << std::endl;

    auto outage_prob = rayleigh_fading_d::calculate_outage_probability(
        -10.0); 

    std::cout << "  Outage probability (10 dB below avg): " << outage_prob << std::endl;

    std::mt19937 rng(42);

    auto samples = rayleigh_fading_d::generate_rayleigh_fading_sequence(100, sigma, rng);

    std::cout << "  Generated " << samples.size() << " fading samples" << std::endl;

    assert(pdf > 0.0 && pdf < 1.0);
    assert(cdf > 0.0 && cdf < 1.0);
    assert(mean > 0.0);
    assert(variance > 0.0);
    assert(samples.size() == 100);

    std::cout << "  Rayleigh fading tests passed" << std::endl;
}

void test_rician_fading() {
    std::cout << "Testing Rician Fading Model..." << std::endl;

    const double k_db = 10.0;
    const double k_linear = rician_fading_d::calculate_k_factor_linear(k_db);
    const double k_db_back = rician_fading_d::calculate_k_factor_db(k_linear);

    std::cout << "  K-factor: " << k_db << " dB = " << k_linear << " linear" << std::endl;
    assert(std::abs(k_db - k_db_back) < 1e-10);

    const double envelope = 2.0;
    const double omega = 4.0; 

    auto pdf = rician_fading_d::calculate_rician_pdf(envelope, k_linear, omega);
    auto cdf = rician_fading_d::calculate_rician_cdf(envelope, k_linear, omega);

    std::cout << "  PDF: " << pdf << ", CDF: " << cdf << std::endl;

    auto mean = rician_fading_d::calculate_rician_mean(k_linear, omega);
    auto variance = rician_fading_d::calculate_rician_variance(k_linear, omega);

    std::cout << "  Mean: " << mean << ", Variance: " << variance << std::endl;

    std::mt19937 rng(42);
    rician_fading_d::rician_parameters params;
    params.k_factor_db = k_db;
    params.total_power_db = 10.0 * std::log10(omega);
    params.doppler_frequency_hz = 50.0;
    params.phase_los_rad = 0.0;

    auto samples = rician_fading_d::generate_rician_fading_sequence(50, params, rng);

    std::cout << "  Generated " << samples.size() << " Rician samples" << std::endl;

    auto estimated_params = rician_fading_d::estimate_rician_parameters(samples);
    std::cout << "  Estimated K-factor: " << estimated_params.k_factor_db << " dB" << std::endl;

    assert(pdf > 0.0);
    assert(cdf >= 0.0 && cdf <= 1.0);
    assert(mean > 0.0);
    assert(samples.size() == 50);

    std::cout << "  Rician fading tests passed" << std::endl;
}

void test_deterministic_ray_tracing() {
    std::cout << "Testing Deterministic Ray Tracing..." << std::endl;

    std::array<double, 3> tx_pos = {0.0, 0.0, 10.0}; 
    std::array<double, 3> rx_pos = {100.0, 50.0, 1.5}; 

    const double frequency = 2.4e9; 

    auto direct_path = deterministic_ray_tracing_d::trace_direct_path(tx_pos, rx_pos, frequency);

    std::cout << "  Direct path length: " << direct_path.total_path_length << " m" << std::endl;
    std::cout << "  Direct path delay: " << direct_path.total_delay_s * 1e9 << " ns" << std::endl;
    std::cout << "  Direct path gain: " << 20.0 * std::log10(std::abs(direct_path.total_path_gain)) << " dB" << std::endl;

    deterministic_ray_tracing_d::geometric_obstacle building;
    building.center = {50.0, 25.0, 5.0};
    building.dimensions = {20.0, 10.0, 15.0}; 
    building.material = deterministic_ray_tracing_d::create_material_properties("concrete");
    building.obstacle_type = 1; 

    auto reflected_path = deterministic_ray_tracing_d::trace_single_reflection_path(
        tx_pos, rx_pos, building, frequency);

    std::cout << "  Reflected path length: " << reflected_path.total_path_length << " m" << std::endl;
    std::cout << "  Reflected path gain: " << 20.0 * std::log10(std::abs(reflected_path.total_path_gain)) << " dB" << std::endl;

    std::vector<deterministic_ray_tracing_d::geometric_obstacle> obstacles = {building};
    auto multipath_profile = deterministic_ray_tracing_d::calculate_multipath_profile(
        tx_pos, rx_pos, obstacles, frequency, 2);

    std::cout << "  Found " << multipath_profile.size() << " propagation paths" << std::endl;

    auto rms_delay_spread = deterministic_ray_tracing_d::calculate_rms_delay_spread_s(multipath_profile);
    std::cout << "  RMS delay spread: " << rms_delay_spread * 1e9 << " ns" << std::endl;

    auto received_power = deterministic_ray_tracing_d::calculate_received_power_dbm(
        20.0, 10.0, 5.0, multipath_profile); 

    std::cout << "  Received power: " << received_power << " dBm" << std::endl;

    assert(direct_path.total_path_length > 0.0);
    assert(direct_path.total_delay_s > 0.0);
    assert(multipath_profile.size() >= 1); 
    assert(rms_delay_spread >= 0.0);

    std::cout << "  Ray tracing tests passed" << std::endl;
}

void test_link_budget() {
    std::cout << "Testing Link Budget Analysis..." << std::endl;

    auto tx = link_budget_d::create_typical_transmitter("cellular_base_station");
    auto rx = link_budget_d::create_typical_receiver("cellular_mobile");

    std::cout << "  Transmitter EIRP: " << tx.eirp_dbm << " dBm" << std::endl;
    std::cout << "  Receiver sensitivity: " << rx.sensitivity_dbm << " dBm" << std::endl;

    link_budget_d::path_parameters path;
    path.frequency_hz = 1900e6; 
    path.distance_m = 2000.0;   

    const double wavelength = 3e8 / path.frequency_hz;
    path.path_loss_db = 20.0 * std::log10(4.0 * M_PI * path.distance_m / wavelength);
    path.atmospheric_loss_db = 0.2; 
    path.rain_loss_db = 1.0;
    path.fade_margin_db = 10.0;
    path.polarization_loss_db = 0.5;

    std::cout << "  Free space path loss: " << path.path_loss_db << " dB" << std::endl;

    const double bandwidth = 5e6; 
    const double required_snr = 10.0; 

    auto result = link_budget_d::analyze_link(tx, rx, path, bandwidth, required_snr);

    std::cout << "  Received power: " << result.received_power_dbm << " dBm" << std::endl;
    std::cout << "  Noise floor: " << result.noise_floor_dbm << " dBm" << std::endl;
    std::cout << "  SNR: " << result.snr_db << " dB" << std::endl;
    std::cout << "  Link margin: " << result.link_margin_db << " dB" << std::endl;
    std::cout << "  Link viable: " << (result.link_viable ? "Yes" : "No") << std::endl;
    std::cout << "  Availability: " << result.availability_percent << "%" << std::endl;

    auto max_range = link_budget_d::calculate_maximum_range_m(
        tx, rx, path.frequency_hz, bandwidth, required_snr);

    std::cout << "  Maximum range: " << max_range / 1000.0 << " km" << std::endl;

    auto data_rate = link_budget_d::calculate_data_rate_bps(bandwidth, result.snr_db, "QPSK");
    std::cout << "  Data rate (QPSK): " << data_rate / 1e6 << " Mbps" << std::endl;

    auto system_noise_temp = link_budget_d::calculate_system_noise_temperature_k(290.0, rx.noise_figure_db);
    auto g_over_t = link_budget_d::calculate_g_over_t_db_per_k(rx.antenna_gain_db, system_noise_temp);

    std::cout << "  System noise temperature: " << system_noise_temp << " K" << std::endl;
    std::cout << "  G/T: " << g_over_t << " dB/K" << std::endl;

    assert(result.received_power_dbm > -150.0 && result.received_power_dbm < 50.0);
    assert(result.noise_floor_dbm < -50.0);
    assert(max_range > 1000.0);
    assert(data_rate > 1e6);
    assert(system_noise_temp > 290.0);

    std::cout << "  Link budget tests passed" << std::endl;
}

void test_frequency_sweep() {
    std::cout << "Testing Frequency Sweep Analysis..." << std::endl;

    auto tx = link_budget_d::create_typical_transmitter("wifi_ap");
    auto rx = link_budget_d::create_typical_receiver("wifi_client");

    link_budget_d::vector_type frequencies = {2.4e9, 3.5e9, 5.0e9, 5.8e9};

    auto results = link_budget_d::analyze_frequency_sweep(
        tx, rx, frequencies, 100.0, 20e6); 

    std::cout << "  Frequency sweep results:" << std::endl;
    for (size_t i = 0; i < frequencies.size(); ++i) {
        std::cout << "    " << frequencies[i] / 1e9 << " GHz: " 
                  << results[i].received_power_dbm << " dBm, "
                  << results[i].snr_db << " dB SNR" << std::endl;
    }

    assert(results.size() == frequencies.size());

    std::cout << "  Frequency sweep tests passed" << std::endl;
}

void test_ham_radio_equipment() {
    std::cout << "Testing Ham Radio Equipment Profiles..." << std::endl;

    auto tx_1w = link_budget_d::create_typical_transmitter("ham_handheld_1w");
    auto tx_100w = link_budget_d::create_typical_transmitter("ham_mobile_100w");
    auto tx_1kw = link_budget_d::create_typical_transmitter("ham_base_1kw");
    auto tx_2kw = link_budget_d::create_typical_transmitter("ham_legal_limit_2kw");

    std::cout << "  Ham equipment EIRP values:" << std::endl;
    std::cout << "    1W handheld: " << tx_1w.eirp_dbm << " dBm" << std::endl;
    std::cout << "    100W mobile: " << tx_100w.eirp_dbm << " dBm" << std::endl;
    std::cout << "    1kW base: " << tx_1kw.eirp_dbm << " dBm" << std::endl;
    std::cout << "    2kW legal limit: " << tx_2kw.eirp_dbm << " dBm" << std::endl;

    auto rx_handheld = link_budget_d::create_typical_receiver("ham_handheld");
    auto rx_mobile = link_budget_d::create_typical_receiver("ham_mobile");
    auto rx_base_hf = link_budget_d::create_typical_receiver("ham_base_hf");
    auto rx_contest = link_budget_d::create_typical_receiver("ham_contest_station");

    std::cout << "  Ham receiver sensitivities:" << std::endl;
    std::cout << "    Handheld: " << rx_handheld.sensitivity_dbm << " dBm" << std::endl;
    std::cout << "    Mobile: " << rx_mobile.sensitivity_dbm << " dBm" << std::endl;
    std::cout << "    HF base: " << rx_base_hf.sensitivity_dbm << " dBm" << std::endl;
    std::cout << "    Contest station: " << rx_contest.sensitivity_dbm << " dBm" << std::endl;

    const double hf_freq = 14.2e6; 
    const double vhf_freq = 144.2e6; 
    const double uhf_freq = 432.1e6; 

    auto hf_noise = link_budget_d::calculate_ham_band_noise_floor_dbm(hf_freq, 2500.0);
    auto vhf_noise = link_budget_d::calculate_ham_band_noise_floor_dbm(vhf_freq, 2500.0);
    auto uhf_noise = link_budget_d::calculate_ham_band_noise_floor_dbm(uhf_freq, 2500.0);

    std::cout << "  Band-specific noise floors (2.5 kHz BW):" << std::endl;
    std::cout << "    20m (14.2 MHz): " << hf_noise << " dBm" << std::endl;
    std::cout << "    2m (144.2 MHz): " << vhf_noise << " dBm" << std::endl;
    std::cout << "    70cm (432.1 MHz): " << uhf_noise << " dBm" << std::endl;

    auto cw_snr = link_budget_d::calculate_required_snr_for_ham_mode("CW");
    auto ssb_snr = link_budget_d::calculate_required_snr_for_ham_mode("SSB");
    auto psk31_snr = link_budget_d::calculate_required_snr_for_ham_mode("PSK31");
    auto ft8_snr = link_budget_d::calculate_required_snr_for_ham_mode("FT8");

    std::cout << "  Required SNR for ham modes:" << std::endl;
    std::cout << "    CW: " << cw_snr << " dB" << std::endl;
    std::cout << "    SSB: " << ssb_snr << " dB" << std::endl;
    std::cout << "    PSK31: " << psk31_snr << " dB" << std::endl;
    std::cout << "    FT8: " << ft8_snr << " dB" << std::endl;

    link_budget_d::path_parameters dx_path;
    dx_path.frequency_hz = hf_freq;
    dx_path.distance_m = 10000e3; 
    dx_path.path_loss_db = 160.0;  
    dx_path.atmospheric_loss_db = 0;
    dx_path.rain_loss_db = 0;
    dx_path.fade_margin_db = 20.0; 
    dx_path.polarization_loss_db = 3.0; 

    auto dx_result = link_budget_d::analyze_link(tx_1kw, rx_contest, dx_path, 2500.0, cw_snr);

    std::cout << "  HF DX link analysis (1kW to contest station, 10,000 km):" << std::endl;
    std::cout << "    Received power: " << dx_result.received_power_dbm << " dBm" << std::endl;
    std::cout << "    SNR: " << dx_result.snr_db << " dB" << std::endl;
    std::cout << "    Link viable: " << (dx_result.link_viable ? "Yes" : "No") << std::endl;

    assert(tx_1w.eirp_dbm < tx_100w.eirp_dbm);
    assert(tx_100w.eirp_dbm < tx_1kw.eirp_dbm);
    assert(tx_1kw.eirp_dbm < tx_2kw.eirp_dbm);
    assert(rx_contest.sensitivity_dbm < rx_base_hf.sensitivity_dbm);
    assert(hf_noise > vhf_noise); 
    assert(ft8_snr < cw_snr); 

    std::cout << "  Ham radio equipment tests passed" << std::endl;
}

int main() {
    std::cout << "=== RadioVectorLib Multipath and Fading Equations Test ===" << std::endl;
    std::cout << std::endl;

    try {
        test_two_ray_model();
        std::cout << std::endl;

        test_rayleigh_fading();
        std::cout << std::endl;

        test_rician_fading();
        std::cout << std::endl;

        test_deterministic_ray_tracing();
        std::cout << std::endl;

        test_link_budget();
        std::cout << std::endl;

        test_frequency_sweep();
        std::cout << std::endl;

        test_ham_radio_equipment();
        std::cout << std::endl;

        std::cout << "=== ALL MULTIPATH/FADING TESTS PASSED! ===" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}