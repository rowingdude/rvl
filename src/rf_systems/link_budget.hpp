/*
================================================================================
||                                                                            ||
||                              RadioVectorLibrary                            ||
||                                   (RVL)                                    ||
||                                                                            ||
||                                                                            ||
||                         Developer: Benjamin Cance, KC8BWS                  ||
||                         Email: kc8bws@kc8bws.com                           ||
||                         MIT License                                        ||
||                                                                            ||
||                                                                            ||
================================================================================
*/

#ifndef RVL_RF_SYSTEMS_LINK_BUDGET_HPP
#define RVL_RF_SYSTEMS_LINK_BUDGET_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <vector>
#include <string>
#include <unordered_map>

namespace rvl {
namespace rf_systems {

template<typename T>
class link_budget {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    struct transmitter_parameters {
        T power_dbm;
        T antenna_gain_db;
        T feedline_loss_db;
        T misc_losses_db;
        T eirp_dbm;
    };

    struct receiver_parameters {
        T antenna_gain_db;
        T feedline_loss_db;
        T noise_figure_db;
        T sensitivity_dbm;
        T misc_losses_db;
    };

    struct path_parameters {
        T frequency_hz;
        T distance_m;
        T path_loss_db;
        T atmospheric_loss_db;
        T rain_loss_db;
        T fade_margin_db;
        T polarization_loss_db;
    };

    struct link_analysis_result {
        T received_power_dbm;
        T noise_floor_dbm;
        T snr_db;
        T link_margin_db;
        T availability_percent;
        bool link_viable;
        std::vector<std::string> limiting_factors;
    };

    static T calculate_eirp_dbm(const transmitter_parameters& tx) {
        return tx.power_dbm + tx.antenna_gain_db - tx.feedline_loss_db - tx.misc_losses_db;
    }

    static T calculate_noise_floor_dbm(T bandwidth_hz, T noise_figure_db = T(0), T temperature_k = T(290)) {
        core::check_positive(bandwidth_hz, "Bandwidth");
        core::check_positive(temperature_k, "Temperature");

        const T k_boltzmann = constants::physical<T>::k_B;
        const T thermal_noise_power_w = k_boltzmann * temperature_k * bandwidth_hz;
        const T thermal_noise_dbm = T(10.0) * std::log10(thermal_noise_power_w * T(1000.0));

        return thermal_noise_dbm + noise_figure_db;
    }

    static T calculate_received_power_dbm(const transmitter_parameters& tx,
                                        const receiver_parameters& rx,
                                        const path_parameters& path) {
        const T eirp = calculate_eirp_dbm(tx);

        return eirp + rx.antenna_gain_db - rx.feedline_loss_db - rx.misc_losses_db -
               path.path_loss_db - path.atmospheric_loss_db - path.rain_loss_db - 
               path.polarization_loss_db;
    }

    static link_analysis_result analyze_link(const transmitter_parameters& tx,
                                           const receiver_parameters& rx,
                                           const path_parameters& path,
                                           T bandwidth_hz,
                                           T required_snr_db = T(10.0)) {
        link_analysis_result result;

        result.received_power_dbm = calculate_received_power_dbm(tx, rx, path);
        result.noise_floor_dbm = calculate_noise_floor_dbm(bandwidth_hz, rx.noise_figure_db);
        result.snr_db = result.received_power_dbm - result.noise_floor_dbm;

        const T sensitivity_margin = result.received_power_dbm - rx.sensitivity_dbm;

        const T snr_margin = result.snr_db - required_snr_db;
        result.link_margin_db = std::min(snr_margin, sensitivity_margin);

        result.link_margin_db -= path.fade_margin_db;

        result.link_viable = result.link_margin_db > T(0);

        if (path.fade_margin_db > T(0)) {
            result.availability_percent = T(99.9) * (T(1.0) - std::exp(-path.fade_margin_db / T(10.0)));
        } else {
            result.availability_percent = T(50.0);
        }

        if (result.snr_db < required_snr_db) {
            result.limiting_factors.push_back("Insufficient SNR");
        }
        if (result.received_power_dbm < rx.sensitivity_dbm) {
            result.limiting_factors.push_back("Below receiver sensitivity");
        }
        if (path.path_loss_db > T(150.0)) {
            result.limiting_factors.push_back("High path loss");
        }
        if (path.fade_margin_db < T(10.0)) {
            result.limiting_factors.push_back("Insufficient fade margin");
        }

        return result;
    }

    static T calculate_required_transmit_power_dbm(const receiver_parameters& rx,
                                                 const path_parameters& path,
                                                 T bandwidth_hz,
                                                 T required_snr_db,
                                                 T tx_antenna_gain_db,
                                                 T tx_feedline_loss_db = T(0)) {
        const T noise_floor = calculate_noise_floor_dbm(bandwidth_hz, rx.noise_figure_db);
        const T required_rx_power = noise_floor + required_snr_db;

        const T required_eirp = required_rx_power - rx.antenna_gain_db + rx.feedline_loss_db +
                               rx.misc_losses_db + path.path_loss_db + path.atmospheric_loss_db +
                               path.rain_loss_db + path.polarization_loss_db + path.fade_margin_db;

        return required_eirp - tx_antenna_gain_db + tx_feedline_loss_db;
    }

    static T calculate_maximum_range_m(const transmitter_parameters& tx,
                                     const receiver_parameters& rx,
                                     T frequency_hz,
                                     T bandwidth_hz,
                                     T required_snr_db,
                                     T fade_margin_db = T(10.0)) {
        const T eirp = calculate_eirp_dbm(tx);
        const T noise_floor = calculate_noise_floor_dbm(bandwidth_hz, rx.noise_figure_db);
        const T required_rx_power = noise_floor + required_snr_db;

        const T available_budget = eirp + rx.antenna_gain_db - rx.feedline_loss_db - 
                                  rx.misc_losses_db - required_rx_power - fade_margin_db;

        const T max_fspl_db = available_budget;

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T max_range = wavelength * std::pow(T(10.0), max_fspl_db / T(20.0)) / (T(4.0) * constants::mathematical<T>::pi);

        return max_range;
    }

    static T calculate_system_noise_temperature_k(T antenna_noise_temp_k, T receiver_noise_figure_db) {
        core::check_positive(antenna_noise_temp_k, "Antenna noise temperature");

        const T receiver_noise_factor = std::pow(T(10.0), receiver_noise_figure_db / T(10.0));
        const T receiver_noise_temp = (receiver_noise_factor - T(1.0)) * constants::physical<T>::standard_temp;

        return antenna_noise_temp_k + receiver_noise_temp;
    }

    static T calculate_g_over_t_db_per_k(T antenna_gain_db, T system_noise_temperature_k) {
        core::check_positive(system_noise_temperature_k, "System noise temperature");

        const T noise_temp_db = T(10.0) * std::log10(system_noise_temperature_k);
        return antenna_gain_db - noise_temp_db;
    }

    static T calculate_data_rate_bps(T bandwidth_hz, T snr_db, const std::string& modulation_type = "QPSK") {
        core::check_positive(bandwidth_hz, "Bandwidth");

        const T snr_linear = std::pow(T(10.0), snr_db / T(10.0));

        T spectral_efficiency; 
        if (modulation_type == "BPSK" || modulation_type == "PSK31") {
            spectral_efficiency = T(1.0);
        } else if (modulation_type == "QPSK" || modulation_type == "PSK63") {
            spectral_efficiency = T(2.0);
        } else if (modulation_type == "16QAM") {
            spectral_efficiency = T(4.0);
        } else if (modulation_type == "64QAM") {
            spectral_efficiency = T(6.0);
        } else if (modulation_type == "FSK" || modulation_type == "RTTY") {
            spectral_efficiency = T(0.5); 
        } else if (modulation_type == "MFSK" || modulation_type == "OLIVIA") {
            spectral_efficiency = T(0.3); 
        } else if (modulation_type == "JT65" || modulation_type == "FT8") {
            spectral_efficiency = T(0.1); 
        } else {
            spectral_efficiency = T(2.0); 
        }

        const T shannon_limit = std::log2(T(1.0) + snr_linear);
        const T practical_efficiency = std::min(spectral_efficiency, shannon_limit * T(0.8)); 

        return bandwidth_hz * practical_efficiency;
    }

    static T calculate_ham_band_noise_floor_dbm(T frequency_hz, T bandwidth_hz = T(2500.0)) {

        const T thermal_noise = calculate_noise_floor_dbm(bandwidth_hz, T(0), constants::physical<T>::standard_temp;

        T atmospheric_noise_db = T(0);
        T man_made_noise_db = T(0);

        if (frequency_hz < 30e6) { 

            atmospheric_noise_db = T(50.0) - T(20.0) * std::log10(frequency_hz / 1e6); 
            man_made_noise_db = T(30.0) - T(10.0) * std::log10(frequency_hz / 1e6);
        } else if (frequency_hz < 300e6) { 
            atmospheric_noise_db = T(5.0);
            man_made_noise_db = T(20.0) - T(15.0) * std::log10(frequency_hz / 100e6);
        } else { 
            atmospheric_noise_db = T(0);
            man_made_noise_db = T(10.0) - T(10.0) * std::log10(frequency_hz / 1e9);
        }

        const T total_external_noise = T(10.0) * std::log10(
            std::pow(T(10.0), atmospheric_noise_db / T(10.0)) + 
            std::pow(T(10.0), man_made_noise_db / T(10.0))
        );

        return std::max(thermal_noise, thermal_noise + total_external_noise);
    }

    static T calculate_required_snr_for_ham_mode(const std::string& mode) {

        if (mode == "CW") {
            return T(3.0); 
        } else if (mode == "SSB" || mode == "AM") {
            return T(10.0); 
        } else if (mode == "FM") {
            return T(12.0); 
        } else if (mode == "PSK31") {
            return T(3.0); 
        } else if (mode == "RTTY") {
            return T(6.0); 
        } else if (mode == "FT8" || mode == "JT65") {
            return T(-24.0); 
        } else if (mode == "WSJT") {
            return T(-28.0); 
        } else if (mode == "PACKET" || mode == "APRS") {
            return T(12.0); 
        } else {
            return T(10.0); 
        }
    }

    static std::vector<link_analysis_result> analyze_frequency_sweep(const transmitter_parameters& tx,
                                                                   const receiver_parameters& rx,
                                                                   const vector_type& frequencies,
                                                                   T distance_m,
                                                                   T bandwidth_hz) {
        std::vector<link_analysis_result> results;
        results.reserve(frequencies.size());

        for (const auto& freq : frequencies) {
            path_parameters path;
            path.frequency_hz = freq;
            path.distance_m = distance_m;

            const T wavelength = constants::physical<T>::c / freq;
            path.path_loss_db = T(20.0) * std::log10(T(4.0) * constants::mathematical<T>::pi * distance_m / wavelength);

            path.atmospheric_loss_db = T(0.1) * distance_m / T(1000.0); 
            path.rain_loss_db = T(0);
            path.fade_margin_db = T(10.0);
            path.polarization_loss_db = T(0);

            results.push_back(analyze_link(tx, rx, path, bandwidth_hz));
        }

        return results;
    }

    static T calculate_diversity_improvement_db(int num_branches, T correlation_coefficient = T(0.3)) {
        core::check_positive(num_branches, "Number of diversity branches");
        core::check_range(correlation_coefficient, T(0), T(1), "Correlation coefficient");

        if (num_branches == 1) {
            return T(0);
        }

        const T ideal_gain = T(10.0) * std::log10(T(num_branches));
        const T correlation_penalty = T(10.0) * std::log10(T(1.0) + correlation_coefficient * (T(num_branches) - T(1.0)));

        return ideal_gain - correlation_penalty;
    }

    static T calculate_interference_margin_db(T desired_signal_dbm, T interference_power_dbm, T required_sir_db = T(10.0)) {
        const T sir_actual = desired_signal_dbm - interference_power_dbm;
        return sir_actual - required_sir_db;
    }

    static transmitter_parameters create_typical_transmitter(const std::string& application) {
        transmitter_parameters tx;

        if (application == "cellular_base_station") {
            tx.power_dbm = T(43.0);      
            tx.antenna_gain_db = T(17.0);
            tx.feedline_loss_db = T(2.0);
            tx.misc_losses_db = T(1.0);
        } else if (application == "wifi_ap") {
            tx.power_dbm = T(20.0);      
            tx.antenna_gain_db = T(5.0);
            tx.feedline_loss_db = T(0.5);
            tx.misc_losses_db = T(0.5);
        } else if (application == "satellite") {
            tx.power_dbm = T(50.0);      
            tx.antenna_gain_db = T(30.0);
            tx.feedline_loss_db = T(1.0);
            tx.misc_losses_db = T(2.0);
        } else if (application == "mobile_phone") {
            tx.power_dbm = T(23.0);      
            tx.antenna_gain_db = T(0.0);
            tx.feedline_loss_db = T(0);
            tx.misc_losses_db = T(1.0);
        } else if (application == "ham_handheld_1w") {
            tx.power_dbm = T(30.0);      
            tx.antenna_gain_db = T(0.0); 
            tx.feedline_loss_db = T(0.2);
            tx.misc_losses_db = T(0.5);
        } else if (application == "ham_mobile_100w") {
            tx.power_dbm = T(50.0);      
            tx.antenna_gain_db = T(3.0); 
            tx.feedline_loss_db = T(1.5);
            tx.misc_losses_db = T(1.0);
        } else if (application == "ham_base_1kw") {
            tx.power_dbm = T(60.0);      
            tx.antenna_gain_db = T(10.0); 
            tx.feedline_loss_db = T(2.0); 
            tx.misc_losses_db = T(1.5);
        } else if (application == "ham_legal_limit_2kw") {
            tx.power_dbm = T(63.0);      
            tx.antenna_gain_db = T(15.0); 
            tx.feedline_loss_db = T(2.5); 
            tx.misc_losses_db = T(2.0);
        } else if (application == "ham_qrp_5w") {
            tx.power_dbm = T(37.0);      
            tx.antenna_gain_db = T(2.0); 
            tx.feedline_loss_db = T(1.0);
            tx.misc_losses_db = T(0.5);
        } else if (application == "ham_microwave_10ghz") {
            tx.power_dbm = T(20.0);      
            tx.antenna_gain_db = T(30.0); 
            tx.feedline_loss_db = T(5.0); 
            tx.misc_losses_db = T(2.0);
        } else if (application == "ham_weak_signal_vhf") {
            tx.power_dbm = T(47.0);      
            tx.antenna_gain_db = T(12.0); 
            tx.feedline_loss_db = T(1.5);
            tx.misc_losses_db = T(1.0);
        } else {

            tx.power_dbm = T(10.0);
            tx.antenna_gain_db = T(2.0);
            tx.feedline_loss_db = T(1.0);
            tx.misc_losses_db = T(1.0);
        }

        tx.eirp_dbm = calculate_eirp_dbm(tx);
        return tx;
    }

    static receiver_parameters create_typical_receiver(const std::string& application) {
        receiver_parameters rx;

        if (application == "cellular_mobile") {
            rx.antenna_gain_db = T(0.0);
            rx.feedline_loss_db = T(0);
            rx.noise_figure_db = T(8.0);
            rx.sensitivity_dbm = T(-110.0);
            rx.misc_losses_db = T(1.0);
        } else if (application == "wifi_client") {
            rx.antenna_gain_db = T(2.0);
            rx.feedline_loss_db = T(0);
            rx.noise_figure_db = T(6.0);
            rx.sensitivity_dbm = T(-85.0);
            rx.misc_losses_db = T(0.5);
        } else if (application == "satellite_terminal") {
            rx.antenna_gain_db = T(35.0);
            rx.feedline_loss_db = T(1.5);
            rx.noise_figure_db = T(1.5);
            rx.sensitivity_dbm = T(-120.0);
            rx.misc_losses_db = T(1.0);
        } else if (application == "ham_handheld") {
            rx.antenna_gain_db = T(0.0);    
            rx.feedline_loss_db = T(0.2);
            rx.noise_figure_db = T(8.0);    
            rx.sensitivity_dbm = T(-120.0); 
            rx.misc_losses_db = T(0.5);
        } else if (application == "ham_mobile") {
            rx.antenna_gain_db = T(3.0);    
            rx.feedline_loss_db = T(1.5);
            rx.noise_figure_db = T(6.0);    
            rx.sensitivity_dbm = T(-125.0); 
            rx.misc_losses_db = T(1.0);
        } else if (application == "ham_base_hf") {
            rx.antenna_gain_db = T(5.0);    
            rx.feedline_loss_db = T(2.0);   
            rx.noise_figure_db = T(12.0);   
            rx.sensitivity_dbm = T(-130.0); 
            rx.misc_losses_db = T(1.5);
        } else if (application == "ham_base_vhf_uhf") {
            rx.antenna_gain_db = T(12.0);   
            rx.feedline_loss_db = T(1.5);
            rx.noise_figure_db = T(2.5);    
            rx.sensitivity_dbm = T(-135.0); 
            rx.misc_losses_db = T(1.0);
        } else if (application == "ham_microwave") {
            rx.antenna_gain_db = T(30.0);   
            rx.feedline_loss_db = T(5.0);   
            rx.noise_figure_db = T(4.0);    
            rx.sensitivity_dbm = T(-110.0); 
            rx.misc_losses_db = T(2.0);
        } else if (application == "ham_sdr_dongle") {
            rx.antenna_gain_db = T(0.0);    
            rx.feedline_loss_db = T(0.5);
            rx.noise_figure_db = T(12.0);   
            rx.sensitivity_dbm = T(-100.0); 
            rx.misc_losses_db = T(1.0);
        } else if (application == "ham_contest_station") {
            rx.antenna_gain_db = T(15.0);   
            rx.feedline_loss_db = T(2.5);   
            rx.noise_figure_db = T(1.5);    
            rx.sensitivity_dbm = T(-140.0); 
            rx.misc_losses_db = T(2.0);
        } else {

            rx.antenna_gain_db = T(2.0);
            rx.feedline_loss_db = T(1.0);
            rx.noise_figure_db = T(10.0);
            rx.sensitivity_dbm = T(-100.0);
            rx.misc_losses_db = T(1.0);
        }

        return rx;
    }

    static constexpr T thermal_noise_floor_dbm_per_hz() { return constants::rf_components<T>::thermal_noise_floor_dbm_per_hz; }
    static constexpr T typical_fade_margin_db() { return constants::rf_components<T>::typical_fade_margin_db; }
    static constexpr T typical_implementation_loss_db() { return constants::rf_components<T>::typical_implementation_loss_db; }
    static constexpr T standard_temperature_k() { return constants::physical<T>::standard_temp; }
};

using link_budget_f = link_budget<float>;
using link_budget_d = link_budget<double>;

} 
} 

#endif 