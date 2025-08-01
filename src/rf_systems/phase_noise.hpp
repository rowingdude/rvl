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

#ifndef RVL_RF_SYSTEMS_PHASE_NOISE_HPP
#define RVL_RF_SYSTEMS_PHASE_NOISE_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace rf_systems {

template<typename T>
class phase_noise {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static T power_spectral_density_dbc_per_hz(T noise_power_w, T carrier_power_w) {
        core::check_positive(noise_power_w, "Noise power");
        core::check_positive(carrier_power_w, "Carrier power");

        const T power_ratio = noise_power_w / carrier_power_w;
        return T(10.0) * std::log10(power_ratio);
    }

    static T power_spectral_density_from_voltage(T noise_voltage_rms, T carrier_voltage_rms) {
        core::check_positive(noise_voltage_rms, "Noise voltage RMS");
        core::check_positive(carrier_voltage_rms, "Carrier voltage RMS");

        const T voltage_ratio_squared = std::pow(noise_voltage_rms / carrier_voltage_rms, T(2.0));
        return T(10.0) * std::log10(voltage_ratio_squared);
    }

    static T rms_phase_jitter_rad(T phase_noise_dbc_per_hz, T integration_bandwidth_hz) {
        core::check_positive(integration_bandwidth_hz, "Integration bandwidth");

        const T phase_noise_linear = std::pow(T(10.0), phase_noise_dbc_per_hz / T(10.0));
        const T integrated_phase_variance = T(2.0) * phase_noise_linear * integration_bandwidth_hz;

        return std::sqrt(integrated_phase_variance);
    }

    static T rms_phase_jitter_deg(T phase_noise_dbc_per_hz, T integration_bandwidth_hz) {
        const T jitter_rad = rms_phase_jitter_rad(phase_noise_dbc_per_hz, integration_bandwidth_hz);
        return jitter_rad * constants::mathematical<T>::rad_to_deg;
    }

    static T timing_jitter_seconds(T phase_noise_dbc_per_hz, T carrier_frequency_hz, T integration_bandwidth_hz) {
        core::check_positive(carrier_frequency_hz, "Carrier frequency");

        const T phase_jitter_rad = rms_phase_jitter_rad(phase_noise_dbc_per_hz, integration_bandwidth_hz);
        const T omega_c = T(2.0) * constants::mathematical<T>::pi * carrier_frequency_hz;

        return phase_jitter_rad / omega_c;
    }

    static T timing_jitter_femtoseconds(T phase_noise_dbc_per_hz, T carrier_frequency_hz, T integration_bandwidth_hz) {
        const T jitter_seconds = timing_jitter_seconds(phase_noise_dbc_per_hz, carrier_frequency_hz, integration_bandwidth_hz);
        return jitter_seconds * T(1e15);
    }

    static T phase_noise_from_q_factor(T q_factor, T carrier_frequency_hz, T offset_frequency_hz,
                                     T flicker_corner_hz = T(1000.0), T thermal_floor_dbc = T(-170.0)) {
        core::check_positive(q_factor, "Q factor");
        core::check_positive(carrier_frequency_hz, "Carrier frequency");
        core::check_positive(offset_frequency_hz, "Offset frequency");

        const T loaded_q = q_factor;
        const T omega_0 = T(2.0) * constants::mathematical<T>::pi * carrier_frequency_hz;
        const T omega_m = T(2.0) * constants::mathematical<T>::pi * offset_frequency_hz;

        T phase_noise_dbc;

        if (offset_frequency_hz < flicker_corner_hz) {
            const T flicker_term = T(10.0) * std::log10(flicker_corner_hz / offset_frequency_hz);
            const T resonator_term = T(20.0) * std::log10(omega_0 / (T(2.0) * loaded_q * omega_m));
            phase_noise_dbc = thermal_floor_dbc + flicker_term + resonator_term;
        } else {
            const T resonator_term = T(20.0) * std::log10(omega_0 / (T(2.0) * loaded_q * omega_m));
            phase_noise_dbc = thermal_floor_dbc + resonator_term;
        }

        return phase_noise_dbc;
    }

    static T additive_phase_noise_dbm_per_hz(T phase_noise_dbc_per_hz, T carrier_power_dbm) {
        return phase_noise_dbc_per_hz + carrier_power_dbm;
    }

    static T spur_to_phase_noise_ratio(T spur_level_dbc, T phase_noise_dbc_per_hz, T spur_bandwidth_hz = T(1.0)) {
        const T spur_power_density = spur_level_dbc - T(10.0) * std::log10(spur_bandwidth_hz);
        return spur_power_density - phase_noise_dbc_per_hz;
    }

    static T close_in_phase_noise_1_over_f3(T offset_frequency_hz, T reference_frequency_hz,
                                           T reference_phase_noise_dbc, T slope_db_per_decade = T(30.0)) {
        core::check_positive(offset_frequency_hz, "Offset frequency");
        core::check_positive(reference_frequency_hz, "Reference frequency");

        const T frequency_ratio = offset_frequency_hz / reference_frequency_hz;
        const T log_ratio = std::log10(frequency_ratio);

        return reference_phase_noise_dbc - slope_db_per_decade * log_ratio;
    }

    static T phase_noise_from_oscillator_specs(T carrier_frequency_hz, T offset_frequency_hz,
                                             T power_consumption_mw, T tuning_voltage_v = T(2.5)) {
        core::check_positive(carrier_frequency_hz, "Carrier frequency");
        core::check_positive(offset_frequency_hz, "Offset frequency");
        core::check_positive(power_consumption_mw, "Power consumption");

        const T base_phase_noise = T(-40.0) - T(20.0) * std::log10(carrier_frequency_hz / T(1e9));
        const T offset_term = T(-20.0) * std::log10(offset_frequency_hz / T(1e3));
        const T power_term = T(10.0) * std::log10(power_consumption_mw / T(10.0));

        return base_phase_noise + offset_term - power_term;
    }

    static void phase_noise_profile_batch(const vector_type& offset_frequencies_hz,
                                        vector_type& phase_noise_dbc_per_hz,
                                        T carrier_frequency_hz, T q_factor,
                                        T flicker_corner_hz = T(1000.0),
                                        T thermal_floor_dbc = T(-170.0)) {
        const size_t n = offset_frequencies_hz.size();
        if (phase_noise_dbc_per_hz.size() != n) {
            phase_noise_dbc_per_hz.resize(n);
        }

        for (size_t i = 0; i < n; ++i) {
            if (offset_frequencies_hz[i] <= T(0)) {
                throw core::invalid_argument_error("All offset frequencies must be positive");
            }

            phase_noise_dbc_per_hz[i] = phase_noise_from_q_factor(q_factor, carrier_frequency_hz,
                                                                 offset_frequencies_hz[i],
                                                                 flicker_corner_hz, thermal_floor_dbc);
        }
    }

    static void timing_jitter_batch(const vector_type& phase_noise_levels_dbc,
                                  const vector_type& carrier_frequencies_hz,
                                  vector_type& timing_jitter_fs,
                                  T integration_bandwidth_hz = T(1e6)) {
        const size_t n = phase_noise_levels_dbc.size();
        if (carrier_frequencies_hz.size() != n || timing_jitter_fs.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (carrier_frequencies_hz[i] <= T(0)) {
                throw core::invalid_argument_error("All carrier frequencies must be positive");
            }

            timing_jitter_fs[i] = timing_jitter_femtoseconds(phase_noise_levels_dbc[i],
                                                           carrier_frequencies_hz[i],
                                                           integration_bandwidth_hz);
        }
    }

    static constexpr T excellent_phase_noise_1khz_dbc() { return T(-120.0); }
    static constexpr T good_phase_noise_1khz_dbc() { return T(-100.0); }
    static constexpr T fair_phase_noise_1khz_dbc() { return T(-80.0); }
    static constexpr T poor_phase_noise_1khz_dbc() { return T(-60.0); }

    static constexpr T crystal_oscillator_q() { return T(1e6); }
    static constexpr T ceramic_resonator_q() { return T(1e3); }
    static constexpr T lc_oscillator_q() { return T(100.0); }
    static constexpr T rc_oscillator_q() { return T(10.0); }

    static constexpr T thermal_noise_floor_dbc() { return T(-174.0); }
    static constexpr T typical_flicker_corner_hz() { return T(1000.0); }

    static bool is_phase_noise_adequate(T phase_noise_dbc_per_hz, T requirement_dbc_per_hz) {
        return phase_noise_dbc_per_hz <= requirement_dbc_per_hz;
    }

    static T required_q_for_phase_noise(T target_phase_noise_dbc, T carrier_frequency_hz,
                                       T offset_frequency_hz, T thermal_floor_dbc = T(-170.0)) {
        const T omega_0 = T(2.0) * constants::mathematical<T>::pi * carrier_frequency_hz;
        const T omega_m = T(2.0) * constants::mathematical<T>::pi * offset_frequency_hz;

        const T numerator = target_phase_noise_dbc - thermal_floor_dbc;
        const T q_term = numerator / T(-20.0);
        const T q_factor = omega_0 / (T(2.0) * omega_m * std::pow(T(10.0), q_term));

        return q_factor;
    }
};

using phase_noise_f = phase_noise<float>;
using phase_noise_d = phase_noise<double>;

} 
} 

#endif 