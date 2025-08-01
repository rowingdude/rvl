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

#ifndef RVL_RF_SYSTEMS_INTERMODULATION_HPP
#define RVL_RF_SYSTEMS_INTERMODULATION_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace rf_systems {

template<typename T>
class intermodulation {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static T input_ip3_from_imd(T input_power_dbm, T imd_level_db) {
        return input_power_dbm + imd_level_db / T(2.0);
    }

    static T output_ip3_from_input_ip3(T input_ip3_dbm, T gain_db) {
        return input_ip3_dbm + gain_db;
    }

    static T input_ip3_from_output_ip3(T output_ip3_dbm, T gain_db) {
        return output_ip3_dbm - gain_db;
    }

    static T imd_level(T input_power_dbm, T input_ip3_dbm) {
        return T(2.0) * (input_ip3_dbm - input_power_dbm);
    }

    static T third_order_distortion_power(T input_power_dbm, T input_ip3_dbm) {
        const T delta = input_power_dbm - input_ip3_dbm;
        return input_power_dbm + T(2.0) * delta;
    }

    static T spurious_free_dynamic_range(T input_ip3_dbm, T noise_floor_dbm, T gain_db = T(0.0)) {
        core::check_finite(input_ip3_dbm, "Input IP3");
        core::check_finite(noise_floor_dbm, "Noise floor");

        const T output_ip3 = input_ip3_dbm + gain_db;
        const T output_noise_floor = noise_floor_dbm + gain_db;

        return (output_ip3 - output_noise_floor) * T(2.0) / T(3.0);
    }

    static T cascade_ip3_two_stages(T ip3_1_dbm, T ip3_2_dbm, T gain_1_db) {
        const T ip3_1_linear = std::pow(T(10.0), ip3_1_dbm / T(10.0));
        const T ip3_2_linear = std::pow(T(10.0), (ip3_2_dbm - gain_1_db) / T(10.0));
        const T gain_1_linear = std::pow(T(10.0), gain_1_db / T(10.0));

        const T ip3_total_linear = T(1.0) / (T(1.0) / ip3_1_linear + gain_1_linear / ip3_2_linear);

        return T(10.0) * std::log10(ip3_total_linear);
    }

    static T cascade_ip3_multiple_stages(const vector_type& ip3_values_dbm, const vector_type& gains_db) {
        if (ip3_values_dbm.size() != gains_db.size()) {
            throw core::dimension_mismatch_error("IP3 and gain vectors must have same size");
        }

        const size_t n = ip3_values_dbm.size();
        if (n == 0) {
            throw core::invalid_argument_error("Empty input vectors");
        }

        if (n == 1) {
            return ip3_values_dbm[0];
        }

        T cumulative_gain = T(0.0);
        T inverse_ip3_sum = T(0.0);

        for (size_t i = 0; i < n; ++i) {
            const T ip3_linear = std::pow(T(10.0), ip3_values_dbm[i] / T(10.0));
            const T gain_to_input_linear = std::pow(T(10.0), cumulative_gain / T(10.0));

            inverse_ip3_sum += gain_to_input_linear / ip3_linear;
            cumulative_gain += gains_db[i];
        }

        const T total_ip3_linear = T(1.0) / inverse_ip3_sum;
        return T(10.0) * std::log10(total_ip3_linear);
    }

    static T second_order_intercept_point(T input_power_dbm, T second_order_distortion_db) {
        return input_power_dbm + second_order_distortion_db;
    }

    static T conversion_compression_point(T input_ip3_dbm, T compression_point_offset_db = T(10.0)) {
        return input_ip3_dbm - compression_point_offset_db;
    }

    static T blocking_dynamic_range(T input_ip3_dbm, T sensitivity_dbm, T desired_signal_dbm) {
        const T interferer_level = input_ip3_dbm - (desired_signal_dbm - sensitivity_dbm) / T(2.0);
        return interferer_level - desired_signal_dbm;
    }

    static T two_tone_test_frequency_separation_mhz(T center_frequency_mhz, T bandwidth_mhz) {
        core::check_positive(center_frequency_mhz, "Center frequency");
        core::check_positive(bandwidth_mhz, "Bandwidth");

        return std::min(bandwidth_mhz / T(10.0), center_frequency_mhz / T(1000.0));
    }

    static void ip3_cascade_batch(const std::vector<vector_type>& stage_ip3_values,
                                const std::vector<vector_type>& stage_gains,
                                vector_type& cascaded_ip3_values) {
        if (stage_ip3_values.empty() || stage_gains.empty()) {
            throw core::invalid_argument_error("Empty input vectors");
        }

        if (stage_ip3_values.size() != stage_gains.size()) {
            throw core::dimension_mismatch_error("IP3 and gain vector arrays must have same size");
        }

        const size_t num_cascades = stage_ip3_values[0].size();

        for (const auto& vec : stage_ip3_values) {
            if (vec.size() != num_cascades) {
                throw core::dimension_mismatch_error("All IP3 vectors must have same size");
            }
        }

        for (const auto& vec : stage_gains) {
            if (vec.size() != num_cascades) {
                throw core::dimension_mismatch_error("All gain vectors must have same size");
            }
        }

        if (cascaded_ip3_values.size() != num_cascades) {
            cascaded_ip3_values.resize(num_cascades);
        }

        const size_t num_stages = stage_ip3_values.size();

        for (size_t cascade = 0; cascade < num_cascades; ++cascade) {
            vector_type ip3_vec(num_stages);
            vector_type gain_vec(num_stages);

            for (size_t stage = 0; stage < num_stages; ++stage) {
                ip3_vec[stage] = stage_ip3_values[stage][cascade];
                gain_vec[stage] = stage_gains[stage][cascade];
            }

            cascaded_ip3_values[cascade] = cascade_ip3_multiple_stages(ip3_vec, gain_vec);
        }
    }

    static void imd_level_batch(const vector_type& input_powers,
                              const vector_type& input_ip3_values,
                              vector_type& imd_levels) {
        const size_t n = input_powers.size();
        if (input_ip3_values.size() != n || imd_levels.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            imd_levels[i] = imd_level(input_powers[i], input_ip3_values[i]);
        }
    }

    static constexpr T typical_mixer_ip3_dbm() { return T(10.0); }
    static constexpr T typical_lna_ip3_dbm() { return T(-5.0); }
    static constexpr T typical_power_amp_ip3_dbm() { return T(30.0); }
    static constexpr T typical_vga_ip3_dbm() { return T(20.0); }

    static constexpr T excellent_ip3_dbm() { return T(40.0); }
    static constexpr T good_ip3_dbm() { return T(20.0); }
    static constexpr T fair_ip3_dbm() { return T(10.0); }
    static constexpr T poor_ip3_dbm() { return T(0.0); }

    static bool is_linearity_adequate(T input_ip3_dbm, T max_input_power_dbm, T margin_db = T(10.0)) {
        return input_ip3_dbm >= (max_input_power_dbm + margin_db);
    }

    static T required_ip3_for_sfdr(T sfdr_db, T noise_floor_dbm, T gain_db = T(0.0)) {
        const T output_noise_floor = noise_floor_dbm + gain_db;
        const T required_output_ip3 = output_noise_floor + T(1.5) * sfdr_db;
        return required_output_ip3 - gain_db;
    }
};

using intermodulation_f = intermodulation<float>;
using intermodulation_d = intermodulation<double>;

} 
} 

#endif 