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

#ifndef RVL_PROPAGATION_GROUP_PHASE_PATH_HPP
#define RVL_PROPAGATION_GROUP_PHASE_PATH_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include "ray_trajectory.hpp"
#include <cmath>
#include <functional>

namespace rvl {
namespace propagation {

template<typename T>
class group_phase_path {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;
    using ray_state_type = typename ray_trajectory<T>::ray_state;

    struct path_result {
        T group_path;          
        T phase_path;          
        T geometric_path;      
        T excess_path;         
        T group_delay_seconds; 
        T phase_delay_seconds; 
    };

    static T calculate_group_refractive_index(T refractive_index, T frequency, T df_dn) {
        core::check_positive(refractive_index, "Refractive index");
        core::check_positive(frequency, "Frequency");

        return refractive_index - frequency * df_dn;
    }

    static T calculate_group_refractive_index_numerical(T frequency, 
                                                      const std::function<T(T)>& refractive_index_func,
                                                      T delta_freq = T(1e3)) {
        core::check_positive(frequency, "Frequency");
        core::check_positive(delta_freq, "Delta frequency");

        const T n_plus = refractive_index_func(frequency + delta_freq);
        const T n_minus = refractive_index_func(frequency - delta_freq);
        const T n_current = refractive_index_func(frequency);

        const T dn_df = (n_plus - n_minus) / (T(2.0) * delta_freq);

        return n_current - frequency * dn_df;
    }

    static path_result integrate_path_simple(const std::vector<ray_state_type>& ray_path,
                                           T frequency,
                                           const std::function<T(T, T, T, T)>& group_index_func) {
        if (ray_path.size() < 2) {
            throw core::invalid_argument_error("Ray path must contain at least 2 points");
        }

        path_result result{};

        for (size_t i = 1; i < ray_path.size(); ++i) {
            const auto& prev_state = ray_path[i-1];
            const auto& curr_state = ray_path[i];

            const T ds = curr_state.path_length - prev_state.path_length;

            const T x_mid = (prev_state.position[0] + curr_state.position[0]) / T(2.0);
            const T y_mid = (prev_state.position[1] + curr_state.position[1]) / T(2.0);
            const T z_mid = (prev_state.position[2] + curr_state.position[2]) / T(2.0);

            const T n_mid = (prev_state.refractive_index + curr_state.refractive_index) / T(2.0);
            const T ng_mid = group_index_func(x_mid, y_mid, z_mid, frequency);

            result.phase_path += n_mid * ds;
            result.group_path += ds / ng_mid;
            result.geometric_path += ds;
        }

        const T c = constants::physical<T>::c;
        result.group_delay_seconds = result.group_path / c;
        result.phase_delay_seconds = result.phase_path / c;
        result.excess_path = result.group_path - result.geometric_path;

        return result;
    }

    static path_result integrate_path_trapezoidal(const std::vector<ray_state_type>& ray_path,
                                                T frequency,
                                                const std::function<T(T, T, T, T)>& group_index_func) {
        if (ray_path.size() < 2) {
            throw core::invalid_argument_error("Ray path must contain at least 2 points");
        }

        path_result result{};

        for (size_t i = 1; i < ray_path.size(); ++i) {
            const auto& prev_state = ray_path[i-1];
            const auto& curr_state = ray_path[i];

            const T ds = curr_state.path_length - prev_state.path_length;

            const T n_prev = prev_state.refractive_index;
            const T n_curr = curr_state.refractive_index;

            const T ng_prev = group_index_func(prev_state.position[0], prev_state.position[1], 
                                             prev_state.position[2], frequency);
            const T ng_curr = group_index_func(curr_state.position[0], curr_state.position[1], 
                                             curr_state.position[2], frequency);

            result.phase_path += ds * (n_prev + n_curr) / T(2.0);
            result.group_path += ds * T(2.0) / (T(1.0)/ng_prev + T(1.0)/ng_curr);
            result.geometric_path += ds;
        }

        const T c = constants::physical<T>::c;
        result.group_delay_seconds = result.group_path / c;
        result.phase_delay_seconds = result.phase_path / c;
        result.excess_path = result.group_path - result.geometric_path;

        return result;
    }

    static T calculate_doppler_shift(T transmit_frequency, T ray_velocity, T propagation_velocity = T(0)) {
        core::check_positive(transmit_frequency, "Transmit frequency");

        const T c = constants::physical<T>::c;
        const T effective_velocity = propagation_velocity != T(0) ? propagation_velocity : c;

        return transmit_frequency * ray_velocity / effective_velocity;
    }

    static T calculate_time_delay_difference(const path_result& result) {
        return result.group_delay_seconds - result.phase_delay_seconds;
    }

    static T calculate_multipath_delay_spread(const std::vector<path_result>& multiple_paths) {
        if (multiple_paths.empty()) {
            return T(0);
        }

        T min_delay = multiple_paths[0].group_delay_seconds;
        T max_delay = multiple_paths[0].group_delay_seconds;

        for (const auto& path : multiple_paths) {
            min_delay = std::min(min_delay, path.group_delay_seconds);
            max_delay = std::max(max_delay, path.group_delay_seconds);
        }

        return max_delay - min_delay;
    }

    static T calculate_path_loss_excess(const path_result& result, T free_space_path_loss_db) {
        const T excess_path_db = T(20.0) * std::log10(T(1.0) + result.excess_path / result.geometric_path);
        return free_space_path_loss_db + excess_path_db;
    }

    static void calculate_batch_delays(const vector_type& group_paths,
                                     const vector_type& phase_paths,
                                     vector_type& group_delays,
                                     vector_type& phase_delays,
                                     T propagation_velocity = T(0)) {
        const size_t n = group_paths.size();
        if (phase_paths.size() != n || group_delays.size() != n || phase_delays.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        const T c = propagation_velocity != T(0) ? propagation_velocity : constants::physical<T>::c;

        for (size_t i = 0; i < n; ++i) {
            if (group_paths[i] < T(0) || phase_paths[i] < T(0)) {
                throw core::invalid_argument_error("Path lengths must be non-negative");
            }

            group_delays[i] = group_paths[i] / c;
            phase_delays[i] = phase_paths[i] / c;
        }
    }

    static bool is_dispersive_medium(const path_result& result, T tolerance = T(1e-9)) {
        const T delay_difference = std::abs(result.group_delay_seconds - result.phase_delay_seconds);
        return delay_difference > tolerance;
    }

    static constexpr T typical_ionospheric_excess_path_km() { return T(0.5); }
    static constexpr T typical_group_delay_ms() { return T(2.0); }
    static constexpr T typical_frequency_hz() { return T(10e6); }
};

using group_phase_path_f = group_phase_path<float>;
using group_phase_path_d = group_phase_path<double>;

} 
} 

#endif 