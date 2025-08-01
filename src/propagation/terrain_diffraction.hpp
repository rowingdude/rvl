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

#ifndef RVL_PROPAGATION_TERRAIN_DIFFRACTION_HPP
#define RVL_PROPAGATION_TERRAIN_DIFFRACTION_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace rvl {
namespace propagation {

template<typename T>
class terrain_diffraction {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    struct terrain_point {
        T distance_m;
        T height_m;
    };

    struct diffraction_result {
        T diffraction_loss_db;
        T fresnel_clearance;
        T knife_edge_parameter;
        bool line_of_sight;
    };

    static T calculate_fresnel_parameter(T obstacle_height_m, T tx_distance_m, T rx_distance_m, 
                                       T frequency_hz, T tx_height_m = T(0), T rx_height_m = T(0)) {
        core::check_positive(tx_distance_m, "Transmitter distance");
        core::check_positive(rx_distance_m, "Receiver distance");
        core::check_positive(frequency_hz, "Frequency");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T total_distance = tx_distance_m + rx_distance_m;

        const T direct_line_height = tx_height_m + (rx_height_m - tx_height_m) * tx_distance_m / total_distance;
        const T clearance_height = obstacle_height_m - direct_line_height;

        const T fresnel_factor = std::sqrt(T(2.0) * total_distance / (wavelength * tx_distance_m * rx_distance_m));

        return clearance_height * fresnel_factor;
    }

    static T calculate_knife_edge_diffraction_loss_db(T fresnel_parameter) {
        const T v = fresnel_parameter;

        if (v <= T(-2.4)) {
            return T(0);
        } else if (v <= T(0)) {
            return T(20.0) * std::log10(T(0.5) - T(0.62) * v);
        } else if (v <= T(1.0)) {
            return T(20.0) * std::log10(T(0.5) * std::exp(-T(0.95) * v));
        } else if (v <= T(2.4)) {
            return T(20.0) * std::log10(T(0.4) - std::sqrt(T(0.1184) - std::pow(v - T(0.1), T(2))));
        } else {
            return T(20.0) * std::log10(T(0.225) / v);
        }
    }

    static T calculate_fresnel_zone_radius(T distance1_m, T distance2_m, T frequency_hz, int zone_number = 1) {
        core::check_positive(distance1_m, "Distance 1");
        core::check_positive(distance2_m, "Distance 2");
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(zone_number, "Zone number");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T total_distance = distance1_m + distance2_m;

        return std::sqrt(T(zone_number) * wavelength * distance1_m * distance2_m / total_distance);
    }

    static diffraction_result analyze_single_knife_edge(T obstacle_height_m, T tx_distance_m, T rx_distance_m,
                                                      T frequency_hz, T tx_height_m, T rx_height_m) {
        diffraction_result result;

        result.knife_edge_parameter = calculate_fresnel_parameter(obstacle_height_m, tx_distance_m, rx_distance_m,
                                                                frequency_hz, tx_height_m, rx_height_m);

        result.diffraction_loss_db = calculate_knife_edge_diffraction_loss_db(result.knife_edge_parameter);

        const T total_distance = tx_distance_m + rx_distance_m;
        const T direct_line_height = tx_height_m + (rx_height_m - tx_height_m) * tx_distance_m / total_distance;

        result.line_of_sight = obstacle_height_m < direct_line_height;
        result.fresnel_clearance = direct_line_height - obstacle_height_m;

        return result;
    }

    static T calculate_multiple_knife_edge_loss_db(const std::vector<terrain_point>& terrain_profile,
                                                 T tx_height_m, T rx_height_m, T frequency_hz) {
        if (terrain_profile.size() < 3) {
            return T(0);
        }

        T total_loss = T(0);

        for (size_t i = 1; i < terrain_profile.size() - 1; ++i) {
            const auto& prev_point = terrain_profile[i-1];
            const auto& current_point = terrain_profile[i];
            const auto& next_point = terrain_profile[i+1];

            const T d1 = current_point.distance_m - prev_point.distance_m;
            const T d2 = next_point.distance_m - current_point.distance_m;

            T effective_tx_height = (i == 1) ? tx_height_m + prev_point.height_m : prev_point.height_m;
            T effective_rx_height = (i == terrain_profile.size() - 2) ? rx_height_m + next_point.height_m : next_point.height_m;

            const T v = calculate_fresnel_parameter(current_point.height_m, d1, d2, frequency_hz, 
                                                  effective_tx_height, effective_rx_height);

            if (v > T(-0.8)) {
                const T edge_loss = calculate_knife_edge_diffraction_loss_db(v);
                total_loss += edge_loss;
            }
        }

        return total_loss;
    }

    static T calculate_bullington_diffraction_loss_db(const std::vector<terrain_point>& terrain_profile,
                                                    T tx_height_m, T rx_height_m, T frequency_hz) {
        if (terrain_profile.empty()) {
            return T(0);
        }

        const T total_distance = terrain_profile.back().distance_m - terrain_profile.front().distance_m;

        T max_v = T(-10.0);
        size_t critical_edge_index = 0;

        for (size_t i = 1; i < terrain_profile.size() - 1; ++i) {
            const T d1 = terrain_profile[i].distance_m - terrain_profile.front().distance_m;
            const T d2 = terrain_profile.back().distance_m - terrain_profile[i].distance_m;

            const T effective_tx_height = tx_height_m + terrain_profile.front().height_m;
            const T effective_rx_height = rx_height_m + terrain_profile.back().height_m;

            const T v = calculate_fresnel_parameter(terrain_profile[i].height_m, d1, d2, frequency_hz,
                                                  effective_tx_height, effective_rx_height);

            if (v > max_v) {
                max_v = v;
                critical_edge_index = i;
            }
        }

        if (max_v <= T(-0.8)) {
            return T(0);
        }

        return calculate_knife_edge_diffraction_loss_db(max_v);
    }

    static T calculate_epstein_peterson_loss_db(const std::vector<terrain_point>& terrain_profile,
                                              T tx_height_m, T rx_height_m, T frequency_hz) {
        if (terrain_profile.size() < 2) {
            return T(0);
        }

        std::vector<T> individual_losses;

        for (size_t i = 1; i < terrain_profile.size() - 1; ++i) {
            const T d1 = terrain_profile[i].distance_m - terrain_profile[0].distance_m;
            const T d2 = terrain_profile.back().distance_m - terrain_profile[i].distance_m;

            const T v = calculate_fresnel_parameter(terrain_profile[i].height_m, d1, d2, frequency_hz,
                                                  tx_height_m + terrain_profile[0].height_m,
                                                  rx_height_m + terrain_profile.back().height_m);

            if (v > T(-1.0)) {
                const T loss = calculate_knife_edge_diffraction_loss_db(v);
                individual_losses.push_back(loss);
            }
        }

        if (individual_losses.empty()) {
            return T(0);
        }

        T total_loss_linear = T(0);
        for (const T& loss_db : individual_losses) {
            total_loss_linear += std::pow(T(10.0), loss_db / T(10.0));
        }

        return T(10.0) * std::log10(total_loss_linear);
    }

    static T calculate_smooth_earth_diffraction_db(T distance_m, T frequency_hz, 
                                                 T tx_height_m, T rx_height_m, T earth_radius_factor = T(1.33)) {
        core::check_positive(distance_m, "Distance");
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(tx_height_m, "Transmitter height");
        core::check_positive(rx_height_m, "Receiver height");

        const T earth_radius = constants::physical<T>::earth_radius * earth_radius_factor;
        const T wavelength = constants::physical<T>::c / frequency_hz;

        const T horizon_tx = T(3.57) * std::sqrt(tx_height_m);
        const T horizon_rx = T(3.57) * std::sqrt(rx_height_m);
        const T radio_horizon = horizon_tx + horizon_rx;

        if (distance_m / T(1000.0) <= radio_horizon) {
            return T(0);
        }

        const T beta = T(2.0) * std::sqrt((distance_m / T(1000.0) - radio_horizon) / radio_horizon);
        const T G = T(0.05751) * beta - T(0.002158) * beta * beta * beta;

        const T F = T(11.0) + T(10.0) * std::log10(wavelength) - T(17.6) * std::log10(distance_m / T(1000.0));

        return std::max(F + G, T(0));
    }

    static void calculate_diffraction_losses_batch(const vector_type& obstacle_heights,
                                                 const vector_type& tx_distances,
                                                 const vector_type& rx_distances,
                                                 const vector_type& frequencies,
                                                 T tx_height_m, T rx_height_m,
                                                 vector_type& diffraction_losses) {
        const size_t n = obstacle_heights.size();
        if (tx_distances.size() != n || rx_distances.size() != n || 
            frequencies.size() != n || diffraction_losses.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            const T v = calculate_fresnel_parameter(obstacle_heights[i], tx_distances[i], rx_distances[i],
                                                  frequencies[i], tx_height_m, rx_height_m);

            diffraction_losses[i] = calculate_knife_edge_diffraction_loss_db(v);
        }
    }

    static std::vector<terrain_point> create_smooth_terrain_profile(T total_distance_m, T start_height_m, 
                                                                  T end_height_m, size_t num_points = 100) {
        std::vector<terrain_point> profile;
        profile.reserve(num_points);

        for (size_t i = 0; i < num_points; ++i) {
            terrain_point point;
            point.distance_m = total_distance_m * T(i) / T(num_points - 1);
            point.height_m = start_height_m + (end_height_m - start_height_m) * T(i) / T(num_points - 1);
            profile.push_back(point);
        }

        return profile;
    }

    static T calculate_clearance_requirement_m(T distance1_m, T distance2_m, T frequency_hz, T clearance_factor = T(0.6)) {
        const T fresnel_radius = calculate_fresnel_zone_radius(distance1_m, distance2_m, frequency_hz, 1);
        return fresnel_radius * clearance_factor;
    }

    static bool has_adequate_clearance(T obstacle_height_m, T tx_distance_m, T rx_distance_m,
                                     T frequency_hz, T tx_height_m, T rx_height_m, T clearance_factor = T(0.6)) {
        const T required_clearance = calculate_clearance_requirement_m(tx_distance_m, rx_distance_m, frequency_hz, clearance_factor);

        const T total_distance = tx_distance_m + rx_distance_m;
        const T direct_line_height = tx_height_m + (rx_height_m - tx_height_m) * tx_distance_m / total_distance;
        const T actual_clearance = direct_line_height - obstacle_height_m;

        return actual_clearance >= required_clearance;
    }

    static T calculate_effective_earth_radius(T atmospheric_refractivity_gradient = T(-40e-6)) {
        const T standard_earth_radius = constants::physical<T>::earth_radius;
        const T k_factor = T(1.0) / (T(1.0) + standard_earth_radius * atmospheric_refractivity_gradient);
        return standard_earth_radius * k_factor;
    }

    static constexpr T standard_atmosphere_k_factor() { return T(1.33); }
    static constexpr T minimum_clearance_factor() { return T(0.6); }
    static constexpr T maximum_fresnel_parameter() { return T(2.4); }
};

using terrain_diffraction_f = terrain_diffraction<float>;
using terrain_diffraction_d = terrain_diffraction<double>;

} 
} 

#endif 