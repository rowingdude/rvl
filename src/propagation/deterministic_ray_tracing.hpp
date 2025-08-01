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

#ifndef RVL_PROPAGATION_DETERMINISTIC_RAY_TRACING_HPP
#define RVL_PROPAGATION_DETERMINISTIC_RAY_TRACING_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <vector>
#include <array>
#include <algorithm>

namespace rvl {
namespace propagation {

template<typename T>
class deterministic_ray_tracing {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using point3d_type = std::array<T, 3>;
    using vector3d_type = std::array<T, 3>;

    struct material_properties {
        complex_type relative_permittivity;
        complex_type relative_permeability;
        T roughness_factor;
        std::string material_name;
    };

    struct ray_segment {
        point3d_type start_point;
        point3d_type end_point;
        T path_length;
        complex_type path_gain;
        T delay_s;
        int interaction_type; 
    };

    struct ray_path {
        std::vector<ray_segment> segments;
        T total_path_length;
        complex_type total_path_gain;
        T total_delay_s;
        int total_interactions;
        T departure_angle_elevation_rad;
        T departure_angle_azimuth_rad;
        T arrival_angle_elevation_rad;
        T arrival_angle_azimuth_rad;
    };

    struct geometric_obstacle {
        point3d_type center;
        vector3d_type dimensions; 
        material_properties material;
        int obstacle_type; 
    };

    static complex_type calculate_reflection_coefficient(T incident_angle_rad, T frequency_hz,
                                                       const material_properties& material,
                                                       bool horizontal_polarization = true) {
        core::check_range(incident_angle_rad, T(0), constants::mathematical<T>::pi / T(2.0), "Incident angle");
        core::check_positive(frequency_hz, "Frequency");

        const T cos_theta = std::cos(incident_angle_rad);
        const T sin_theta = std::sin(incident_angle_rad);

        const complex_type eps_r = material.relative_permittivity;
        const complex_type mu_r = material.relative_permeability;

        const complex_type sqrt_term = std::sqrt(eps_r * mu_r - complex_type(sin_theta * sin_theta, T(0)));

        complex_type reflection_coeff;
        if (horizontal_polarization) {

            const complex_type numerator = complex_type(cos_theta, T(0)) - sqrt_term / mu_r;
            const complex_type denominator = complex_type(cos_theta, T(0)) + sqrt_term / mu_r;
            reflection_coeff = numerator / denominator;
        } else {

            const complex_type numerator = eps_r * complex_type(cos_theta, T(0)) - sqrt_term;
            const complex_type denominator = eps_r * complex_type(cos_theta, T(0)) + sqrt_term;
            reflection_coeff = numerator / denominator;
        }

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T roughness_factor = std::exp(-T(2.0) * std::pow(T(2.0) * constants::mathematical<T>::pi * 
                                                              material.roughness_factor * cos_theta / wavelength, T(2)));

        return reflection_coeff * complex_type(roughness_factor, T(0));
    }

    static complex_type calculate_transmission_coefficient(T incident_angle_rad, T frequency_hz,
                                                         const material_properties& material,
                                                         bool horizontal_polarization = true) {
        core::check_range(incident_angle_rad, T(0), constants::mathematical<T>::pi / T(2.0), "Incident angle");
        core::check_positive(frequency_hz, "Frequency");

        const T cos_theta_i = std::cos(incident_angle_rad);
        const T sin_theta_i = std::sin(incident_angle_rad);

        const complex_type eps_r = material.relative_permittivity;
        const complex_type mu_r = material.relative_permeability;

        const complex_type n_rel = std::sqrt(eps_r * mu_r);
        const complex_type sin_theta_t_squared = (sin_theta_i * sin_theta_i) / (n_rel * n_rel);
        const complex_type cos_theta_t = std::sqrt(T(1.0) - sin_theta_t_squared);

        complex_type transmission_coeff;
        if (horizontal_polarization) {

            const complex_type numerator = T(2.0) * complex_type(cos_theta_i, T(0));
            const complex_type denominator = complex_type(cos_theta_i, T(0)) + cos_theta_t / mu_r;
            transmission_coeff = numerator / denominator;
        } else {

            const complex_type numerator = T(2.0) * complex_type(cos_theta_i, T(0));
            const complex_type denominator = eps_r * complex_type(cos_theta_i, T(0)) + cos_theta_t;
            transmission_coeff = numerator / denominator;
        }

        return transmission_coeff;
    }

    static T calculate_path_length(const point3d_type& point1, const point3d_type& point2) {
        const T dx = point2[0] - point1[0];
        const T dy = point2[1] - point1[1];
        const T dz = point2[2] - point1[2];

        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    static vector3d_type calculate_direction_vector(const point3d_type& from, const point3d_type& to) {
        const T dx = to[0] - from[0];
        const T dy = to[1] - from[1];
        const T dz = to[2] - from[2];

        const T length = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (length <= T(0)) {
            return {T(0), T(0), T(0)};
        }

        return {dx/length, dy/length, dz/length};
    }

    static T calculate_angle_between_vectors(const vector3d_type& v1, const vector3d_type& v2) {
        const T dot_product = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
        const T v1_mag = std::sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
        const T v2_mag = std::sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);

        if (v1_mag <= T(0) || v2_mag <= T(0)) {
            return T(0);
        }

        const T cos_angle = std::max(T(-1.0), std::min(T(1.0), dot_product / (v1_mag * v2_mag)));
        return std::acos(cos_angle);
    }

    static ray_path trace_direct_path(const point3d_type& transmitter, const point3d_type& receiver,
                                    T frequency_hz) {
        ray_path path;

        ray_segment direct_segment;
        direct_segment.start_point = transmitter;
        direct_segment.end_point = receiver;
        direct_segment.path_length = calculate_path_length(transmitter, receiver);
        direct_segment.interaction_type = 0; 

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const complex_type free_space_gain = complex_type(wavelength / (T(4.0) * constants::mathematical<T>::pi * direct_segment.path_length), T(0));
        direct_segment.path_gain = free_space_gain;
        direct_segment.delay_s = direct_segment.path_length / constants::physical<T>::c;

        path.segments.push_back(direct_segment);
        path.total_path_length = direct_segment.path_length;
        path.total_path_gain = direct_segment.path_gain;
        path.total_delay_s = direct_segment.delay_s;
        path.total_interactions = 0;

        const auto direction = calculate_direction_vector(transmitter, receiver);
        path.departure_angle_elevation_rad = std::asin(direction[2]);
        path.departure_angle_azimuth_rad = std::atan2(direction[1], direction[0]);
        path.arrival_angle_elevation_rad = -path.departure_angle_elevation_rad;
        path.arrival_angle_azimuth_rad = path.departure_angle_azimuth_rad + constants::mathematical<T>::pi;

        return path;
    }

    static ray_path trace_single_reflection_path(const point3d_type& transmitter, const point3d_type& receiver,
                                               const geometric_obstacle& reflector, T frequency_hz) {
        ray_path path;

        point3d_type reflection_point = reflector.center;

        ray_segment segment1;
        segment1.start_point = transmitter;
        segment1.end_point = reflection_point;
        segment1.path_length = calculate_path_length(transmitter, reflection_point);
        segment1.interaction_type = 0; 

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const complex_type gain1 = complex_type(wavelength / (T(4.0) * constants::mathematical<T>::pi * segment1.path_length), T(0));
        segment1.path_gain = gain1;
        segment1.delay_s = segment1.path_length / constants::physical<T>::c;

        ray_segment segment2;
        segment2.start_point = reflection_point;
        segment2.end_point = receiver;
        segment2.path_length = calculate_path_length(reflection_point, receiver);
        segment2.interaction_type = 1; 

        const auto incident_dir = calculate_direction_vector(transmitter, reflection_point);
        const auto surface_normal = vector3d_type{T(0), T(0), T(1)}; 
        const T incident_angle_full = calculate_angle_between_vectors(incident_dir, surface_normal);

        const T grazing_angle = std::min(incident_angle_full, constants::mathematical<T>::pi - incident_angle_full);
        const T grazing_angle_clamped = std::min(grazing_angle, constants::mathematical<T>::pi / T(2.0));

        const complex_type reflection_coeff = calculate_reflection_coefficient(
            grazing_angle_clamped, frequency_hz, reflector.material);

        const complex_type gain2 = complex_type(wavelength / (T(4.0) * constants::mathematical<T>::pi * segment2.path_length), T(0));
        segment2.path_gain = gain2 * reflection_coeff;
        segment2.delay_s = segment2.path_length / constants::physical<T>::c;

        path.segments.push_back(segment1);
        path.segments.push_back(segment2);
        path.total_path_length = segment1.path_length + segment2.path_length;
        path.total_path_gain = segment1.path_gain * segment2.path_gain;
        path.total_delay_s = segment1.delay_s + segment2.delay_s;
        path.total_interactions = 1;

        return path;
    }

    static T calculate_received_power_dbm(T transmit_power_dbm, T tx_antenna_gain_db, T rx_antenna_gain_db,
                                        const std::vector<ray_path>& paths) {
        if (paths.empty()) {
            return T(-300.0); 
        }

        complex_type total_field(T(0), T(0));

        for (const auto& path : paths) {

            const T phase_delay = T(2.0) * constants::mathematical<T>::pi * path.total_delay_s;
            const complex_type phase_term = std::exp(complex_type(T(0), -phase_delay));

            total_field += path.total_path_gain * phase_term;
        }

        const T field_magnitude_squared = std::abs(total_field) * std::abs(total_field);

        if (field_magnitude_squared <= T(0)) {
            return T(-300.0);
        }

        const T path_gain_db = T(10.0) * std::log10(field_magnitude_squared);

        return transmit_power_dbm + tx_antenna_gain_db + rx_antenna_gain_db + path_gain_db;
    }

    static std::vector<ray_path> calculate_multipath_profile(const point3d_type& transmitter,
                                                           const point3d_type& receiver,
                                                           const std::vector<geometric_obstacle>& obstacles,
                                                           T frequency_hz,
                                                           int max_reflections = 3) {
        std::vector<ray_path> paths;

        paths.push_back(trace_direct_path(transmitter, receiver, frequency_hz));

        for (const auto& obstacle : obstacles) {
            auto reflected_path = trace_single_reflection_path(transmitter, receiver, obstacle, frequency_hz);

            if (std::abs(reflected_path.total_path_gain) > T(1e-6)) {
                paths.push_back(reflected_path);
            }
        }

        std::sort(paths.begin(), paths.end(), 
                 [](const ray_path& a, const ray_path& b) {
                     return a.total_delay_s < b.total_delay_s;
                 });

        return paths;
    }

    static T calculate_rms_delay_spread_s(const std::vector<ray_path>& paths) {
        if (paths.empty()) {
            return T(0);
        }

        T total_power = T(0);
        T weighted_delay_sum = T(0);

        for (const auto& path : paths) {
            const T power = std::abs(path.total_path_gain) * std::abs(path.total_path_gain);
            total_power += power;
            weighted_delay_sum += power * path.total_delay_s;
        }

        if (total_power <= T(0)) {
            return T(0);
        }

        const T mean_delay = weighted_delay_sum / total_power;

        T weighted_delay_variance = T(0);
        for (const auto& path : paths) {
            const T power = std::abs(path.total_path_gain) * std::abs(path.total_path_gain);
            const T delay_diff = path.total_delay_s - mean_delay;
            weighted_delay_variance += power * delay_diff * delay_diff;
        }

        return std::sqrt(weighted_delay_variance / total_power);
    }

    static material_properties create_material_properties(const std::string& material_name) {
        material_properties props;
        props.material_name = material_name;

        if (material_name == "concrete") {
            props.relative_permittivity = complex_type(T(5.5), T(-0.5));
            props.relative_permeability = complex_type(T(1.0), T(0));
            props.roughness_factor = T(0.02);
        } else if (material_name == "brick") {
            props.relative_permittivity = complex_type(T(4.0), T(-0.3));
            props.relative_permeability = complex_type(T(1.0), T(0));
            props.roughness_factor = T(0.05);
        } else if (material_name == "glass") {
            props.relative_permittivity = complex_type(T(6.0), T(-0.1));
            props.relative_permeability = complex_type(T(1.0), T(0));
            props.roughness_factor = T(0.001);
        } else if (material_name == "metal") {
            props.relative_permittivity = complex_type(T(1.0), T(-1000.0));
            props.relative_permeability = complex_type(T(1.0), T(0));
            props.roughness_factor = T(0.01);
        } else if (material_name == "water") {
            props.relative_permittivity = complex_type(T(81.0), T(-10.0));
            props.relative_permeability = complex_type(T(1.0), T(0));
            props.roughness_factor = T(0.1);
        } else {

            props.relative_permittivity = complex_type(T(3.0), T(-0.1));
            props.relative_permeability = complex_type(T(1.0), T(0));
            props.roughness_factor = T(0.01);
        }

        return props;
    }

    static bool is_line_of_sight(const point3d_type& transmitter, const point3d_type& receiver,
                                const std::vector<geometric_obstacle>& obstacles) {

        const auto direction = calculate_direction_vector(transmitter, receiver);
        const T path_length = calculate_path_length(transmitter, receiver);

        for (const auto& obstacle : obstacles) {

            const T obstacle_distance = calculate_path_length(transmitter, obstacle.center);

            if (obstacle_distance < path_length && obstacle_distance > T(1.0)) {

                return false;
            }
        }

        return true;
    }

    static constexpr T typical_concrete_permittivity() { return T(5.5); }
    static constexpr T typical_glass_permittivity() { return T(6.0); }
    static constexpr T typical_metal_conductivity() { return T(10000000.0); }
    static constexpr int maximum_ray_bounces() { return 5; }
};

using deterministic_ray_tracing_f = deterministic_ray_tracing<float>;
using deterministic_ray_tracing_d = deterministic_ray_tracing<double>;

} 
} 

#endif 