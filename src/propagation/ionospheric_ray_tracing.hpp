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

#ifndef RVL_PROPAGATION_IONOSPHERIC_RAY_TRACING_HPP
#define RVL_PROPAGATION_IONOSPHERIC_RAY_TRACING_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include "ionospheric_refractive_index.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>

namespace rvl {
namespace propagation {

template<typename T>
class ionospheric_ray_tracing {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;
    using vector3_type = std::array<T, 3>;

    struct ray_state {
        vector3_type position;      
        vector3_type direction;     
        T group_path;              
        T phase_path;              
        T time_delay;              
        T doppler_shift;           
        T absorption_loss;         
        bool is_valid;             
    };

    struct ionospheric_profile {
        vector_type height_km;     
        vector_type electron_density; 
        vector_type collision_frequency; 
        T earth_magnetic_field_t;  
        T magnetic_inclination_rad; 
        T magnetic_declination_rad; 
    };

    struct ray_hop {
        T launch_angle_rad;        
        T azimuth_angle_rad;       
        T hop_distance_km;         
        T maximum_height_km;       
        T group_delay_ms;          
        T phase_delay_ms;          
        T absorption_loss_db;      
        vector3_type reflection_point; 
        bool is_valid_hop;         
    };

    struct muf_analysis {
        T maximum_usable_frequency_hz; 
        T optimum_working_frequency_hz; 
        T lowest_usable_frequency_hz;  
        T critical_frequency_hz;       
        T penetration_frequency_hz;    
        vector_type frequency_range;   
        vector_type reliability;       
    };

    static ray_state integrate_ray_step(const ray_state& current_state,
                                      const ionospheric_profile& profile,
                                      T frequency_hz,
                                      T step_size_km = T(1.0)) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(step_size_km, "Step size");

        ray_state next_state = current_state;

        const auto& r = current_state.position;
        const auto& k = current_state.direction;

        const T height = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) - (constants::physical<T>::earth_radius / T(1000.0));

        auto n_squared = ionospheric_refractive_index<T>::calculate_n_squared(
            interpolate_electron_density(profile, height),
            frequency_hz,
            interpolate_collision_frequency(profile, height)
        );

        const T n = std::sqrt(std::max(n_squared.real(), T(0.001))); 

        vector3_type grad_n = calculate_refractive_index_gradient(profile, r, frequency_hz);

        next_state.position[0] = r[0] + step_size_km * k[0] / n;
        next_state.position[1] = r[1] + step_size_km * k[1] / n;
        next_state.position[2] = r[2] + step_size_km * k[2] / n;

        next_state.direction[0] = k[0] + step_size_km * grad_n[0];
        next_state.direction[1] = k[1] + step_size_km * grad_n[1];
        next_state.direction[2] = k[2] + step_size_km * grad_n[2];

        T k_magnitude = std::sqrt(next_state.direction[0] * next_state.direction[0] +
                                next_state.direction[1] * next_state.direction[1] +
                                next_state.direction[2] * next_state.direction[2]);

        if (k_magnitude > T(0)) {
            T norm_factor = n / k_magnitude;
            next_state.direction[0] *= norm_factor;
            next_state.direction[1] *= norm_factor;
            next_state.direction[2] *= norm_factor;
        }

        const T ds = step_size_km;
        const T group_velocity_factor = calculate_group_velocity_factor(n_squared);

        next_state.group_path += ds * group_velocity_factor;
        next_state.phase_path += ds * n;
        const T c_km_per_ms = constants::physical<T>::c / T(1e6); 
        next_state.time_delay += ds * group_velocity_factor / c_km_per_ms;

        const T new_height = calculate_height_from_position(next_state.position);
        next_state.is_valid = (new_height >= T(0)) && (new_height <= T(2000.0)) && (n > T(0.1));

        return next_state;
    }

    static std::vector<ray_state> trace_ray_path(const vector3_type& transmitter_pos,
                                                const vector3_type& receiver_pos,
                                                T launch_elevation_rad,
                                                T launch_azimuth_rad,
                                                T frequency_hz,
                                                const ionospheric_profile& profile,
                                                T step_size_km = T(1.0),
                                                size_t max_steps = 10000) {

        std::vector<ray_state> ray_path;
        ray_path.reserve(max_steps);

        ray_state current_state;
        current_state.position = transmitter_pos;
        current_state.group_path = T(0);
        current_state.phase_path = T(0);
        current_state.time_delay = T(0);
        current_state.doppler_shift = T(0);
        current_state.absorption_loss = T(0);
        current_state.is_valid = true;

        current_state.direction[0] = std::cos(launch_elevation_rad) * std::cos(launch_azimuth_rad);
        current_state.direction[1] = std::cos(launch_elevation_rad) * std::sin(launch_azimuth_rad);
        current_state.direction[2] = std::sin(launch_elevation_rad);

        for (size_t step = 0; step < max_steps && current_state.is_valid; ++step) {
            ray_path.push_back(current_state);

            current_state = integrate_ray_step(current_state, profile, frequency_hz, step_size_km);

            const T height = calculate_height_from_position(current_state.position);
            if (height <= T(0.1)) {

                current_state.is_valid = false;
            }

            const T distance_to_receiver = calculate_distance(current_state.position, receiver_pos);
            if (distance_to_receiver < step_size_km * T(2.0)) {

                break;
            }
        }

        return ray_path;
    }

    static muf_analysis calculate_muf(const vector3_type& transmitter_pos,
                                    const vector3_type& receiver_pos,
                                    const ionospheric_profile& profile,
                                    size_t num_elevation_angles = 180) {

        muf_analysis analysis;

        T max_electron_density = *std::max_element(profile.electron_density.begin(), 
                                                 profile.electron_density.end());

        analysis.critical_frequency_hz = calculate_plasma_frequency(max_electron_density);

        const T path_distance = calculate_great_circle_distance(transmitter_pos, receiver_pos);

        T max_frequency = T(0);
        T optimal_elevation = T(0);

        const T elevation_step = constants::mathematical<T>::half_pi / T(num_elevation_angles);

        for (size_t i = 1; i < num_elevation_angles; ++i) {
            const T elevation = T(i) * elevation_step;

            T test_frequency = analysis.critical_frequency_hz * T(3.0); 
            T frequency_step = test_frequency * T(0.1);

            for (int attempts = 0; attempts < 20; ++attempts) {
                auto ray_path = trace_ray_path(transmitter_pos, receiver_pos,
                                             elevation, T(0), 
                                             test_frequency, profile);

                if (!ray_path.empty() && ray_path.back().is_valid) {

                    if (test_frequency > max_frequency) {
                        max_frequency = test_frequency;
                        optimal_elevation = elevation;
                    }
                    test_frequency += frequency_step;
                } else {

                    test_frequency -= frequency_step;
                    frequency_step *= T(0.5);
                }

                if (frequency_step < analysis.critical_frequency_hz * T(0.001)) {
                    break; 
                }
            }
        }

        analysis.maximum_usable_frequency_hz = max_frequency;
        analysis.optimum_working_frequency_hz = max_frequency * T(0.85);
        analysis.lowest_usable_frequency_hz = max_frequency * T(0.1); 
        analysis.penetration_frequency_hz = max_frequency * T(1.05);

        return analysis;
    }

    static std::vector<ray_hop> calculate_multi_hop_path(const vector3_type& transmitter_pos,
                                                       const vector3_type& receiver_pos,
                                                       T frequency_hz,
                                                       const ionospheric_profile& profile,
                                                       size_t max_hops = 5) {

        std::vector<ray_hop> hops;

        const T total_distance = calculate_great_circle_distance(transmitter_pos, receiver_pos);
        const T typical_hop_distance = T(2000.0); 
        const size_t estimated_hops = static_cast<size_t>(std::ceil(total_distance / typical_hop_distance));

        if (estimated_hops > max_hops) {
            return hops; 
        }

        vector3_type current_pos = transmitter_pos;
        T remaining_distance = total_distance;

        for (size_t hop = 0; hop < estimated_hops && remaining_distance > T(100.0); ++hop) {
            ray_hop current_hop;

            current_hop.hop_distance_km = std::min(typical_hop_distance, remaining_distance);

            vector3_type target_pos = calculate_intermediate_position(current_pos, receiver_pos, 
                                                                   current_hop.hop_distance_km / remaining_distance);

            T optimal_elevation = find_optimal_elevation_angle(current_pos, target_pos, frequency_hz, profile);
            current_hop.launch_angle_rad = optimal_elevation;

            auto ray_path = trace_ray_path(current_pos, target_pos, optimal_elevation, T(0),
                                         frequency_hz, profile);

            if (!ray_path.empty() && ray_path.back().is_valid) {
                current_hop.is_valid_hop = true;
                current_hop.group_delay_ms = ray_path.back().time_delay;
                current_hop.absorption_loss_db = ray_path.back().absorption_loss;
                current_hop.maximum_height_km = find_maximum_height(ray_path);
                current_hop.reflection_point = find_reflection_point(ray_path);

                current_pos = target_pos;
                remaining_distance -= current_hop.hop_distance_km;
            } else {
                current_hop.is_valid_hop = false;
                break; 
            }

            hops.push_back(current_hop);
        }

        return hops;
    }

private:

    static T interpolate_electron_density(const ionospheric_profile& profile, T height_km) {

        if (height_km <= profile.height_km.front()) {
            return std::max(profile.electron_density.front(), T(1e6)); 
        }
        if (height_km >= profile.height_km.back()) {
            return std::max(profile.electron_density.back(), T(1e6));
        }

        for (size_t i = 0; i < profile.height_km.size() - 1; ++i) {
            if (height_km >= profile.height_km[i] && height_km <= profile.height_km[i + 1]) {
                T fraction = (height_km - profile.height_km[i]) / 
                           (profile.height_km[i + 1] - profile.height_km[i]);
                T interpolated = profile.electron_density[i] + 
                               fraction * (profile.electron_density[i + 1] - profile.electron_density[i]);
                return std::max(interpolated, T(1e6)); 
            }
        }
        return T(1e6); 
    }

    static T interpolate_collision_frequency(const ionospheric_profile& profile, T height_km) {

        if (height_km <= profile.height_km.front()) {
            return profile.collision_frequency.front();
        }
        if (height_km >= profile.height_km.back()) {
            return profile.collision_frequency.back();
        }

        for (size_t i = 0; i < profile.height_km.size() - 1; ++i) {
            if (height_km >= profile.height_km[i] && height_km <= profile.height_km[i + 1]) {
                T fraction = (height_km - profile.height_km[i]) / 
                           (profile.height_km[i + 1] - profile.height_km[i]);
                return profile.collision_frequency[i] + 
                       fraction * (profile.collision_frequency[i + 1] - profile.collision_frequency[i]);
            }
        }
        return T(1e6); 
    }

    static vector3_type calculate_refractive_index_gradient(const ionospheric_profile& profile,
                                                          const vector3_type& position,
                                                          T frequency_hz) {
        vector3_type gradient;
        const T delta = T(0.1); 

        T height_center = calculate_height_from_position(position);
        T n_center = std::sqrt(ionospheric_refractive_index<T>::calculate_n_squared(
            interpolate_electron_density(profile, height_center),
            frequency_hz,
            interpolate_collision_frequency(profile, height_center)
        ).real());

        for (int i = 0; i < 3; ++i) {
            vector3_type pos_plus = position;
            pos_plus[i] += delta;
            T height_plus = calculate_height_from_position(pos_plus);
            T n_plus = std::sqrt(ionospheric_refractive_index<T>::calculate_n_squared(
                interpolate_electron_density(profile, height_plus),
                frequency_hz,
                interpolate_collision_frequency(profile, height_plus)
            ).real());

            gradient[i] = (n_plus - n_center) / delta;
        }

        return gradient;
    }

    static T calculate_height_from_position(const vector3_type& position) {
        T radius = std::sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
        return radius - (constants::physical<T>::earth_radius / T(1000.0)); 
    }

    static T calculate_distance(const vector3_type& pos1, const vector3_type& pos2) {
        return std::sqrt(std::pow(pos1[0] - pos2[0], 2) +
                        std::pow(pos1[1] - pos2[1], 2) +
                        std::pow(pos1[2] - pos2[2], 2));
    }

    static T calculate_great_circle_distance(const vector3_type& pos1, const vector3_type& pos2) {

        T distance_3d = calculate_distance(pos1, pos2);
        T chord_to_arc_factor = constants::mathematical<T>::pi / T(2.0);
        return distance_3d * chord_to_arc_factor; 
    }

    static T calculate_plasma_frequency(T electron_density) {

        const T e_charge = T(1.602e-19);  
        const T m_electron = T(9.109e-31); 
        const T epsilon_0 = constants::physical<T>::epsilon_0;

        return T(1.0) / (T(2.0) * constants::mathematical<T>::pi) *
               std::sqrt(electron_density * e_charge * e_charge / (epsilon_0 * m_electron));
    }

    static T calculate_group_velocity_factor(const std::complex<T>& n_squared) {

        T n = std::sqrt(std::max(n_squared.real(), T(0.001)));

        return T(1.0) / n; 
    }

    static T find_optimal_elevation_angle(const vector3_type& start_pos,
                                         const vector3_type& end_pos,
                                         T frequency_hz,
                                         const ionospheric_profile& profile) {

        T best_elevation = constants::mathematical<T>::pi / T(6.0); 
        T min_distance = std::numeric_limits<T>::max();

        for (int angle_idx = 1; angle_idx < 90; ++angle_idx) {
            T test_elevation = T(angle_idx) * constants::mathematical<T>::pi / T(180.0);

            auto ray_path = trace_ray_path(start_pos, end_pos, test_elevation, T(0),
                                         frequency_hz, profile, T(5.0), 1000);

            if (!ray_path.empty()) {
                T final_distance = calculate_distance(ray_path.back().position, end_pos);
                if (final_distance < min_distance) {
                    min_distance = final_distance;
                    best_elevation = test_elevation;
                }
            }
        }

        return best_elevation;
    }

    static T find_maximum_height(const std::vector<ray_state>& ray_path) {
        T max_height = T(0);
        for (const auto& state : ray_path) {
            T height = calculate_height_from_position(state.position);
            max_height = std::max(max_height, height);
        }
        return max_height;
    }

    static vector3_type find_reflection_point(const std::vector<ray_state>& ray_path) {

        size_t max_height_idx = 0;
        T max_height = T(0);

        for (size_t i = 0; i < ray_path.size(); ++i) {
            T height = calculate_height_from_position(ray_path[i].position);
            if (height > max_height) {
                max_height = height;
                max_height_idx = i;
            }
        }

        return ray_path.empty() ? vector3_type{T(0), T(0), T(0)} : ray_path[max_height_idx].position;
    }

    static vector3_type calculate_intermediate_position(const vector3_type& start,
                                                      const vector3_type& end,
                                                      T fraction) {
        vector3_type result;
        result[0] = start[0] + fraction * (end[0] - start[0]);
        result[1] = start[1] + fraction * (end[1] - start[1]);
        result[2] = start[2] + fraction * (end[2] - start[2]);
        return result;
    }

public:

    static ionospheric_profile create_chapman_profile(T foF2_hz = T(10e6),
                                                     T hmF2_km = T(300.0),
                                                     T scale_height_km = T(50.0)) {
        ionospheric_profile profile;

        const size_t num_points = 100;
        profile.height_km.resize(num_points);
        profile.electron_density.resize(num_points);
        profile.collision_frequency.resize(num_points);

        for (size_t i = 0; i < num_points; ++i) {
            T height = T(60.0) + T(i) * T(940.0) / T(num_points - 1);
            profile.height_km[i] = height;

            T z = (height - hmF2_km) / scale_height_km;
            T ne_peak = std::pow(foF2_hz * T(2.0) * constants::mathematical<T>::pi / 
                               T(8.98), T(2.0)) * constants::physical<T>::epsilon_0 * 
                               T(9.109e-31) / std::pow(T(1.602e-19), T(2.0));

            profile.electron_density[i] = ne_peak * std::exp(T(1.0) - z - std::exp(-z));

            profile.collision_frequency[i] = T(1e7) * std::exp(-height / T(20.0));
        }

        profile.earth_magnetic_field_t = T(5e-5); 
        profile.magnetic_inclination_rad = T(1.2); 
        profile.magnetic_declination_rad = T(0.2); 

        return profile;
    }

    static constexpr T typical_foF2_mhz() { return T(10.0); }
    static constexpr T typical_hmF2_km() { return T(300.0); }
    static constexpr T typical_scale_height_km() { return T(50.0); }
    static constexpr T earth_radius_km() { return constants::physical<T>::earth_radius / T(1000.0); }
};

using ionospheric_ray_tracing_f = ionospheric_ray_tracing<float>;
using ionospheric_ray_tracing_d = ionospheric_ray_tracing<double>;

} 
} 

#endif 