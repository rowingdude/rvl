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

#ifndef RVL_PROPAGATION_RAY_TRAJECTORY_HPP
#define RVL_PROPAGATION_RAY_TRAJECTORY_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <array>
#include <functional>

namespace rvl {
namespace propagation {

template<typename T>
class ray_trajectory {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;
    using point3d_type = std::array<T, 3>;
    using vector3d_type = std::array<T, 3>;

    struct ray_state {
        point3d_type position;     
        vector3d_type wave_vector; 
        T path_length;             
        T refractive_index;        
    };

    static ray_state hamilton_step(const ray_state& current_state, 
                                 T step_size,
                                 const std::function<T(T, T, T)>& refractive_index_func,
                                 const std::function<vector3d_type(T, T, T)>& refractive_index_gradient_func) {
        const T K = std::sqrt(current_state.wave_vector[0] * current_state.wave_vector[0] +
                             current_state.wave_vector[1] * current_state.wave_vector[1] +
                             current_state.wave_vector[2] * current_state.wave_vector[2]);

        if (K <= T(0)) {
            throw core::invalid_argument_error("Wave vector magnitude must be positive");
        }

        const T x = current_state.position[0];
        const T y = current_state.position[1];
        const T z = current_state.position[2];
        const T kx = current_state.wave_vector[0];
        const T ky = current_state.wave_vector[1];
        const T kz = current_state.wave_vector[2];

        const T N = refractive_index_func(x, y, z);
        const vector3d_type grad_N = refractive_index_gradient_func(x, y, z);

        ray_state next_state;

        next_state.position[0] = x + step_size * (kx / K);
        next_state.position[1] = y + step_size * (ky / K);
        next_state.position[2] = z + step_size * (kz / K);

        next_state.wave_vector[0] = kx - step_size * (grad_N[0] * N / K);
        next_state.wave_vector[1] = ky - step_size * (grad_N[1] * N / K);
        next_state.wave_vector[2] = kz - step_size * (grad_N[2] * N / K);

        next_state.path_length = current_state.path_length + step_size;
        next_state.refractive_index = refractive_index_func(next_state.position[0], 
                                                          next_state.position[1], 
                                                          next_state.position[2]);

        return next_state;
    }

    static ray_state runge_kutta_4_step(const ray_state& current_state, 
                                       T step_size,
                                       const std::function<T(T, T, T)>& refractive_index_func,
                                       const std::function<vector3d_type(T, T, T)>& refractive_index_gradient_func) {

        auto derivative_func = [&](const ray_state& state) -> std::array<T, 6> {
            const T K = std::sqrt(state.wave_vector[0] * state.wave_vector[0] +
                                 state.wave_vector[1] * state.wave_vector[1] +
                                 state.wave_vector[2] * state.wave_vector[2]);

            const T N = refractive_index_func(state.position[0], state.position[1], state.position[2]);
            const vector3d_type grad_N = refractive_index_gradient_func(state.position[0], state.position[1], state.position[2]);

            return std::array<T, 6>{
                state.wave_vector[0] / K,                    
                state.wave_vector[1] / K,                    
                state.wave_vector[2] / K,                    
                -grad_N[0] * N / K,                          
                -grad_N[1] * N / K,                          
                -grad_N[2] * N / K                           
            };
        };

        const auto k1 = derivative_func(current_state);

        ray_state temp_state1 = current_state;
        temp_state1.position[0] += step_size * k1[0] / T(2.0);
        temp_state1.position[1] += step_size * k1[1] / T(2.0);
        temp_state1.position[2] += step_size * k1[2] / T(2.0);
        temp_state1.wave_vector[0] += step_size * k1[3] / T(2.0);
        temp_state1.wave_vector[1] += step_size * k1[4] / T(2.0);
        temp_state1.wave_vector[2] += step_size * k1[5] / T(2.0);

        const auto k2 = derivative_func(temp_state1);

        ray_state temp_state2 = current_state;
        temp_state2.position[0] += step_size * k2[0] / T(2.0);
        temp_state2.position[1] += step_size * k2[1] / T(2.0);
        temp_state2.position[2] += step_size * k2[2] / T(2.0);
        temp_state2.wave_vector[0] += step_size * k2[3] / T(2.0);
        temp_state2.wave_vector[1] += step_size * k2[4] / T(2.0);
        temp_state2.wave_vector[2] += step_size * k2[5] / T(2.0);

        const auto k3 = derivative_func(temp_state2);

        ray_state temp_state3 = current_state;
        temp_state3.position[0] += step_size * k3[0];
        temp_state3.position[1] += step_size * k3[1];
        temp_state3.position[2] += step_size * k3[2];
        temp_state3.wave_vector[0] += step_size * k3[3];
        temp_state3.wave_vector[1] += step_size * k3[4];
        temp_state3.wave_vector[2] += step_size * k3[5];

        const auto k4 = derivative_func(temp_state3);

        ray_state next_state;
        next_state.position[0] = current_state.position[0] + step_size * (k1[0] + T(2.0)*k2[0] + T(2.0)*k3[0] + k4[0]) / T(6.0);
        next_state.position[1] = current_state.position[1] + step_size * (k1[1] + T(2.0)*k2[1] + T(2.0)*k3[1] + k4[1]) / T(6.0);
        next_state.position[2] = current_state.position[2] + step_size * (k1[2] + T(2.0)*k2[2] + T(2.0)*k3[2] + k4[2]) / T(6.0);

        next_state.wave_vector[0] = current_state.wave_vector[0] + step_size * (k1[3] + T(2.0)*k2[3] + T(2.0)*k3[3] + k4[3]) / T(6.0);
        next_state.wave_vector[1] = current_state.wave_vector[1] + step_size * (k1[4] + T(2.0)*k2[4] + T(2.0)*k3[4] + k4[4]) / T(6.0);
        next_state.wave_vector[2] = current_state.wave_vector[2] + step_size * (k1[5] + T(2.0)*k2[5] + T(2.0)*k3[5] + k4[5]) / T(6.0);

        next_state.path_length = current_state.path_length + step_size;
        next_state.refractive_index = refractive_index_func(next_state.position[0], 
                                                          next_state.position[1], 
                                                          next_state.position[2]);

        return next_state;
    }

    static T calculate_ray_distance(const ray_state& state1, const ray_state& state2) {
        const T dx = state2.position[0] - state1.position[0];
        const T dy = state2.position[1] - state1.position[1];
        const T dz = state2.position[2] - state1.position[2];

        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    static T calculate_wave_vector_magnitude(const ray_state& state) {
        return std::sqrt(state.wave_vector[0] * state.wave_vector[0] +
                        state.wave_vector[1] * state.wave_vector[1] +
                        state.wave_vector[2] * state.wave_vector[2]);
    }

    static vector3d_type calculate_ray_direction(const ray_state& state) {
        const T K = calculate_wave_vector_magnitude(state);
        if (K <= T(0)) {
            throw core::invalid_argument_error("Wave vector magnitude must be positive");
        }

        return vector3d_type{
            state.wave_vector[0] / K,
            state.wave_vector[1] / K,
            state.wave_vector[2] / K
        };
    }

    static ray_state create_initial_state(const point3d_type& position, 
                                        const vector3d_type& initial_direction,
                                        T wave_number) {
        core::check_positive(wave_number, "Wave number");

        const T magnitude = std::sqrt(initial_direction[0] * initial_direction[0] +
                                     initial_direction[1] * initial_direction[1] +
                                     initial_direction[2] * initial_direction[2]);

        if (magnitude <= T(0)) {
            throw core::invalid_argument_error("Initial direction must have positive magnitude");
        }

        ray_state state;
        state.position = position;
        state.wave_vector[0] = wave_number * initial_direction[0] / magnitude;
        state.wave_vector[1] = wave_number * initial_direction[1] / magnitude;
        state.wave_vector[2] = wave_number * initial_direction[2] / magnitude;
        state.path_length = T(0);
        state.refractive_index = T(1);

        return state;
    }

    static bool is_ray_turning_point(const ray_state& state, T tolerance = T(1e-6)) {
        const T kz = state.wave_vector[2];
        return std::abs(kz) < tolerance;
    }

    static T calculate_elevation_angle_rad(const ray_state& state) {
        const vector3d_type direction = calculate_ray_direction(state);
        return std::asin(direction[2]);
    }

    static T calculate_elevation_angle_deg(const ray_state& state) {
        return calculate_elevation_angle_rad(state) * constants::mathematical<T>::rad_to_deg;
    }

    static constexpr T typical_step_size_km() { return T(1.0); }
    static constexpr T typical_wave_number() { return T(0.1); }
};

using ray_trajectory_f = ray_trajectory<float>;
using ray_trajectory_d = ray_trajectory<double>;

} 
} 

#endif 