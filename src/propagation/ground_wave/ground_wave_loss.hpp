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

#ifndef RVL_PROPAGATION_GROUND_WAVE_GROUND_WAVE_LOSS_HPP
#define RVL_PROPAGATION_GROUND_WAVE_GROUND_WAVE_LOSS_HPP

#include "../../core/constants.hpp"
#include "../../core/error.hpp"
#include "../../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace propagation {
namespace ground_wave {

template<typename T>
class ground_wave_loss {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    enum class ground_type {
        POOR,       
        AVERAGE,    
        GOOD,       
        SEAWATER    
    };

    static T terrain_constant(ground_type type, T frequency_hz) {
        core::check_positive(frequency_hz, "Frequency");

        const T freq_mhz = frequency_hz / T(1e6);

        switch (type) {
            case ground_type::POOR:
                return T(40.0) + T(5.0) * std::log10(freq_mhz);
            case ground_type::AVERAGE:
                return T(32.0) + T(3.0) * std::log10(freq_mhz);
            case ground_type::GOOD:
                return T(25.0) + T(2.0) * std::log10(freq_mhz);
            case ground_type::SEAWATER:
                return T(5.0) + T(1.0) * std::log10(freq_mhz);
            default:
                return T(32.0) + T(3.0) * std::log10(freq_mhz);
        }
    }

    static T calculate_db(T distance_m, T frequency_hz, ground_type type = ground_type::AVERAGE) {
        core::check_positive(distance_m, "Distance");
        core::check_positive(frequency_hz, "Frequency");

        const T distance_km = distance_m / T(1000.0);
        const T freq_mhz = frequency_hz / T(1e6);
        const T K = terrain_constant(type, frequency_hz);

        return T(20.0) * std::log10(distance_km) + T(20.0) * std::log10(freq_mhz) + K;
    }

    static T calculate_with_custom_constant(T distance_m, T frequency_hz, T terrain_constant_db) {
        core::check_positive(distance_m, "Distance");
        core::check_positive(frequency_hz, "Frequency");

        const T distance_km = distance_m / T(1000.0);
        const T freq_mhz = frequency_hz / T(1e6);

        return T(20.0) * std::log10(distance_km) + T(20.0) * std::log10(freq_mhz) + terrain_constant_db;
    }

    static T calculate_with_conductivity(T distance_m, T frequency_hz, 
                                       T conductivity_s_per_m, T relative_permittivity) {
        core::check_positive(distance_m, "Distance");
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(conductivity_s_per_m, "Conductivity");
        core::check_positive(relative_permittivity, "Relative permittivity");

        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;
        const T epsilon_0 = constants::physical<T>::epsilon_0;

        const T sigma_over_omega_epsilon = conductivity_s_per_m / (omega * epsilon_0 * relative_permittivity);

        T correction_factor;
        if (sigma_over_omega_epsilon > T(10.0)) {
            correction_factor = T(-5.0);
        } else if (sigma_over_omega_epsilon > T(1.0)) {
            correction_factor = T(-2.0) * std::log10(sigma_over_omega_epsilon);
        } else {
            correction_factor = T(5.0) * std::log10(sigma_over_omega_epsilon);
        }

        const T basic_loss = calculate_db(distance_m, frequency_hz, ground_type::AVERAGE);
        return basic_loss + correction_factor;
    }

    static T diffraction_loss_knife_edge(T distance1_m, T distance2_m, T height_m, T frequency_hz) {
        core::check_positive(distance1_m, "Distance 1");
        core::check_positive(distance2_m, "Distance 2");
        core::check_positive(frequency_hz, "Frequency");

        const T total_distance = distance1_m + distance2_m;
        const T wavelength = constants::physical<T>::c / frequency_hz;

        const T fresnel_parameter = height_m * std::sqrt(T(2.0) * total_distance / 
                                                        (wavelength * distance1_m * distance2_m));

        T diffraction_loss;
        if (fresnel_parameter <= T(-2.4)) {
            diffraction_loss = T(0.0);
        } else if (fresnel_parameter <= T(0.0)) {
            diffraction_loss = T(20.0) * fresnel_parameter + T(0.1) * fresnel_parameter * fresnel_parameter;
        } else if (fresnel_parameter <= T(2.4)) {
            diffraction_loss = T(20.0) * std::log10(T(0.5) + T(0.62) * fresnel_parameter);
        } else {
            diffraction_loss = T(20.0) * std::log10(T(0.225) / fresnel_parameter);
        }

        return std::max(T(0.0), diffraction_loss);
    }

    static void calculate_batch(const vector_type& distances,
                              const vector_type& frequencies,
                              vector_type& losses,
                              ground_type type = ground_type::AVERAGE) {
        const size_t n = distances.size();
        if (frequencies.size() != n || losses.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (distances[i] <= T(0) || frequencies[i] <= T(0)) {
                throw core::invalid_argument_error("All distances and frequencies must be positive");
            }

            losses[i] = calculate_db(distances[i], frequencies[i], type);
        }
    }

    static void calculate_with_constants_batch(const vector_type& distances,
                                             const vector_type& frequencies,
                                             const vector_type& terrain_constants,
                                             vector_type& losses) {
        const size_t n = distances.size();
        if (frequencies.size() != n || terrain_constants.size() != n || losses.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (distances[i] <= T(0) || frequencies[i] <= T(0)) {
                throw core::invalid_argument_error("All distances and frequencies must be positive");
            }

            losses[i] = calculate_with_custom_constant(distances[i], frequencies[i], terrain_constants[i]);
        }
    }

    static constexpr T typical_vhf_terrain_constant() {
        return T(35.0);
    }

    static constexpr T typical_uhf_terrain_constant() {
        return T(40.0);
    }
};

using ground_wave_loss_f = ground_wave_loss<float>;
using ground_wave_loss_d = ground_wave_loss<double>;

} 
} 
} 

#endif 