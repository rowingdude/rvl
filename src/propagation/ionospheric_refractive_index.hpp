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

#ifndef RVL_PROPAGATION_IONOSPHERIC_REFRACTIVE_INDEX_HPP
#define RVL_PROPAGATION_IONOSPHERIC_REFRACTIVE_INDEX_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>

namespace rvl {
namespace propagation {

template<typename T>
class ionospheric_refractive_index {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    static complex_type calculate_n_squared(T electron_density, T wave_frequency, T collision_frequency = T(0)) {
        core::check_positive(electron_density, "Electron density");
        core::check_positive(wave_frequency, "Wave frequency");
        core::check_non_negative(collision_frequency, "Collision frequency");

        const T e = constants::physical<T>::q_e;
        const T m_e = constants::physical<T>::m_e;
        const T eps_0 = constants::physical<T>::epsilon_0;

        const T omega = T(2.0) * constants::mathematical<T>::pi * wave_frequency;
        const T omega_p_squared = (electron_density * e * e) / (m_e * eps_0);

        const T X = omega_p_squared / (omega * omega);
        const T Z = collision_frequency / omega;

        const complex_type denominator(T(1.0) - X / T(2.0), -Z);
        const complex_type n_squared = T(1.0) - X / denominator;

        return n_squared;
    }

    static T calculate_refractive_index_magnitude(T electron_density, T wave_frequency, T collision_frequency = T(0)) {
        const complex_type n_squared = calculate_n_squared(electron_density, wave_frequency, collision_frequency);
        return std::sqrt(std::abs(n_squared));
    }

    static complex_type calculate_plasma_frequency(T electron_density) {
        core::check_positive(electron_density, "Electron density");

        const T e = constants::physical<T>::q_e;
        const T m_e = constants::physical<T>::m_e;
        const T eps_0 = constants::physical<T>::epsilon_0;

        const T omega_p_squared = (electron_density * e * e) / (m_e * eps_0);
        const T omega_p = std::sqrt(omega_p_squared);

        return complex_type(omega_p, T(0));
    }

    static T calculate_critical_frequency(T electron_density) {
        const complex_type omega_p = calculate_plasma_frequency(electron_density);
        return omega_p.real() / (T(2.0) * constants::mathematical<T>::pi);
    }

    static T calculate_muf_factor(T electron_density, T wave_frequency, T elevation_angle_rad) {
        core::check_positive(electron_density, "Electron density");
        core::check_positive(wave_frequency, "Wave frequency");
        core::check_range(elevation_angle_rad, T(0), constants::mathematical<T>::pi / T(2.0), "Elevation angle");

        const T critical_freq = calculate_critical_frequency(electron_density);
        const T secant_angle = T(1.0) / std::sin(elevation_angle_rad);

        return critical_freq * secant_angle;
    }

    static T calculate_group_refractive_index(T electron_density, T wave_frequency, T collision_frequency = T(0)) {
        const T h = wave_frequency * T(1e-6);
        const T n1 = calculate_refractive_index_magnitude(electron_density, wave_frequency, collision_frequency);
        const T n2 = calculate_refractive_index_magnitude(electron_density, wave_frequency + h, collision_frequency);

        const T dn_df = (n2 - n1) / h;
        return n1 - wave_frequency * dn_df;
    }

    static complex_type calculate_absorption_coefficient(T electron_density, T wave_frequency, T collision_frequency) {
        core::check_positive(electron_density, "Electron density");
        core::check_positive(wave_frequency, "Wave frequency");
        core::check_positive(collision_frequency, "Collision frequency");

        const complex_type n_squared = calculate_n_squared(electron_density, wave_frequency, collision_frequency);
        const complex_type n = std::sqrt(n_squared);

        const T k0 = T(2.0) * constants::mathematical<T>::pi * wave_frequency / constants::physical<T>::c;
        return k0 * n.imag();
    }

    static void calculate_n_squared_batch(const vector_type& electron_densities,
                                        const vector_type& wave_frequencies,
                                        const vector_type& collision_frequencies,
                                        complex_vector_type& n_squared_results) {
        const size_t n = electron_densities.size();
        if (wave_frequencies.size() != n || collision_frequencies.size() != n || n_squared_results.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (electron_densities[i] <= T(0) || wave_frequencies[i] <= T(0) || collision_frequencies[i] < T(0)) {
                throw core::invalid_argument_error("Invalid parameters at index");
            }

            n_squared_results[i] = calculate_n_squared(electron_densities[i], wave_frequencies[i], collision_frequencies[i]);
        }
    }

    static void calculate_critical_frequency_batch(const vector_type& electron_densities,
                                                 vector_type& critical_frequencies) {
        const size_t n = electron_densities.size();
        if (critical_frequencies.size() != n) {
            throw core::dimension_mismatch_error("Vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (electron_densities[i] <= T(0)) {
                throw core::invalid_argument_error("Electron density must be positive");
            }

            critical_frequencies[i] = calculate_critical_frequency(electron_densities[i]);
        }
    }

    static bool is_wave_reflected(T electron_density, T wave_frequency, T elevation_angle_rad) {
        const T muf = calculate_muf_factor(electron_density, wave_frequency, elevation_angle_rad);
        return wave_frequency < muf;
    }

    static constexpr T typical_electron_density_f2() { return T(1e12); }
    static constexpr T typical_critical_frequency_mhz() { return T(10.0); }
    static constexpr T typical_collision_frequency() { return T(1e3); }
};

using ionospheric_refractive_index_f = ionospheric_refractive_index<float>;
using ionospheric_refractive_index_d = ionospheric_refractive_index<double>;

} 
} 

#endif 