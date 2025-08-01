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

#ifndef RVL_PROPAGATION_IONOSPHERIC_SCINTILLATION_HPP
#define RVL_PROPAGATION_IONOSPHERIC_SCINTILLATION_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <random>

namespace rvl {
namespace propagation {

template<typename T>
class ionospheric_scintillation {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    static T calculate_rytov_variance(T structure_constant, T path_length, T wavenumber) {
        core::check_positive(structure_constant, "Structure constant");
        core::check_positive(path_length, "Path length");
        core::check_positive(wavenumber, "Wavenumber");

        const T k_power = std::pow(wavenumber, T(7.0/6.0));
        const T L_power = std::pow(path_length, T(5.0/3.0));

        return T(1.23) * structure_constant * L_power * k_power;
    }

    static T calculate_rytov_variance_frequency(T structure_constant, T path_length, T frequency) {
        core::check_positive(frequency, "Frequency");

        const T c = constants::physical<T>::c;
        const T wavenumber = T(2.0) * constants::mathematical<T>::pi * frequency / c;

        return calculate_rytov_variance(structure_constant, path_length, wavenumber);
    }

    static T scintillation_index_s4(const complex_vector_type& electric_field) {
        if (electric_field.empty()) {
            throw core::invalid_argument_error("Electric field vector cannot be empty");
        }

        const size_t n = electric_field.size();

        T sum_magnitude_squared = T(0);
        T sum_magnitude_fourth = T(0);
        T sum_magnitude = T(0);

        for (size_t i = 0; i < n; ++i) {
            const T magnitude = std::abs(electric_field[i]);
            const T magnitude_squared = magnitude * magnitude;

            sum_magnitude += magnitude;
            sum_magnitude_squared += magnitude_squared;
            sum_magnitude_fourth += magnitude_squared * magnitude_squared;
        }

        const T mean_magnitude = sum_magnitude / T(n);
        const T mean_magnitude_squared = sum_magnitude_squared / T(n);
        const T mean_magnitude_fourth = sum_magnitude_fourth / T(n);

        const T variance_magnitude_squared = mean_magnitude_fourth - mean_magnitude_squared * mean_magnitude_squared;
        const T mean_squared_magnitude_squared = mean_magnitude_squared * mean_magnitude_squared;

        if (mean_squared_magnitude_squared <= T(0)) {
            return T(0);
        }

        return std::sqrt(variance_magnitude_squared / mean_squared_magnitude_squared);
    }

    static T phase_scintillation_index(const complex_vector_type& electric_field) {
        if (electric_field.size() < 2) {
            throw core::invalid_argument_error("Need at least 2 field samples");
        }

        const size_t n = electric_field.size();
        T sum_phase_variance = T(0);
        T sum_count = T(0);

        for (size_t i = 1; i < n; ++i) {
            const T phase1 = std::arg(electric_field[i-1]);
            const T phase2 = std::arg(electric_field[i]);

            T phase_diff = phase2 - phase1;

            while (phase_diff > constants::mathematical<T>::pi) {
                phase_diff -= T(2.0) * constants::mathematical<T>::pi;
            }
            while (phase_diff < -constants::mathematical<T>::pi) {
                phase_diff += T(2.0) * constants::mathematical<T>::pi;
            }

            sum_phase_variance += phase_diff * phase_diff;
            sum_count += T(1.0);
        }

        return std::sqrt(sum_phase_variance / sum_count);
    }

    static T fresnel_scale_length(T path_length, T wavenumber) {
        core::check_positive(path_length, "Path length");
        core::check_positive(wavenumber, "Wavenumber");

        return std::sqrt(T(2.0) * path_length / wavenumber);
    }

    static T weak_scattering_threshold() {
        return T(1.0);
    }

    static bool is_weak_scattering_regime(T rytov_variance) {
        return rytov_variance < weak_scattering_threshold();
    }

    static T strong_scattering_s4_saturated() {
        return T(1.0);
    }

    static T transition_regime_s4(T rytov_variance) {
        if (rytov_variance <= T(0)) {
            return T(0);
        }

        if (rytov_variance < weak_scattering_threshold()) {
            return std::sqrt(rytov_variance);
        } else {
            const T saturation_level = strong_scattering_s4_saturated();
            const T transition_factor = T(1.0) / (T(1.0) + rytov_variance);
            return saturation_level * (T(1.0) - transition_factor);
        }
    }

    static complex_vector_type generate_phase_screen(size_t num_points, T correlation_length, 
                                                   T phase_variance, std::mt19937& rng) {
        if (num_points == 0) {
            throw core::invalid_argument_error("Number of points must be positive");
        }

        core::check_positive(correlation_length, "Correlation length");
        core::check_positive(phase_variance, "Phase variance");

        complex_vector_type phase_screen(num_points);
        std::normal_distribution<T> normal_dist(T(0), std::sqrt(phase_variance));

        for (size_t i = 0; i < num_points; ++i) {
            const T real_part = normal_dist(rng);
            const T imag_part = normal_dist(rng);
            phase_screen[i] = std::exp(complex_type(T(0), real_part));
        }

        return phase_screen;
    }

    static T calculate_decorrelation_time(T velocity, T correlation_length) {
        core::check_positive(velocity, "Velocity");
        core::check_positive(correlation_length, "Correlation length");

        return correlation_length / velocity;
    }

    static T calculate_scintillation_bandwidth(T decorrelation_time) {
        core::check_positive(decorrelation_time, "Decorrelation time");

        return T(1.0) / (T(2.0) * constants::mathematical<T>::pi * decorrelation_time);
    }

    static void calculate_s4_batch(const vector_type& structure_constants,
                                 const vector_type& path_lengths,
                                 const vector_type& frequencies,
                                 vector_type& s4_indices) {
        const size_t n = structure_constants.size();
        if (path_lengths.size() != n || frequencies.size() != n || s4_indices.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (structure_constants[i] <= T(0) || path_lengths[i] <= T(0) || frequencies[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            const T rytov_var = calculate_rytov_variance_frequency(structure_constants[i], 
                                                                path_lengths[i], 
                                                                frequencies[i]);
            s4_indices[i] = transition_regime_s4(rytov_var);
        }
    }

    static T estimate_structure_constant_equatorial(T solar_activity_index = T(100.0)) {
        core::check_non_negative(solar_activity_index, "Solar activity index");

        const T base_ck = T(1e-50);
        const T activity_factor = T(1.0) + solar_activity_index / T(200.0);

        return base_ck * activity_factor;
    }

    static T estimate_structure_constant_auroral(T geomagnetic_index = T(3.0)) {
        core::check_non_negative(geomagnetic_index, "Geomagnetic index");

        const T base_ck = T(5e-50);
        const T magnetic_factor = T(1.0) + geomagnetic_index / T(5.0);

        return base_ck * magnetic_factor;
    }

    static constexpr T typical_structure_constant() { return T(1e-50); }
    static constexpr T typical_correlation_length_km() { return T(1.0); }
    static constexpr T typical_s4_strong() { return T(0.8); }
    static constexpr T typical_s4_weak() { return T(0.1); }
};

using ionospheric_scintillation_f = ionospheric_scintillation<float>;
using ionospheric_scintillation_d = ionospheric_scintillation<double>;

} 
} 

#endif 