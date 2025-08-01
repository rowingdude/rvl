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

#ifndef RVL_PROPAGATION_GROUND_WAVE_HPP
#define RVL_PROPAGATION_GROUND_WAVE_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>

namespace rvl {
namespace propagation {

template<typename T>
class ground_wave {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;

    struct ground_parameters {
        T relative_permittivity;
        T conductivity_s_per_m;
        T roughness_factor;
    };

    static complex_type calculate_ground_wave_attenuation_factor(T distance_m, T frequency_hz,
                                                               const ground_parameters& ground) {
        core::check_positive(distance_m, "Distance");
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(ground.relative_permittivity, "Relative permittivity");
        core::check_non_negative(ground.conductivity_s_per_m, "Ground conductivity");

        const T c = constants::physical<T>::c;
        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;
        const T wavelength = c / frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;

        const T epsilon_0 = constants::physical<T>::epsilon_0;
        const complex_type epsilon_r(ground.relative_permittivity, -ground.conductivity_s_per_m / (omega * epsilon_0));

        const complex_type numerical_distance = k * distance_m * std::sqrt(epsilon_r - T(1.0)) / T(2.0);

        if (std::abs(numerical_distance) < T(0.1)) {
            return complex_type(T(1.0), T(0));
        }

        const complex_type sqrt_pi_z = std::sqrt(constants::mathematical<T>::pi * numerical_distance);
        const complex_type erfc_arg = sqrt_pi_z;

        const complex_type w_function = calculate_complex_error_function(erfc_arg);
        const complex_type attenuation_factor = T(2.0) * sqrt_pi_z * std::exp(-numerical_distance) * w_function;

        return T(1.0) + attenuation_factor;
    }

    static complex_type calculate_complex_error_function(const complex_type& z) {
        const T real_part = z.real();
        const T imag_part = z.imag();

        if (std::abs(z) < T(1.0)) {
            complex_type sum = z;
            complex_type term = z;

            for (int n = 1; n <= 50; ++n) {
                term *= -z * z / T(n);
                const complex_type factorial_term = term / T(2 * n + 1);
                sum += factorial_term;

                if (std::abs(factorial_term) < T(1e-15)) {
                    break;
                }
            }

            return T(2.0) / std::sqrt(constants::mathematical<T>::pi) * sum;
        } else {
            const complex_type continued_fraction = calculate_continued_fraction_erfc(z);
            return continued_fraction;
        }
    }

    static complex_type calculate_continued_fraction_erfc(const complex_type& z) {
        const int max_iterations = 100;
        const T tolerance = T(1e-15);

        complex_type b = z * z + T(0.5);
        complex_type c = T(1e30);
        complex_type d = T(1.0) / b;
        complex_type h = d;

        for (int i = 1; i <= max_iterations; ++i) {
            const T a = -T(i) + T(0.5);
            b += T(2.0);
            d = T(1.0) / (a * d + b);
            c = b + a / c;
            const complex_type del = c * d;
            h *= del;

            if (std::abs(del - T(1.0)) < tolerance) {
                break;
            }
        }

        const complex_type gamma = std::sqrt(constants::mathematical<T>::pi) / (z * std::exp(z * z));
        return gamma * h;
    }

    static T calculate_ground_wave_field_strength_db(T transmit_power_w, T frequency_hz, T distance_m,
                                                   const ground_parameters& ground,
                                                   T transmitter_height_m = T(1.0),
                                                   T receiver_height_m = T(1.0)) {
        core::check_positive(transmit_power_w, "Transmit power");
        core::check_positive(transmitter_height_m, "Transmitter height");
        core::check_positive(receiver_height_m, "Receiver height");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T free_space_field = std::sqrt(T(30.0) * transmit_power_w) / distance_m;

        const complex_type attenuation_factor = calculate_ground_wave_attenuation_factor(distance_m, frequency_hz, ground);

        const T height_factor = T(2.0) * std::sin(T(2.0) * constants::mathematical<T>::pi * 
                                                transmitter_height_m * receiver_height_m / 
                                                (wavelength * distance_m));

        const T field_strength = free_space_field * std::abs(attenuation_factor) * std::max(height_factor, T(0.1));

        return T(20.0) * std::log10(field_strength * T(1e6));
    }

    static T calculate_ground_wave_path_loss_db(T frequency_hz, T distance_m, const ground_parameters& ground) {
        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T free_space_loss = T(20.0) * std::log10(T(4.0) * constants::mathematical<T>::pi * distance_m / wavelength);

        const complex_type attenuation_factor = calculate_ground_wave_attenuation_factor(distance_m, frequency_hz, ground);
        const T ground_loss = -T(20.0) * std::log10(std::abs(attenuation_factor));

        return free_space_loss + ground_loss;
    }

    static T calculate_surface_wave_range_km(T frequency_hz, T transmit_power_dbm, 
                                           T receiver_sensitivity_dbm, const ground_parameters& ground) {
        core::check_positive(frequency_hz, "Frequency");

        const T available_margin_db = transmit_power_dbm - receiver_sensitivity_dbm;

        T test_distance = T(1000.0);
        const T step_factor = T(1.2);

        for (int iteration = 0; iteration < 100; ++iteration) {
            const T path_loss = calculate_ground_wave_path_loss_db(frequency_hz, test_distance, ground);

            if (path_loss > available_margin_db) {
                return test_distance / T(1000.0);
            }

            test_distance *= step_factor;

            if (test_distance > T(1000000.0)) {
                break;
            }
        }

        return test_distance / T(1000.0);
    }

    static T calculate_ground_conductivity_from_soil_type(const std::string& soil_type) {
        if (soil_type == "seawater") return T(5.0);
        else if (soil_type == "freshwater") return T(0.01);
        else if (soil_type == "wet_ground") return T(0.02);
        else if (soil_type == "medium_dry_ground") return T(0.005);
        else if (soil_type == "dry_ground") return T(0.001);
        else if (soil_type == "very_dry_ground") return T(0.0002);
        else return T(0.005);
    }

    static T calculate_relative_permittivity_from_soil_type(const std::string& soil_type) {
        if (soil_type == "seawater") return T(81.0);
        else if (soil_type == "freshwater") return T(81.0);
        else if (soil_type == "wet_ground") return T(30.0);
        else if (soil_type == "medium_dry_ground") return T(15.0);
        else if (soil_type == "dry_ground") return T(4.0);
        else if (soil_type == "very_dry_ground") return T(3.0);
        else return T(15.0);
    }

    static ground_parameters create_ground_parameters(const std::string& soil_type, T roughness_factor = T(1.0)) {
        ground_parameters params;
        params.conductivity_s_per_m = calculate_ground_conductivity_from_soil_type(soil_type);
        params.relative_permittivity = calculate_relative_permittivity_from_soil_type(soil_type);
        params.roughness_factor = roughness_factor;
        return params;
    }

    static T calculate_skip_distance_km(T frequency_mhz, T critical_frequency_mhz, T layer_height_km) {
        core::check_positive(frequency_mhz, "Frequency");
        core::check_positive(critical_frequency_mhz, "Critical frequency");
        core::check_positive(layer_height_km, "Layer height");

        if (frequency_mhz <= critical_frequency_mhz) {
            return T(0);
        }

        const T muf_factor = frequency_mhz / critical_frequency_mhz;
        const T sin_angle = T(1.0) / muf_factor;

        if (sin_angle >= T(1.0)) {
            return T(0);
        }

        const T angle_rad = std::asin(sin_angle);
        const T earth_radius = constants::physical<T>::earth_radius / T(1000.0);

        return T(2.0) * earth_radius * std::sin(angle_rad) * std::cos(angle_rad) * layer_height_km / earth_radius;
    }

    static void calculate_path_losses_batch(const vector_type& frequencies,
                                          const vector_type& distances,
                                          const vector_type& conductivities,
                                          const vector_type& permittivities,
                                          vector_type& path_losses) {
        const size_t n = frequencies.size();
        if (distances.size() != n || conductivities.size() != n || 
            permittivities.size() != n || path_losses.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            ground_parameters ground;
            ground.conductivity_s_per_m = conductivities[i];
            ground.relative_permittivity = permittivities[i];
            ground.roughness_factor = T(1.0);

            path_losses[i] = calculate_ground_wave_path_loss_db(frequencies[i], distances[i], ground);
        }
    }

    static T calculate_ground_wave_delay_excess_ns(T distance_m, T frequency_hz, const ground_parameters& ground) {
        const T c = constants::physical<T>::c;
        const T free_space_delay = distance_m / c;

        const complex_type attenuation_factor = calculate_ground_wave_attenuation_factor(distance_m, frequency_hz, ground);
        const T phase_delay_factor = std::arg(attenuation_factor) / (T(2.0) * constants::mathematical<T>::pi);

        const T wavelength = c / frequency_hz;
        const T excess_delay = phase_delay_factor * wavelength / c;

        return excess_delay * T(1e9);
    }

    static bool is_ground_wave_dominant(T frequency_mhz, T distance_km) {
        if (frequency_mhz > T(30.0)) {
            return false;
        }

        const T vhf_threshold_km = T(50.0) * std::exp(-frequency_mhz / T(5.0));
        return distance_km < vhf_threshold_km;
    }

    static constexpr ground_parameters seawater_ground() {
        return {T(81.0), T(5.0), T(1.0)};
    }

    static constexpr ground_parameters average_ground() {
        return {T(15.0), T(0.005), T(1.0)};
    }

    static constexpr ground_parameters poor_ground() {
        return {T(4.0), T(0.001), T(1.0)};
    }
};

using ground_wave_f = ground_wave<float>;
using ground_wave_d = ground_wave<double>;

} 
} 

#endif 