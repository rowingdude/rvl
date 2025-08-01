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

#ifndef RVL_ANTENNA_ARRAY_FACTOR_HPP
#define RVL_ANTENNA_ARRAY_FACTOR_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <vector>
#include <array>

namespace rvl {
namespace antenna {

template<typename T>
class array_factor {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    struct array_element {
        T x_m;                   
        T y_m;                   
        T z_m;                   
        complex_type amplitude;  
        T element_gain_db;       
    };

    struct linear_array_params {
        size_t num_elements;     
        T element_spacing_m;     
        T frequency_hz;          
        T progressive_phase_rad; 
        T beam_steering_angle_rad; 
    };

    struct planar_array_params {
        size_t num_x_elements;   
        size_t num_y_elements;   
        T spacing_x_m;           
        T spacing_y_m;           
        T frequency_hz;          
        T steering_theta_rad;    
        T steering_phi_rad;      
    };

    struct circular_array_params {
        size_t num_elements;     
        T radius_m;              
        T frequency_hz;          
        T beam_steering_angle_rad; 
    };

    static complex_type calculate_linear_array_factor(T observation_angle_rad,
                                                    const linear_array_params& params,
                                                    const complex_vector_type& element_amplitudes = {}) {
        core::check_positive(params.num_elements, "Number of elements");
        core::check_positive(params.frequency_hz, "Frequency");

        const T wavelength = constants::physical<T>::c / params.frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T kd_cos_theta = k * params.element_spacing_m * std::cos(observation_angle_rad);

        complex_type array_factor(T(0.0), T(0.0));

        const bool use_provided_amplitudes = !element_amplitudes.empty() && 
                                           element_amplitudes.size() == params.num_elements;

        for (size_t n = 0; n < params.num_elements; ++n) {
            const T phase = T(n) * (kd_cos_theta + params.progressive_phase_rad);
            const complex_type phase_term = std::exp(complex_type(T(0.0), phase));

            complex_type amplitude(T(1.0), T(0.0)); 
            if (use_provided_amplitudes) {
                amplitude = element_amplitudes[n];
            }

            array_factor += amplitude * phase_term;
        }

        return array_factor;
    }

    static complex_type calculate_planar_array_factor(T theta_rad, T phi_rad,
                                                    const planar_array_params& params,
                                                    const std::vector<complex_vector_type>& element_amplitudes = {}) {
        core::check_positive(params.num_x_elements, "Number of X elements");
        core::check_positive(params.num_y_elements, "Number of Y elements");
        core::check_positive(params.frequency_hz, "Frequency");

        const T wavelength = constants::physical<T>::c / params.frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;

        const T sin_theta = std::sin(theta_rad);
        const T cos_theta = std::cos(theta_rad);
        const T sin_phi = std::sin(phi_rad);
        const T cos_phi = std::cos(phi_rad);

        const T kx = k * sin_theta * cos_phi;
        const T ky = k * sin_theta * sin_phi;

        const T kx_steer = k * std::sin(params.steering_theta_rad) * std::cos(params.steering_phi_rad);
        const T ky_steer = k * std::sin(params.steering_theta_rad) * std::sin(params.steering_phi_rad);

        complex_type array_factor(T(0.0), T(0.0));

        const bool use_provided_amplitudes = !element_amplitudes.empty() && 
                                           element_amplitudes.size() == params.num_y_elements &&
                                           element_amplitudes[0].size() == params.num_x_elements;

        for (size_t m = 0; m < params.num_y_elements; ++m) {
            for (size_t n = 0; n < params.num_x_elements; ++n) {
                const T x_pos = T(n) * params.spacing_x_m;
                const T y_pos = T(m) * params.spacing_y_m;

                const T phase = (kx - kx_steer) * x_pos + (ky - ky_steer) * y_pos;
                const complex_type phase_term = std::exp(complex_type(T(0.0), phase));

                complex_type amplitude(T(1.0), T(0.0)); 
                if (use_provided_amplitudes) {
                    amplitude = element_amplitudes[m][n];
                }

                array_factor += amplitude * phase_term;
            }
        }

        return array_factor;
    }

    static complex_type calculate_circular_array_factor(T observation_angle_rad,
                                                       const circular_array_params& params,
                                                       const complex_vector_type& element_amplitudes = {}) {
        core::check_positive(params.num_elements, "Number of elements");
        core::check_positive(params.radius_m, "Array radius");
        core::check_positive(params.frequency_hz, "Frequency");

        const T wavelength = constants::physical<T>::c / params.frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T kr = k * params.radius_m;

        complex_type array_factor(T(0.0), T(0.0));

        const bool use_provided_amplitudes = !element_amplitudes.empty() && 
                                           element_amplitudes.size() == params.num_elements;

        for (size_t n = 0; n < params.num_elements; ++n) {
            const T element_angle = T(2.0) * constants::mathematical<T>::pi * T(n) / T(params.num_elements);
            const T phase_diff = observation_angle_rad - element_angle;
            const T phase = kr * std::cos(phase_diff) - T(n) * params.beam_steering_angle_rad;

            const complex_type phase_term = std::exp(complex_type(T(0.0), phase));

            complex_type amplitude(T(1.0), T(0.0)); 
            if (use_provided_amplitudes) {
                amplitude = element_amplitudes[n];
            }

            array_factor += amplitude * phase_term;
        }

        return array_factor;
    }

    static complex_vector_type generate_uniform_amplitudes(size_t num_elements) {
        complex_vector_type amplitudes(num_elements);
        const complex_type unit_amplitude(T(1.0), T(0.0));

        for (size_t i = 0; i < num_elements; ++i) {
            amplitudes[i] = unit_amplitude;
        }

        return amplitudes;
    }

    static complex_vector_type generate_dolph_chebyshev_amplitudes(size_t num_elements,
                                                                 T sidelobe_level_db) {
        core::check_positive(num_elements, "Number of elements");
        core::check_range(sidelobe_level_db, T(-60.0), T(0.0), "Sidelobe level");

        complex_vector_type amplitudes(num_elements);

        const T R = std::pow(T(10.0), -sidelobe_level_db / T(20.0));
        const T x0 = std::cosh(std::acosh(R) / T(num_elements - 1));

        for (size_t n = 0; n < num_elements; ++n) {
            T sum = T(0.0);

            for (size_t k = 1; k <= (num_elements - 1) / 2; ++k) {
                const T x = x0 * std::cos(constants::mathematical<T>::pi * T(k) / T(num_elements));
                const T T_k = std::cos(T(num_elements - 1) * std::acos(x));
                const T cos_term = std::cos(T(2.0) * constants::mathematical<T>::pi * T(k) * T(n) / T(num_elements));
                sum += T_k * cos_term;
            }

            amplitudes[n] = complex_type(T(1.0) + T(2.0) * sum / T(num_elements), T(0.0));
        }

        T max_amplitude = T(0.0);
        for (const auto& amp : amplitudes) {
            max_amplitude = std::max(max_amplitude, std::abs(amp));
        }

        for (auto& amp : amplitudes) {
            amp /= max_amplitude;
        }

        return amplitudes;
    }

    static complex_vector_type generate_taylor_amplitudes(size_t num_elements,
                                                        T sidelobe_level_db,
                                                        size_t num_equal_sidelobes = 5) {
        core::check_positive(num_elements, "Number of elements");
        core::check_range(sidelobe_level_db, T(-60.0), T(0.0), "Sidelobe level");

        complex_vector_type amplitudes(num_elements);

        const T R = std::pow(T(10.0), -sidelobe_level_db / T(20.0));
        const T A = std::acosh(R) / constants::mathematical<T>::pi;
        const T sigma = T(num_equal_sidelobes) / std::sqrt(A * A + T(num_equal_sidelobes - 0.5) * T(num_equal_sidelobes - 0.5));

        for (size_t n = 0; n < num_elements; ++n) {
            const T m = T(n) - T(num_elements - 1) / T(2.0);
            const T x = m / (T(num_elements) / T(2.0));

            T amplitude = T(1.0);

            for (size_t k = 1; k <= num_equal_sidelobes; ++k) {
                const T F_k = std::sqrt((A * A + T(k - 0.5) * T(k - 0.5)) / 
                                      (A * A + T(k) * T(k)));
                const T numerator = T(1.0) - x * x / (sigma * sigma * (A * A + T(k - 0.5) * T(k - 0.5)));
                amplitude *= numerator * F_k;
            }

            amplitudes[n] = complex_type(amplitude, T(0.0));
        }

        return amplitudes;
    }

    static complex_vector_type generate_binomial_amplitudes(size_t num_elements) {
        core::check_positive(num_elements, "Number of elements");

        complex_vector_type amplitudes(num_elements);

        for (size_t n = 0; n < num_elements; ++n) {
            T binomial_coeff = T(1.0);

            for (size_t k = 1; k <= n; ++k) {
                binomial_coeff *= T(num_elements - 1 - k + 1) / T(k);
            }

            amplitudes[n] = complex_type(binomial_coeff, T(0.0));
        }

        T sum_magnitude = T(0.0);
        for (const auto& amp : amplitudes) {
            sum_magnitude += std::abs(amp);
        }

        for (auto& amp : amplitudes) {
            amp /= sum_magnitude;
        }

        return amplitudes;
    }

    static std::vector<T> calculate_grating_lobe_angles(const linear_array_params& params) {
        std::vector<T> grating_angles;

        const T wavelength = constants::physical<T>::c / params.frequency_hz;
        const T d_over_lambda = params.element_spacing_m / wavelength;

        if (d_over_lambda <= T(1.0)) {
            return grating_angles; 
        }

        for (int n = -5; n <= 5; ++n) {
            if (n == 0) continue; 

            const T sin_theta = T(n) / d_over_lambda;

            if (std::abs(sin_theta) <= T(1.0)) {
                const T theta = std::asin(sin_theta);
                grating_angles.push_back(theta);
            }
        }

        return grating_angles;
    }

    static T calculate_array_factor_db(T observation_angle_rad,
                                     const linear_array_params& params,
                                     const complex_vector_type& element_amplitudes = {}) {
        const complex_type af = calculate_linear_array_factor(observation_angle_rad, params, element_amplitudes);
        const T magnitude = std::abs(af);

        if (magnitude <= T(0)) {
            return T(-100.0); 
        }

        return T(20.0) * std::log10(magnitude);
    }

    static T calculate_array_directivity_db(const linear_array_params& params,
                                          const complex_vector_type& element_amplitudes = {}) {

        const T wavelength = constants::physical<T>::c / params.frequency_hz;
        const T d_over_lambda = params.element_spacing_m / wavelength;

        T directivity_db;

        if (d_over_lambda <= T(0.5)) {

            directivity_db = T(10.0) * std::log10(T(params.num_elements));
        } else {

            directivity_db = T(10.0) * std::log10(T(params.num_elements)) + 
                           T(10.0) * std::log10(T(2.0) * d_over_lambda);
        }

        return directivity_db;
    }

    static T calculate_array_beamwidth_rad(const linear_array_params& params) {
        const T wavelength = constants::physical<T>::c / params.frequency_hz;
        const T array_length = T(params.num_elements - 1) * params.element_spacing_m;

        const T beamwidth_rad = T(0.886) * wavelength / array_length;

        return beamwidth_rad;
    }

    static void calculate_array_pattern_batch(const linear_array_params& params,
                                            const vector_type& observation_angles,
                                            vector_type& pattern_db,
                                            const complex_vector_type& element_amplitudes = {}) {
        if (observation_angles.size() != pattern_db.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = observation_angles.size();
        for (size_t i = 0; i < n; ++i) {
            pattern_db[i] = calculate_array_factor_db(observation_angles[i], params, element_amplitudes);
        }
    }

    static constexpr T maximum_element_spacing_wavelengths() { return constants::antenna_design<T>::maximum_element_spacing_wavelengths; }
    static constexpr T typical_element_spacing_wavelengths() { return constants::antenna_design<T>::typical_element_spacing_wavelengths; }
    static constexpr T minimum_element_spacing_wavelengths() { return constants::antenna_design<T>::minimum_element_spacing_wavelengths; }
    static constexpr size_t maximum_practical_elements() { return constants::antenna_design<T>::maximum_practical_elements; }
    static constexpr T typical_sidelobe_level_db() { return constants::antenna_design<T>::typical_sidelobe_level_db; }
    static constexpr T low_sidelobe_level_db() { return constants::antenna_design<T>::low_sidelobe_level_db; }
};

using array_factor_f = array_factor<float>;
using array_factor_d = array_factor<double>;

} 
} 

#endif 