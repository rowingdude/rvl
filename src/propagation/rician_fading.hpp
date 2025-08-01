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

#ifndef RVL_PROPAGATION_RICIAN_FADING_HPP
#define RVL_PROPAGATION_RICIAN_FADING_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <random>

namespace rvl {
namespace propagation {

template<typename T>
class rician_fading {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    struct rician_parameters {
        T k_factor_db;           
        T total_power_db;        
        T doppler_frequency_hz;  
        T phase_los_rad;         
    };

    static T calculate_k_factor_linear(T k_factor_db) {
        return std::pow(T(10.0), k_factor_db / T(10.0));
    }

    static T calculate_k_factor_db(T k_factor_linear) {
        core::check_non_negative(k_factor_linear, "K-factor (linear)");
        return T(10.0) * std::log10(k_factor_linear);
    }

    static T calculate_rician_pdf(T envelope, T k_factor_linear, T omega) {
        core::check_non_negative(envelope, "Envelope");
        core::check_non_negative(k_factor_linear, "K-factor");
        core::check_positive(omega, "Omega parameter");

        if (envelope <= T(0)) {
            return T(0);
        }

        const T sigma_squared = omega / (T(2.0) * (T(1.0) + k_factor_linear));
        const T s = std::sqrt(k_factor_linear * omega);

        const T exp_term = std::exp(-(envelope * envelope + s * s) / (T(2.0) * sigma_squared));
        const T bessel_arg = envelope * s / sigma_squared;
        const T bessel_i0 = calculate_modified_bessel_i0(bessel_arg);

        return (envelope / sigma_squared) * exp_term * bessel_i0;
    }

    static T calculate_rician_cdf(T envelope, T k_factor_linear, T omega) {

        core::check_non_negative(envelope, "Envelope");
        core::check_non_negative(k_factor_linear, "K-factor");
        core::check_positive(omega, "Omega parameter");

        if (envelope <= T(0)) {
            return T(0);
        }

        const T sigma_squared = omega / (T(2.0) * (T(1.0) + k_factor_linear));
        const T a = std::sqrt(k_factor_linear * omega);
        const T x = envelope;

        const T alpha = a / std::sqrt(sigma_squared);
        const T beta = x / std::sqrt(sigma_squared);

        return calculate_marcum_q_function(alpha, beta);
    }

    static T calculate_modified_bessel_i0(T x) {

        if (std::abs(x) < T(3.75)) {
            const T t = x / T(3.75);
            const T t2 = t * t;
            return T(1.0) + T(3.5156229) * t2 + T(3.0899424) * t2 * t2 +
                   T(1.2067492) * t2 * t2 * t2 + T(0.2659732) * t2 * t2 * t2 * t2 +
                   T(0.0360768) * t2 * t2 * t2 * t2 * t2 + T(0.0045813) * t2 * t2 * t2 * t2 * t2 * t2;
        } else {
            const T t = T(3.75) / std::abs(x);
            const T exp_x = std::exp(std::abs(x));
            const T sqrt_x = std::sqrt(std::abs(x));

            const T polynomial = T(0.39894228) + T(0.01328592) * t + T(0.00225319) * t * t -
                               T(0.00157565) * t * t * t + T(0.00916281) * t * t * t * t -
                               T(0.02057706) * t * t * t * t * t + T(0.02635537) * t * t * t * t * t * t -
                               T(0.01647633) * t * t * t * t * t * t * t + T(0.00392377) * t * t * t * t * t * t * t * t;

            return (exp_x / sqrt_x) * polynomial;
        }
    }

    static T calculate_marcum_q_function(T alpha, T beta) {

        if (alpha <= T(0)) {

            return std::exp(-beta * beta / T(2.0));
        }

        if (beta <= T(0)) {
            return T(1.0);
        }

        if (alpha < T(5.0) && beta < T(5.0)) {

            T sum = T(0);
            const T alpha_squared = alpha * alpha;
            const T beta_squared = beta * beta;
            const T exp_term = std::exp(-(alpha_squared + beta_squared) / T(2.0));

            for (int n = 0; n < 50; ++n) {
                const T term = std::pow(alpha * beta / T(2.0), T(n)) * 
                              calculate_modified_bessel_i0(alpha * beta) / factorial_approx(n);
                sum += term;

                if (term < T(1e-10)) {
                    break;
                }
            }

            return T(1.0) - exp_term * sum;
        } else {

            const T u = (beta - alpha) / std::sqrt(T(2.0));
            return T(0.5) * (T(1.0) + std::erf(u / std::sqrt(T(2.0))));
        }
    }

    static T factorial_approx(int n) {
        if (n <= 1) return T(1.0);
        if (n <= 12) {
            T result = T(1.0);
            for (int i = 2; i <= n; ++i) {
                result *= T(i);
            }
            return result;
        }

        const T n_real = T(n);
        return std::sqrt(T(2.0) * constants::mathematical<T>::pi * n_real) * 
               std::pow(n_real / constants::mathematical<T>::euler, n_real);
    }

    static T calculate_rician_mean(T k_factor_linear, T omega) {
        core::check_non_negative(k_factor_linear, "K-factor");
        core::check_positive(omega, "Omega parameter");

        const T sigma = std::sqrt(omega / (T(2.0) * (T(1.0) + k_factor_linear)));
        const T laguerre_half = T(1.25331414); 

        return sigma * laguerre_half * std::sqrt(constants::mathematical<T>::pi / T(2.0));
    }

    static T calculate_rician_variance(T k_factor_linear, T omega) {
        core::check_non_negative(k_factor_linear, "K-factor");
        core::check_positive(omega, "Omega parameter");

        const T mean_squared = calculate_rician_mean(k_factor_linear, omega);
        return omega - mean_squared * mean_squared;
    }

    static complex_vector_type generate_rician_fading_sequence(size_t num_samples, 
                                                             const rician_parameters& params,
                                                             std::mt19937& rng) {
        const T k_linear = calculate_k_factor_linear(params.k_factor_db);
        const T total_power_linear = std::pow(T(10.0), params.total_power_db / T(10.0));

        const T los_power = k_linear * total_power_linear / (T(1.0) + k_linear);
        const T scattered_power = total_power_linear / (T(1.0) + k_linear);

        const T sigma = std::sqrt(scattered_power / T(2.0));
        const T los_amplitude = std::sqrt(los_power);

        complex_vector_type fading_sequence(num_samples);
        std::normal_distribution<T> normal_dist(T(0), sigma);

        const complex_type los_component = std::polar(los_amplitude, params.phase_los_rad);

        for (size_t i = 0; i < num_samples; ++i) {
            const T real_scattered = normal_dist(rng);
            const T imag_scattered = normal_dist(rng);
            const complex_type scattered_component(real_scattered, imag_scattered);

            fading_sequence[i] = los_component + scattered_component;
        }

        return fading_sequence;
    }

    static rician_parameters estimate_rician_parameters(const complex_vector_type& signal_samples) {
        if (signal_samples.empty()) {
            throw core::invalid_argument_error("Signal samples cannot be empty");
        }

        rician_parameters params;

        T sum_magnitude = T(0);
        T sum_magnitude_squared = T(0);
        complex_type sum_samples(T(0), T(0));

        for (const auto& sample : signal_samples) {
            const T magnitude = std::abs(sample);
            sum_magnitude += magnitude;
            sum_magnitude_squared += magnitude * magnitude;
            sum_samples += sample;
        }

        const size_t n = signal_samples.size();
        const T mean_magnitude = sum_magnitude / T(n);
        const T mean_power = sum_magnitude_squared / T(n);
        const complex_type mean_complex = sum_samples / T(n);

        const T los_power = std::abs(mean_complex) * std::abs(mean_complex);
        const T scattered_power = mean_power - los_power;

        if (scattered_power > T(0)) {
            const T k_linear = los_power / scattered_power;
            params.k_factor_db = calculate_k_factor_db(k_linear);
        } else {
            params.k_factor_db = T(40.0); 
        }

        params.total_power_db = T(10.0) * std::log10(mean_power);
        params.phase_los_rad = std::arg(mean_complex);
        params.doppler_frequency_hz = T(0); 

        return params;
    }

    static T calculate_outage_probability_rician(T threshold_db, T k_factor_db, T average_power_db) {
        const T k_linear = calculate_k_factor_linear(k_factor_db);
        const T threshold_power = std::pow(T(10.0), (threshold_db - average_power_db) / T(10.0));
        const T omega = T(1.0); 

        if (threshold_power <= T(0.01)) {
            return std::exp(-k_linear) * std::pow(threshold_power, T(1.0) + k_linear);
        }

        const T envelope_threshold = std::sqrt(threshold_power);
        return T(1.0) - calculate_rician_cdf(envelope_threshold, k_linear, omega);
    }

    static T calculate_diversity_gain_rician_db(int num_branches, T k_factor_db, T correlation = T(0)) {
        const T k_linear = calculate_k_factor_linear(k_factor_db);

        T rayleigh_gain = T(10.0) * std::log10(T(num_branches));
        T k_factor_reduction = T(10.0) * std::log10(T(1.0) + k_linear) / T(num_branches);
        T correlation_loss = T(10.0) * std::log10(T(1.0) + correlation * (T(num_branches) - T(1.0)));

        return rayleigh_gain - k_factor_reduction - correlation_loss;
    }

    static void calculate_rician_statistics_batch(const vector_type& envelope_values,
                                                T k_factor_db, T total_power_db,
                                                vector_type& pdf_values,
                                                vector_type& cdf_values) {
        const size_t n = envelope_values.size();
        if (pdf_values.size() != n || cdf_values.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        const T k_linear = calculate_k_factor_linear(k_factor_db);
        const T omega = std::pow(T(10.0), total_power_db / T(10.0));

        for (size_t i = 0; i < n; ++i) {
            pdf_values[i] = calculate_rician_pdf(envelope_values[i], k_linear, omega);
            cdf_values[i] = calculate_rician_cdf(envelope_values[i], k_linear, omega);
        }
    }

    static bool is_dominant_los_condition(T k_factor_db) {
        return k_factor_db > T(10.0); 
    }

    static T convert_to_rayleigh_equivalent_sigma(T k_factor_db, T total_power_db) {

        const T k_linear = calculate_k_factor_linear(k_factor_db);
        const T total_power_linear = std::pow(T(10.0), total_power_db / T(10.0));
        const T scattered_power = total_power_linear / (T(1.0) + k_linear);

        return std::sqrt(scattered_power / T(2.0));
    }

    static constexpr T typical_urban_k_factor_db() { return T(3.0); }
    static constexpr T typical_suburban_k_factor_db() { return T(7.0); }
    static constexpr T typical_rural_k_factor_db() { return T(15.0); }
    static constexpr T minimum_k_factor_db() { return T(-10.0); } 
    static constexpr T maximum_practical_k_factor_db() { return T(30.0); }
};

using rician_fading_f = rician_fading<float>;
using rician_fading_d = rician_fading<double>;

} 
} 

#endif 