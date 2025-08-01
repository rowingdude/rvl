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

#ifndef RVL_PROPAGATION_RAYLEIGH_FADING_HPP
#define RVL_PROPAGATION_RAYLEIGH_FADING_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <random>

namespace rvl {
namespace propagation {

template<typename T>
class rayleigh_fading {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    struct fading_parameters {
        T rms_delay_spread_s;
        T coherence_bandwidth_hz;
        T doppler_spread_hz;
        T coherence_time_s;
        T power_delay_profile_decay_db_per_us;
    };

    static T calculate_rayleigh_pdf(T envelope, T sigma) {
        core::check_positive(sigma, "Sigma parameter");
        core::check_non_negative(envelope, "Envelope");

        if (envelope <= T(0)) {
            return T(0);
        }

        const T sigma_squared = sigma * sigma;
        return (envelope / sigma_squared) * std::exp(-envelope * envelope / (T(2.0) * sigma_squared));
    }

    static T calculate_rayleigh_cdf(T envelope, T sigma) {
        core::check_positive(sigma, "Sigma parameter");
        core::check_non_negative(envelope, "Envelope");

        if (envelope <= T(0)) {
            return T(0);
        }

        const T sigma_squared = sigma * sigma;
        return T(1.0) - std::exp(-envelope * envelope / (T(2.0) * sigma_squared));
    }

    static T calculate_rayleigh_mean(T sigma) {
        core::check_positive(sigma, "Sigma parameter");
        return sigma * std::sqrt(constants::mathematical<T>::pi / T(2.0));
    }

    static T calculate_rayleigh_variance(T sigma) {
        core::check_positive(sigma, "Sigma parameter");
        const T sigma_squared = sigma * sigma;
        return (T(4.0) - constants::mathematical<T>::pi) * sigma_squared / T(2.0);
    }

    static T calculate_rayleigh_rms(T sigma) {
        core::check_positive(sigma, "Sigma parameter");
        return sigma * std::sqrt(T(2.0));
    }

    static T calculate_fade_depth_db(T probability) {
        core::check_range(probability, T(0), T(1), "Probability");

        if (probability >= T(1.0)) {
            return T(-300.0); 
        }

        if (probability <= T(1e-10)) {
            return T(0);
        }

        const T envelope = std::sqrt(-T(2.0) * std::log(T(1.0) - probability));
        return T(20.0) * std::log10(envelope);
    }

    static T calculate_outage_probability(T fade_depth_db) {

        const T envelope_ratio = std::pow(T(10.0), -fade_depth_db / T(20.0));

        if (envelope_ratio >= T(10.0)) {
            return T(1.0);
        }

        return T(1.0) - std::exp(-envelope_ratio * envelope_ratio / T(2.0));
    }

    static complex_vector_type generate_rayleigh_fading_sequence(size_t num_samples, T sigma, 
                                                               std::mt19937& rng) {
        core::check_positive(sigma, "Sigma parameter");

        complex_vector_type fading_sequence(num_samples);
        std::normal_distribution<T> normal_dist(T(0), sigma);

        for (size_t i = 0; i < num_samples; ++i) {
            const T real_part = normal_dist(rng);
            const T imag_part = normal_dist(rng);
            fading_sequence[i] = complex_type(real_part, imag_part);
        }

        return fading_sequence;
    }

    static complex_vector_type generate_correlated_rayleigh_fading(size_t num_samples, T sigma,
                                                                 T sample_rate_hz, T doppler_hz,
                                                                 std::mt19937& rng) {
        core::check_positive(sigma, "Sigma parameter");
        core::check_positive(sample_rate_hz, "Sample rate");
        core::check_positive(doppler_hz, "Doppler frequency");

        auto independent_sequence = generate_rayleigh_fading_sequence(num_samples, sigma, rng);

        const T fc = doppler_hz / sample_rate_hz; 

        if (fc >= T(0.5)) {
            return independent_sequence; 
        }

        const T alpha = T(2.0) * constants::mathematical<T>::pi * fc;
        const T filter_coeff = std::exp(-alpha);

        complex_vector_type filtered_sequence(num_samples);
        filtered_sequence[0] = independent_sequence[0];

        for (size_t i = 1; i < num_samples; ++i) {
            filtered_sequence[i] = filter_coeff * filtered_sequence[i-1] + 
                                 (T(1.0) - filter_coeff) * independent_sequence[i];
        }

        return filtered_sequence;
    }

    static fading_parameters calculate_fading_parameters_from_environment(const std::string& environment_type) {
        fading_parameters params;

        if (environment_type == "urban_microcell") {
            params.rms_delay_spread_s = T(0.5e-6);    
            params.coherence_bandwidth_hz = T(200000); 
            params.doppler_spread_hz = T(50.0);       
            params.coherence_time_s = T(0.01);        
            params.power_delay_profile_decay_db_per_us = T(20.0);
        } else if (environment_type == "urban_macrocell") {
            params.rms_delay_spread_s = T(2.0e-6);    
            params.coherence_bandwidth_hz = T(50000);  
            params.doppler_spread_hz = T(100.0);      
            params.coherence_time_s = T(0.005);       
            params.power_delay_profile_decay_db_per_us = T(15.0);
        } else if (environment_type == "suburban") {
            params.rms_delay_spread_s = T(1.0e-6);    
            params.coherence_bandwidth_hz = T(100000); 
            params.doppler_spread_hz = T(80.0);       
            params.coherence_time_s = T(0.006);       
            params.power_delay_profile_decay_db_per_us = T(18.0);
        } else if (environment_type == "rural") {
            params.rms_delay_spread_s = T(0.2e-6);    
            params.coherence_bandwidth_hz = T(500000); 
            params.doppler_spread_hz = T(120.0);      
            params.coherence_time_s = T(0.004);       
            params.power_delay_profile_decay_db_per_us = T(25.0);
        } else if (environment_type == "indoor") {
            params.rms_delay_spread_s = T(0.05e-6);   
            params.coherence_bandwidth_hz = T(2000000); 
            params.doppler_spread_hz = T(5.0);        
            params.coherence_time_s = T(0.1);         
            params.power_delay_profile_decay_db_per_us = T(30.0);
        } else {

            params.rms_delay_spread_s = T(1.0e-6);
            params.coherence_bandwidth_hz = T(100000);
            params.doppler_spread_hz = T(50.0);
            params.coherence_time_s = T(0.01);
            params.power_delay_profile_decay_db_per_us = T(20.0);
        }

        return params;
    }

    static T calculate_coherence_bandwidth_hz(T rms_delay_spread_s) {
        core::check_positive(rms_delay_spread_s, "RMS delay spread");
        return T(0.2) / rms_delay_spread_s; 
    }

    static T calculate_coherence_time_s(T doppler_spread_hz) {
        core::check_positive(doppler_spread_hz, "Doppler spread");
        return T(0.423) / doppler_spread_hz; 
    }

    static T calculate_level_crossing_rate_hz(T threshold_db, T doppler_hz, T average_power_db = T(0)) {
        core::check_positive(doppler_hz, "Doppler frequency");

        const T threshold_linear = std::pow(T(10.0), (threshold_db - average_power_db) / T(20.0));

        if (threshold_linear <= T(0)) {
            return doppler_hz * std::sqrt(T(2.0) * constants::mathematical<T>::pi);
        }

        const T rho_squared = threshold_linear * threshold_linear / T(2.0);

        return doppler_hz * std::sqrt(T(2.0) * constants::mathematical<T>::pi * rho_squared) * 
               std::exp(-rho_squared);
    }

    static T calculate_average_fade_duration_s(T threshold_db, T doppler_hz, T average_power_db = T(0)) {
        const T outage_prob = calculate_outage_probability(threshold_db - average_power_db);
        const T crossing_rate = calculate_level_crossing_rate_hz(threshold_db, doppler_hz, average_power_db);

        if (crossing_rate <= T(0) || outage_prob <= T(0)) {
            return T(0);
        }

        return outage_prob / crossing_rate;
    }

    static void calculate_fading_statistics_batch(const vector_type& envelope_values,
                                                T sigma,
                                                vector_type& pdf_values,
                                                vector_type& cdf_values) {
        const size_t n = envelope_values.size();
        if (pdf_values.size() != n || cdf_values.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            pdf_values[i] = calculate_rayleigh_pdf(envelope_values[i], sigma);
            cdf_values[i] = calculate_rayleigh_cdf(envelope_values[i], sigma);
        }
    }

    static T calculate_diversity_gain_db(int num_branches, T correlation_coefficient = T(0)) {
        core::check_positive(num_branches, "Number of diversity branches");
        core::check_range(correlation_coefficient, T(0), T(1), "Correlation coefficient");

        if (num_branches == 1) {
            return T(0);
        }

        T ideal_gain_db = T(10.0) * std::log10(T(num_branches));
        T correlation_loss_db = T(10.0) * std::log10(T(1.0) + correlation_coefficient * (T(num_branches) - T(1.0)));

        return ideal_gain_db - correlation_loss_db;
    }

    static bool is_fast_fading(T symbol_duration_s, T coherence_time_s) {
        return symbol_duration_s > coherence_time_s;
    }

    static bool is_frequency_selective(T signal_bandwidth_hz, T coherence_bandwidth_hz) {
        return signal_bandwidth_hz > coherence_bandwidth_hz;
    }

    static constexpr T typical_urban_rms_delay_spread_us() { return T(1.0); }
    static constexpr T typical_doppler_at_walking_speed_hz() { return T(5.0); }
    static constexpr T typical_doppler_at_vehicular_speed_hz() { return T(100.0); }
    static constexpr T rayleigh_distribution_mode() { return T(0); } 
};

using rayleigh_fading_f = rayleigh_fading<float>;
using rayleigh_fading_d = rayleigh_fading<double>;

} 
} 

#endif 