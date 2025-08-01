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

#ifndef RVL_ANTENNA_CURRENT_DISTRIBUTION_HPP
#define RVL_ANTENNA_CURRENT_DISTRIBUTION_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <vector>

namespace rvl {
namespace antenna {

template<typename T>
class current_distribution {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    struct wire_geometry {
        T length_m;              
        T radius_m;              
        T frequency_hz;          
        complex_type feed_impedance; 
        T feed_position_m;       
    };

    struct current_sample {
        T position_m;            
        complex_type current;    
        T magnitude;             
        T phase_rad;             
    };

    static complex_type calculate_dipole_current(T position_m, const wire_geometry& wire) {
        core::check_positive(wire.length_m, "Wire length");
        core::check_positive(wire.frequency_hz, "Frequency");
        core::check_range(position_m, -wire.length_m/T(2.0), wire.length_m/T(2.0), "Position along wire");

        const T wavelength = constants::physical<T>::c / wire.frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T half_length = wire.length_m / T(2.0);

        const T distance_from_end = half_length - std::abs(position_m);
        const T numerator = std::sin(k * distance_from_end);
        const T denominator = std::sin(k * half_length);

        if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
            throw core::numerical_error("Denominator too small in current distribution calculation");
        }

        const complex_type i_zero(T(1.0), T(0.0));
        return i_zero * complex_type(numerator / denominator, T(0.0));
    }

    static complex_type calculate_center_fed_current(T position_m, const wire_geometry& wire) {
        return calculate_dipole_current(position_m, wire);
    }

    static complex_type calculate_off_center_fed_current(T position_m, const wire_geometry& wire) {
        core::check_range(wire.feed_position_m, -wire.length_m/T(2.0), wire.length_m/T(2.0), 
                         "Feed position");

        const T wavelength = constants::physical<T>::c / wire.frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T half_length = wire.length_m / T(2.0);

        if (position_m >= wire.feed_position_m) {

            const T distance_to_end = half_length - position_m;
            const T feed_to_end = half_length - wire.feed_position_m;
            const T numerator = std::sin(k * distance_to_end);
            const T denominator = std::sin(k * feed_to_end);

            if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
                return complex_type(T(0.0), T(0.0));
            }

            return complex_type(numerator / denominator, T(0.0));
        } else {

            const T distance_to_end = half_length + position_m;
            const T end_to_feed = half_length + wire.feed_position_m;
            const T numerator = std::sin(k * distance_to_end);
            const T denominator = std::sin(k * end_to_feed);

            if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
                return complex_type(T(0.0), T(0.0));
            }

            return complex_type(numerator / denominator, T(0.0));
        }
    }

    static complex_type calculate_monopole_current(T position_m, const wire_geometry& wire) {
        core::check_non_negative(position_m, "Monopole position");
        core::check_range(position_m, T(0.0), wire.length_m, "Position along monopole");

        const T wavelength = constants::physical<T>::c / wire.frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;

        const T distance_from_end = wire.length_m - position_m;
        const T numerator = std::sin(k * distance_from_end);
        const T denominator = std::sin(k * wire.length_m);

        if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
            throw core::numerical_error("Denominator too small in monopole current calculation");
        }

        return complex_type(numerator / denominator, T(0.0));
    }

    static complex_type calculate_loop_current(T angle_rad, const wire_geometry& wire) {
        core::check_range(angle_rad, T(0.0), T(2.0) * constants::mathematical<T>::pi, "Angle around loop");

        const T wavelength = constants::physical<T>::c / wire.frequency_hz;
        const T circumference = wire.length_m;
        const T radius = circumference / (T(2.0) * constants::mathematical<T>::pi);

        if (circumference < wavelength / T(10.0)) {
            return complex_type(T(1.0), T(0.0)); 
        }

        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T phase_variation = k * radius * std::cos(angle_rad);

        return complex_type(std::cos(phase_variation), std::sin(phase_variation));
    }

    static complex_type apply_wire_radius_correction(const complex_type& thin_wire_current,
                                                   const wire_geometry& wire) {
        core::check_positive(wire.radius_m, "Wire radius");

        if (wire.radius_m >= wire.length_m / T(100.0)) {

            const T wavelength = constants::physical<T>::c / wire.frequency_hz;
            const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
            const T ka = k * wire.radius_m;

            const T correction_factor = T(1.0) + T(0.1) * ka * ka;
            return thin_wire_current * correction_factor;
        }

        return thin_wire_current; 
    }

    static std::vector<current_sample> generate_current_samples(const wire_geometry& wire,
                                                              size_t num_samples = 100,
                                                              bool center_fed = true) {
        core::check_positive(num_samples, "Number of samples");

        std::vector<current_sample> samples;
        samples.reserve(num_samples);

        const T half_length = wire.length_m / T(2.0);
        const T step = wire.length_m / T(num_samples - 1);

        for (size_t i = 0; i < num_samples; ++i) {
            current_sample sample;
            sample.position_m = -half_length + i * step;

            if (center_fed) {
                sample.current = calculate_center_fed_current(sample.position_m, wire);
            } else {
                sample.current = calculate_off_center_fed_current(sample.position_m, wire);
            }

            sample.current = apply_wire_radius_correction(sample.current, wire);

            sample.magnitude = std::abs(sample.current);
            sample.phase_rad = std::arg(sample.current);

            samples.push_back(sample);
        }

        return samples;
    }

    static T calculate_rms_current(const std::vector<current_sample>& samples) {
        if (samples.empty()) {
            throw core::invalid_argument_error("Sample vector cannot be empty");
        }

        T sum_squared = T(0.0);
        for (const auto& sample : samples) {
            const T magnitude = std::abs(sample.current);
            sum_squared += magnitude * magnitude;
        }

        return std::sqrt(sum_squared / T(samples.size()));
    }

    static T calculate_peak_current(const std::vector<current_sample>& samples) {
        if (samples.empty()) {
            throw core::invalid_argument_error("Sample vector cannot be empty");
        }

        T max_current = T(0.0);
        for (const auto& sample : samples) {
            const T magnitude = std::abs(sample.current);
            max_current = std::max(max_current, magnitude);
        }

        return max_current;
    }

    static std::vector<T> find_current_nulls(const std::vector<current_sample>& samples, 
                                           T threshold = T(0.01)) {
        std::vector<T> nulls;

        for (size_t i = 1; i < samples.size() - 1; ++i) {
            const T prev_mag = std::abs(samples[i-1].current);
            const T curr_mag = std::abs(samples[i].current);
            const T next_mag = std::abs(samples[i+1].current);

            if (curr_mag < threshold && curr_mag < prev_mag && curr_mag < next_mag) {
                nulls.push_back(samples[i].position_m);
            }
        }

        return nulls;
    }

    static complex_type calculate_loaded_wire_current(T position_m, const wire_geometry& wire,
                                                    const std::vector<T>& load_positions,
                                                    const std::vector<complex_type>& load_impedances) {
        if (load_positions.size() != load_impedances.size()) {
            throw core::dimension_mismatch_error("Load positions and impedances must have same size");
        }

        complex_type current = calculate_center_fed_current(position_m, wire);

        const T wavelength = constants::physical<T>::c / wire.frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;

        for (size_t i = 0; i < load_positions.size(); ++i) {
            const T distance_to_load = std::abs(position_m - load_positions[i]);
            const T phase_shift = k * distance_to_load;

            const T load_magnitude = std::abs(load_impedances[i]);
            const T loading_factor = T(1.0) / (T(1.0) + load_magnitude / T(50.0)); 

            current *= complex_type(loading_factor * std::cos(phase_shift), 
                                  loading_factor * std::sin(phase_shift));
        }

        return current;
    }

    static void calculate_current_distribution_batch(const wire_geometry& base_wire,
                                                   const vector_type& frequencies,
                                                   const vector_type& positions,
                                                   std::vector<complex_vector_type>& current_matrix) {
        const size_t num_freq = frequencies.size();
        const size_t num_pos = positions.size();

        current_matrix.resize(num_freq);
        for (auto& row : current_matrix) {
            row.resize(num_pos);
        }

        for (size_t f = 0; f < num_freq; ++f) {
            wire_geometry wire = base_wire;
            wire.frequency_hz = frequencies[f];

            for (size_t p = 0; p < num_pos; ++p) {
                current_matrix[f][p] = calculate_center_fed_current(positions[p], wire);
            }
        }
    }

    static constexpr T typical_wire_radius_awg12() { return constants::antenna_design<T>::typical_wire_radius_awg12; }
    static constexpr T typical_wire_radius_awg14() { return constants::antenna_design<T>::typical_wire_radius_awg14; }
    static constexpr T typical_wire_radius_awg16() { return constants::antenna_design<T>::typical_wire_radius_awg16; }
    static constexpr T minimum_samples_per_wavelength() { return constants::antenna_design<T>::minimum_samples_per_wavelength; }
    static constexpr T maximum_useful_length_wavelengths() { return constants::antenna_design<T>::maximum_useful_length_wavelengths; }
};

using current_distribution_f = current_distribution<float>;
using current_distribution_d = current_distribution<double>;

} 
} 

#endif 