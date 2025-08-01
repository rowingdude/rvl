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

#ifndef RVL_FEEDLINE_VSWR_HPP
#define RVL_FEEDLINE_VSWR_HPP

#include "../core/constants.hpp"
#include "../core/units.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <complex>
#include <cmath>

namespace rvl {
namespace feedline {

template<typename T>
class vswr {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    static T from_reflection_coefficient(const complex_type& gamma) {
        const T gamma_mag = std::abs(gamma);

        if (gamma_mag >= T(1.0)) {
            return std::numeric_limits<T>::infinity();
        }

        return (T(1.0) + gamma_mag) / (T(1.0) - gamma_mag);
    }

    static T from_reflection_coefficient(T gamma_magnitude) {
        core::check_range(gamma_magnitude, T(0.0), T(1.0), "Reflection coefficient magnitude");

        if (gamma_magnitude >= T(1.0)) {
            return std::numeric_limits<T>::infinity();
        }

        return (T(1.0) + gamma_magnitude) / (T(1.0) - gamma_magnitude);
    }

    static T from_impedances(const complex_type& z_load, T z_characteristic = T(50.0)) {
        core::check_positive(z_characteristic, "Characteristic impedance");

        const complex_type gamma = (z_load - complex_type(z_characteristic, T(0.0))) /
                                 (z_load + complex_type(z_characteristic, T(0.0)));

        return from_reflection_coefficient(gamma);
    }

    static complex_type reflection_coefficient_from_vswr(T vswr_value, T phase_radians = T(0.0)) {
        core::check_range(vswr_value, T(1.0), std::numeric_limits<T>::max(), "VSWR");

        const T gamma_mag = (vswr_value - T(1.0)) / (vswr_value + T(1.0));

        return gamma_mag * std::exp(complex_type(T(0.0), phase_radians));
    }

    static T return_loss_db(T vswr_value) {
        core::check_range(vswr_value, T(1.0), std::numeric_limits<T>::max(), "VSWR");

        const T gamma = (vswr_value - T(1.0)) / (vswr_value + T(1.0));
        return -T(20.0) * std::log10(gamma);
    }

    static T from_return_loss_db(T return_loss_db) {
        core::check_positive(return_loss_db, "Return loss");

        const T gamma = std::pow(T(10.0), -return_loss_db / T(20.0));
        return from_reflection_coefficient(gamma);
    }

    static T mismatch_loss_db(T vswr_value) {
        core::check_range(vswr_value, T(1.0), std::numeric_limits<T>::max(), "VSWR");

        const T gamma = (vswr_value - T(1.0)) / (vswr_value + T(1.0));
        const T mismatch_factor = T(1.0) - gamma * gamma;

        return -T(10.0) * std::log10(mismatch_factor);
    }

    static T power_delivered_fraction(T vswr_value) {
        core::check_range(vswr_value, T(1.0), std::numeric_limits<T>::max(), "VSWR");

        const T gamma = (vswr_value - T(1.0)) / (vswr_value + T(1.0));
        return T(1.0) - gamma * gamma;
    }

    static T power_reflected_fraction(T vswr_value) {
        core::check_range(vswr_value, T(1.0), std::numeric_limits<T>::max(), "VSWR");

        const T gamma = (vswr_value - T(1.0)) / (vswr_value + T(1.0));
        return gamma * gamma;
    }

    static void from_reflection_coefficient_batch(const complex_vector_type& gamma_values,
                                                vector_type& vswr_values) {
        if (gamma_values.size() != vswr_values.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = gamma_values.size();
        for (size_t i = 0; i < n; ++i) {
            vswr_values[i] = from_reflection_coefficient(gamma_values[i]);
        }
    }

    static void from_impedances_batch(const complex_vector_type& z_loads,
                                    vector_type& vswr_values,
                                    T z_characteristic = T(50.0)) {
        if (z_loads.size() != vswr_values.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = z_loads.size();
        const complex_type z_char(z_characteristic, T(0.0));

        for (size_t i = 0; i < n; ++i) {
            const complex_type gamma = (z_loads[i] - z_char) / (z_loads[i] + z_char);
            vswr_values[i] = from_reflection_coefficient(gamma);
        }
    }

    static void return_loss_batch(const vector_type& vswr_values, vector_type& return_loss_values) {
        if (vswr_values.size() != return_loss_values.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = vswr_values.size();
        for (size_t i = 0; i < n; ++i) {
            if (vswr_values[i] < T(1.0)) {
                throw core::out_of_range_error("VSWR values must be >= 1.0");
            }

            return_loss_values[i] = return_loss_db(vswr_values[i]);
        }
    }

    static void mismatch_loss_batch(const vector_type& vswr_values, vector_type& mismatch_loss_values) {
        if (vswr_values.size() != mismatch_loss_values.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = vswr_values.size();
        for (size_t i = 0; i < n; ++i) {
            if (vswr_values[i] < T(1.0)) {
                throw core::out_of_range_error("VSWR values must be >= 1.0");
            }

            mismatch_loss_values[i] = mismatch_loss_db(vswr_values[i]);
        }
    }

    static constexpr T perfect_match() { return constants::vswr_thresholds<T>::perfect_match; }
    static constexpr T acceptable_amateur() { return constants::vswr_thresholds<T>::acceptable_amateur; }
    static constexpr T acceptable_commercial() { return constants::vswr_thresholds<T>::acceptable_commercial; }

    static T voltage_nodes_distance(T vswr_value, T wavelength) {
        core::check_range(vswr_value, T(1.0), std::numeric_limits<T>::max(), "VSWR");
        core::check_positive(wavelength, "Wavelength");

        return wavelength / T(2.0);
    }
};

using vswr_f = vswr<float>;
using vswr_d = vswr<double>;

} 
} 

#endif 