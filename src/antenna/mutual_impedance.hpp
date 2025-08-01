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

#ifndef RVL_ANTENNA_MUTUAL_IMPEDANCE_HPP
#define RVL_ANTENNA_MUTUAL_IMPEDANCE_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <complex>
#include <cmath>

namespace rvl {
namespace antenna {

template<typename T>
class mutual_impedance {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    static complex_type calculate_simple(T distance, T frequency, T length1 = T(0.25), T length2 = T(0.25)) {
        core::check_positive(distance, "Distance");
        core::check_positive(frequency, "Frequency");
        core::check_positive(length1, "Antenna 1 length");
        core::check_positive(length2, "Antenna 2 length");

        const T wavelength = constants::physical<T>::c / frequency;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency;

        const T kr = k * distance;
        const T factor = constants::physical<T>::mu_0 * omega / (T(4.0) * constants::mathematical<T>::pi);

        const complex_type exp_term = std::exp(complex_type(T(0.0), -kr));
        const complex_type impedance_factor = exp_term / distance;

        const T geometry_factor = length1 * length2;

        return complex_type(T(0.0), factor) * impedance_factor * geometry_factor;
    }

    static complex_type calculate_parallel_dipoles(T separation, T frequency, T length = T(0.25)) {
        core::check_positive(separation, "Separation distance");
        core::check_positive(frequency, "Frequency");
        core::check_positive(length, "Dipole length");

        const T wavelength = constants::physical<T>::c / frequency;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T ks = k * separation;

        const complex_type j(T(0.0), T(1.0));

        if (separation < wavelength / T(10.0)) {
            const T real_part = T(30.0) * (T(2.0) * std::log(separation / (length * T(0.001))) - T(1.0));
            const T imag_part = T(30.0) * ks;
            return complex_type(real_part, imag_part);
        }

        const T sin_ks = std::sin(ks);
        const T cos_ks = std::cos(ks);
        const T cin_ks = (T(1.0) - cos_ks) / ks;
        const T si_ks = sin_ks / ks;

        const T real_part = T(30.0) * (T(2.0) * cin_ks - sin_ks);
        const T imag_part = T(30.0) * (T(2.0) * si_ks - cos_ks);

        return complex_type(real_part, imag_part);
    }

    static void calculate_array_mutual_impedances(const vector_type& x_positions,
                                                const vector_type& y_positions,
                                                T frequency,
                                                std::vector<std::vector<complex_type>>& z_matrix) {
        const size_t n = x_positions.size();
        if (y_positions.size() != n) {
            throw core::dimension_mismatch_error("Position vectors must have same size");
        }

        z_matrix.resize(n);
        for (size_t i = 0; i < n; ++i) {
            z_matrix[i].resize(n);
        }

        const T wavelength = constants::physical<T>::c / frequency;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency;
        const T factor = constants::physical<T>::mu_0 * omega / (T(4.0) * constants::mathematical<T>::pi);

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    z_matrix[i][j] = complex_type(constants::antenna<T>::dipole_input_impedance, constants::antenna<T>::dipole_reactance);
                } else {
                    const T dx = x_positions[i] - x_positions[j];
                    const T dy = y_positions[i] - y_positions[j];
                    const T distance = std::sqrt(dx * dx + dy * dy);

                    if (distance < std::numeric_limits<T>::epsilon()) {
                        z_matrix[i][j] = complex_type(constants::antenna<T>::dipole_input_impedance, constants::antenna<T>::dipole_reactance);
                    } else {
                        const T kr = k * distance;
                        const complex_type exp_term = std::exp(complex_type(T(0.0), -kr));
                        z_matrix[i][j] = complex_type(T(0.0), factor) * exp_term / distance;
                    }
                }
            }
        }
    }

    static T coupling_coefficient(const complex_type& z_mutual, const complex_type& z_self) {
        const T z_mut_mag = std::abs(z_mutual);
        const T z_self_real = std::real(z_self);

        if (z_self_real <= T(0)) {
            throw core::invalid_argument_error("Self impedance real part must be positive");
        }

        return z_mut_mag / z_self_real;
    }
};

using mutual_impedance_f = mutual_impedance<float>;
using mutual_impedance_d = mutual_impedance<double>;

} 
} 

#endif 