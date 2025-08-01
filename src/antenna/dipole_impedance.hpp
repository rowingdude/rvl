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

#ifndef RVL_ANTENNA_DIPOLE_IMPEDANCE_HPP
#define RVL_ANTENNA_DIPOLE_IMPEDANCE_HPP

#include "../core/constants.hpp"
#include "../core/memory_allocator.hpp"
#include <complex>

namespace rvl {
namespace antenna {

template<typename T>
class dipole_impedance {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<complex_type>;

    static constexpr T resistance_halfwave = constants::antenna<T>::dipole_input_impedance;
    static constexpr T reactance_halfwave = constants::antenna<T>::dipole_reactance;

    static constexpr complex_type ideal_halfwave_impedance() {
        return complex_type(resistance_halfwave, reactance_halfwave);
    }

    static complex_type calculate_impedance(T length, T frequency, T wire_radius = T(0.001)) {
        core::check_positive(length, "Antenna length");
        core::check_positive(frequency, "Frequency");
        core::check_positive(wire_radius, "Wire radius");

        const T wavelength = constants::physical<T>::c / frequency;
        const T electrical_length = length / wavelength;

        if (std::abs(electrical_length - T(0.5)) < T(0.01)) {
            return ideal_halfwave_impedance();
        }

        const T beta_l = T(2.0) * constants::mathematical<T>::pi * electrical_length;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T ln_term = std::log(length / wire_radius);

        const T resistance = T(30.0) * (T(1.0) - std::cos(beta_l)) / std::sin(beta_l);
        const T reactance = T(30.0) * (T(2.0) * ln_term - T(1.0)) * 
                           (std::sin(beta_l) - beta_l * std::cos(beta_l)) / std::sin(beta_l);

        return complex_type(resistance, reactance);
    }

    static void calculate_impedance_batch(const core::memory::simd_vector<T>& lengths,
                                        const core::memory::simd_vector<T>& frequencies,
                                        vector_type& impedances,
                                        T wire_radius = T(0.001)) {
        const size_t n = lengths.size();
        if (frequencies.size() != n || impedances.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            impedances[i] = calculate_impedance(lengths[i], frequencies[i], wire_radius);
        }
    }
};

using dipole_impedance_f = dipole_impedance<float>;
using dipole_impedance_d = dipole_impedance<double>;

} 
} 

#endif 