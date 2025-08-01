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

#ifndef RVL_ANTENNA_ANTENNA_GAIN_HPP
#define RVL_ANTENNA_ANTENNA_GAIN_HPP

#include "../core/constants.hpp"
#include "../core/units.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace antenna {

template<typename T>
class antenna_gain {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static constexpr T numeric_to_dbi(T gain_numeric) {
        core::check_positive(gain_numeric, "Numeric gain");
        return T(10.0) * std::log10(gain_numeric);
    }

    static constexpr T dbi_to_numeric(T gain_dbi) {
        return std::pow(T(10.0), gain_dbi / T(10.0));
    }

    static constexpr T dbi_to_dbd(T gain_dbi) {
        return gain_dbi - T(2.15);
    }

    static constexpr T dbd_to_dbi(T gain_dbd) {
        return gain_dbd + T(2.15);
    }

    static constexpr T halfwave_dipole_gain_dbi() {
        return T(2.15);
    }

    static constexpr T halfwave_dipole_gain_numeric() {
        return T(1.64);
    }

    static void numeric_to_dbi_batch(const vector_type& gain_numeric, vector_type& gain_dbi) {
        if (gain_numeric.size() != gain_dbi.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = gain_numeric.size();
        for (size_t i = 0; i < n; ++i) {
            if (gain_numeric[i] <= T(0)) {
                throw core::invalid_argument_error("All numeric gains must be positive");
            }
            gain_dbi[i] = T(10.0) * std::log10(gain_numeric[i]);
        }
    }

    static void dbi_to_numeric_batch(const vector_type& gain_dbi, vector_type& gain_numeric) {
        if (gain_dbi.size() != gain_numeric.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = gain_dbi.size();
        for (size_t i = 0; i < n; ++i) {
            gain_numeric[i] = std::pow(T(10.0), gain_dbi[i] / T(10.0));
        }
    }

    static T effective_area(T gain_numeric, T frequency) {
        core::check_positive(gain_numeric, "Gain");
        core::check_positive(frequency, "Frequency");

        const T wavelength = constants::physical<T>::c / frequency;
        const T wavelength_sq = wavelength * wavelength;
        return gain_numeric * wavelength_sq / (T(4.0) * constants::mathematical<T>::pi);
    }

    static T gain_from_effective_area(T effective_area, T frequency) {
        core::check_positive(effective_area, "Effective area");
        core::check_positive(frequency, "Frequency");

        const T wavelength = constants::physical<T>::c / frequency;
        const T wavelength_sq = wavelength * wavelength;
        return effective_area * T(4.0) * constants::mathematical<T>::pi / wavelength_sq;
    }
};

using antenna_gain_f = antenna_gain<float>;
using antenna_gain_d = antenna_gain<double>;

} 
} 

#endif 