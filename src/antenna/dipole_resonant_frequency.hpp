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

#ifndef RVL_ANTENNA_DIPOLE_RESONANT_FREQUENCY_HPP
#define RVL_ANTENNA_DIPOLE_RESONANT_FREQUENCY_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace antenna {

template<typename T>
class dipole_resonant_frequency {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static constexpr T calculate_frequency(T length) {
        core::check_positive(length, "Antenna length");
        return constants::physical<T>::c / (T(2.0) * length);
    }

    static constexpr T calculate_length(T frequency) {
        core::check_positive(frequency, "Frequency");
        return constants::physical<T>::c / (T(2.0) * frequency);
    }

    static void calculate_frequency_batch(const vector_type& lengths, vector_type& frequencies) {
        if (lengths.size() != frequencies.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = lengths.size();
        const T c_over_2 = constants::physical<T>::c / T(2.0);

        for (size_t i = 0; i < n; ++i) {
            if (lengths[i] <= T(0)) {
                throw core::invalid_argument_error("All antenna lengths must be positive");
            }
            frequencies[i] = c_over_2 / lengths[i];
        }
    }

    static void calculate_length_batch(const vector_type& frequencies, vector_type& lengths) {
        if (frequencies.size() != lengths.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = frequencies.size();
        const T c_over_2 = constants::physical<T>::c / T(2.0);

        for (size_t i = 0; i < n; ++i) {
            if (frequencies[i] <= T(0)) {
                throw core::invalid_argument_error("All frequencies must be positive");
            }
            lengths[i] = c_over_2 / frequencies[i];
        }
    }
};

using dipole_resonant_frequency_f = dipole_resonant_frequency<float>;
using dipole_resonant_frequency_d = dipole_resonant_frequency<double>;

} 
} 

#endif 