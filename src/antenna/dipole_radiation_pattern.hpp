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

#ifndef RVL_ANTENNA_DIPOLE_RADIATION_PATTERN_HPP
#define RVL_ANTENNA_DIPOLE_RADIATION_PATTERN_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace antenna {

template<typename T>
class dipole_radiation_pattern {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static T electric_field_pattern(T theta) {
        core::check_finite(theta, "Theta angle");

        const T sin_theta = std::sin(theta);
        if (std::abs(sin_theta) < std::numeric_limits<T>::epsilon()) {
            return T(0.0);
        }

        const T cos_theta = std::cos(theta);
        const T numerator = std::cos(constants::mathematical<T>::half_pi * cos_theta);

        return numerator / sin_theta;
    }

    static T power_pattern(T theta) {
        const T e_field = electric_field_pattern(theta);
        return e_field * e_field;
    }

    static T power_pattern_db(T theta) {
        const T power = power_pattern(theta);
        if (power <= T(0)) {
            return T(-100.0);
        }
        return T(10.0) * std::log10(power);
    }

    static void electric_field_pattern_batch(const vector_type& theta_angles, vector_type& patterns) {
        if (theta_angles.size() != patterns.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = theta_angles.size();
        const T half_pi = constants::mathematical<T>::half_pi;

        for (size_t i = 0; i < n; ++i) {
            const T theta = theta_angles[i];
            const T sin_theta = std::sin(theta);

            if (std::abs(sin_theta) < std::numeric_limits<T>::epsilon()) {
                patterns[i] = T(0.0);
            } else {
                const T cos_theta = std::cos(theta);
                const T numerator = std::cos(half_pi * cos_theta);
                patterns[i] = numerator / sin_theta;
            }
        }
    }

    static void power_pattern_batch(const vector_type& theta_angles, vector_type& patterns) {
        electric_field_pattern_batch(theta_angles, patterns);

        const size_t n = patterns.size();
        for (size_t i = 0; i < n; ++i) {
            patterns[i] = patterns[i] * patterns[i];
        }
    }

    static void power_pattern_db_batch(const vector_type& theta_angles, vector_type& patterns_db) {
        power_pattern_batch(theta_angles, patterns_db);

        const size_t n = patterns_db.size();
        for (size_t i = 0; i < n; ++i) {
            if (patterns_db[i] <= T(0)) {
                patterns_db[i] = T(-100.0);
            } else {
                patterns_db[i] = T(10.0) * std::log10(patterns_db[i]);
            }
        }
    }

    static T directivity() {
        return T(1.64);
    }

    static T directivity_db() {
        return T(2.15);
    }

    static T beamwidth_e_plane_deg() {
        return T(78.0);
    }

    static T beamwidth_h_plane_deg() {
        return T(360.0);
    }
};

using dipole_radiation_pattern_f = dipole_radiation_pattern<float>;
using dipole_radiation_pattern_d = dipole_radiation_pattern<double>;

} 
} 

#endif 