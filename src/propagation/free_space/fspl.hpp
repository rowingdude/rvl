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

#ifndef RVL_PROPAGATION_FREE_SPACE_FSPL_HPP
#define RVL_PROPAGATION_FREE_SPACE_FSPL_HPP

#include "../../core/constants.hpp"
#include "../../core/units.hpp"
#include "../../core/error.hpp"
#include "../../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace propagation {
namespace free_space {

template<typename T>
class fspl {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static T calculate_db(T distance_m, T frequency_hz) {
        core::check_positive(distance_m, "Distance");
        core::check_positive(frequency_hz, "Frequency");

        const T log_d = std::log10(distance_m);
        const T log_f = std::log10(frequency_hz);
        const T log_c = std::log10(constants::physical<T>::c);
        const T log_4pi = std::log10(T(4.0) * constants::mathematical<T>::pi);

        return T(20.0) * log_d + T(20.0) * log_f + T(20.0) * (log_4pi - log_c);
    }

    static T calculate_linear(T distance_m, T frequency_hz) {
        core::check_positive(distance_m, "Distance");
        core::check_positive(frequency_hz, "Frequency");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T factor = T(4.0) * constants::mathematical<T>::pi * distance_m / wavelength;

        return factor * factor;
    }

    static T calculate_from_wavelength_db(T distance_m, T wavelength_m) {
        core::check_positive(distance_m, "Distance");
        core::check_positive(wavelength_m, "Wavelength");

        const T ratio = distance_m / wavelength_m;
        const T log_ratio = std::log10(ratio);
        const T log_4pi = std::log10(T(4.0) * constants::mathematical<T>::pi);

        return T(20.0) * log_ratio + T(20.0) * log_4pi;
    }

    static T distance_from_fspl(T fspl_db, T frequency_hz) {
        core::check_positive(frequency_hz, "Frequency");

        const T log_f = std::log10(frequency_hz);
        const T log_c = std::log10(constants::physical<T>::c);
        const T log_4pi = std::log10(T(4.0) * constants::mathematical<T>::pi);

        const T log_d = (fspl_db / T(20.0)) - log_f - (log_4pi - log_c);

        return std::pow(T(10.0), log_d);
    }

    static T frequency_from_fspl(T fspl_db, T distance_m) {
        core::check_positive(distance_m, "Distance");

        const T log_d = std::log10(distance_m);
        const T log_c = std::log10(constants::physical<T>::c);
        const T log_4pi = std::log10(T(4.0) * constants::mathematical<T>::pi);

        const T log_f = (fspl_db / T(20.0)) - log_d - (log_4pi - log_c);

        return std::pow(T(10.0), log_f);
    }

    static void calculate_db_batch(const vector_type& distances, 
                                 const vector_type& frequencies,
                                 vector_type& fspl_db) {
        const size_t n = distances.size();
        if (frequencies.size() != n || fspl_db.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        const T log_c = std::log10(constants::physical<T>::c);
        const T log_4pi = std::log10(T(4.0) * constants::mathematical<T>::pi);
        const T constant_term = T(20.0) * (log_4pi - log_c);

        for (size_t i = 0; i < n; ++i) {
            if (distances[i] <= T(0) || frequencies[i] <= T(0)) {
                throw core::invalid_argument_error("All distances and frequencies must be positive");
            }

            const T log_d = std::log10(distances[i]);
            const T log_f = std::log10(frequencies[i]);
            fspl_db[i] = T(20.0) * log_d + T(20.0) * log_f + constant_term;
        }
    }

    static void calculate_linear_batch(const vector_type& distances,
                                     const vector_type& frequencies,
                                     vector_type& fspl_linear) {
        const size_t n = distances.size();
        if (frequencies.size() != n || fspl_linear.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        const T c = constants::physical<T>::c;
        const T four_pi = T(4.0) * constants::mathematical<T>::pi;

        for (size_t i = 0; i < n; ++i) {
            if (distances[i] <= T(0) || frequencies[i] <= T(0)) {
                throw core::invalid_argument_error("All distances and frequencies must be positive");
            }

            const T wavelength = c / frequencies[i];
            const T factor = four_pi * distances[i] / wavelength;
            fspl_linear[i] = factor * factor;
        }
    }

    static T link_margin_db(T transmit_power_dbm, T transmit_gain_dbi, T receive_gain_dbi,
                           T distance_m, T frequency_hz, T receiver_sensitivity_dbm) {
        const T fspl_db = calculate_db(distance_m, frequency_hz);
        const T received_power_dbm = transmit_power_dbm + transmit_gain_dbi + receive_gain_dbi - fspl_db;

        return received_power_dbm - receiver_sensitivity_dbm;
    }

    static constexpr T fspl_1km_1ghz_db() {
        return T(92.45);
    }

    static constexpr T fspl_1km_1mhz_db() {
        return T(32.45);
    }
};

using fspl_f = fspl<float>;
using fspl_d = fspl<double>;

} 
} 
} 

#endif 