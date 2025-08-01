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

#ifndef RVL_PROPAGATION_FREE_SPACE_FRIIS_POLARIZATION_HPP
#define RVL_PROPAGATION_FREE_SPACE_FRIIS_POLARIZATION_HPP

#include "../../core/constants.hpp"
#include "../../core/error.hpp"
#include "../../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace propagation {
namespace free_space {

template<typename T>
class friis_polarization {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    enum class polarization_type {
        LINEAR_VERTICAL,
        LINEAR_HORIZONTAL,
        CIRCULAR_RIGHT,
        CIRCULAR_LEFT,
        ELLIPTICAL
    };

    static T received_power_watts(T transmit_power_w, T transmit_gain_numeric,
                                T receive_gain_numeric, T wavelength_m, T distance_m,
                                T polarization_mismatch_angle_rad) {
        core::check_positive(transmit_power_w, "Transmit power");
        core::check_positive(transmit_gain_numeric, "Transmit gain");
        core::check_positive(receive_gain_numeric, "Receive gain");
        core::check_positive(wavelength_m, "Wavelength");
        core::check_positive(distance_m, "Distance");

        const T wavelength_factor = wavelength_m / (T(4.0) * constants::mathematical<T>::pi * distance_m);
        const T path_gain = wavelength_factor * wavelength_factor;
        const T polarization_factor = std::cos(polarization_mismatch_angle_rad);
        const T polarization_loss = polarization_factor * polarization_factor;

        return transmit_power_w * transmit_gain_numeric * receive_gain_numeric * 
               path_gain * polarization_loss;
    }

    static T received_power_dbm(T transmit_power_dbm, T transmit_gain_dbi,
                              T receive_gain_dbi, T frequency_hz, T distance_m,
                              T polarization_mismatch_angle_rad) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(distance_m, "Distance");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T path_loss_db = T(20.0) * std::log10(T(4.0) * constants::mathematical<T>::pi * distance_m / wavelength);
        const T polarization_loss_db = polarization_mismatch_loss_db(polarization_mismatch_angle_rad);

        return transmit_power_dbm + transmit_gain_dbi + receive_gain_dbi - path_loss_db - polarization_loss_db;
    }

    static T polarization_mismatch_loss_db(T mismatch_angle_rad) {
        const T cos_angle = std::cos(mismatch_angle_rad);
        const T loss_linear = cos_angle * cos_angle;

        if (loss_linear <= T(0)) {
            return T(100.0);
        }

        return -T(10.0) * std::log10(loss_linear);
    }

    static T polarization_mismatch_angle(polarization_type tx_pol, polarization_type rx_pol) {
        if (tx_pol == rx_pol) {
            return T(0.0);
        }

        if ((tx_pol == polarization_type::LINEAR_VERTICAL && rx_pol == polarization_type::LINEAR_HORIZONTAL) ||
            (tx_pol == polarization_type::LINEAR_HORIZONTAL && rx_pol == polarization_type::LINEAR_VERTICAL)) {
            return constants::mathematical<T>::half_pi;
        }

        if ((tx_pol == polarization_type::CIRCULAR_RIGHT && rx_pol == polarization_type::CIRCULAR_LEFT) ||
            (tx_pol == polarization_type::CIRCULAR_LEFT && rx_pol == polarization_type::CIRCULAR_RIGHT)) {
            return constants::mathematical<T>::pi;
        }

        if (tx_pol == polarization_type::LINEAR_VERTICAL || tx_pol == polarization_type::LINEAR_HORIZONTAL) {
            if (rx_pol == polarization_type::CIRCULAR_LEFT || rx_pol == polarization_type::CIRCULAR_RIGHT) {
                return constants::mathematical<T>::pi / T(4.0);
            }
        }

        if (rx_pol == polarization_type::LINEAR_VERTICAL || rx_pol == polarization_type::LINEAR_HORIZONTAL) {
            if (tx_pol == polarization_type::CIRCULAR_LEFT || tx_pol == polarization_type::CIRCULAR_RIGHT) {
                return constants::mathematical<T>::pi / T(4.0);
            }
        }

        return constants::mathematical<T>::pi / T(4.0);
    }

    static T cross_polarization_discrimination_db(T mismatch_angle_rad) {
        return polarization_mismatch_loss_db(mismatch_angle_rad);
    }

    static T axial_ratio_to_ellipticity_angle(T axial_ratio_db) {
        const T axial_ratio_linear = std::pow(T(10.0), axial_ratio_db / T(20.0));
        return std::atan(T(1.0) / axial_ratio_linear);
    }

    static T ellipticity_angle_to_axial_ratio_db(T ellipticity_angle_rad) {
        const T tan_angle = std::tan(ellipticity_angle_rad);
        if (tan_angle <= T(0)) {
            return T(100.0);
        }
        const T axial_ratio_linear = T(1.0) / tan_angle;
        return T(20.0) * std::log10(axial_ratio_linear);
    }

    static void received_power_polarization_batch(const vector_type& transmit_powers,
                                                const vector_type& transmit_gains,
                                                const vector_type& receive_gains,
                                                const vector_type& frequencies,
                                                const vector_type& distances,
                                                const vector_type& mismatch_angles,
                                                vector_type& received_powers) {
        const size_t n = transmit_powers.size();
        if (transmit_gains.size() != n || receive_gains.size() != n || frequencies.size() != n ||
            distances.size() != n || mismatch_angles.size() != n || received_powers.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        const T c = constants::physical<T>::c;
        const T four_pi = T(4.0) * constants::mathematical<T>::pi;

        for (size_t i = 0; i < n; ++i) {
            if (transmit_powers[i] <= T(0) || transmit_gains[i] <= T(0) || receive_gains[i] <= T(0) ||
                frequencies[i] <= T(0) || distances[i] <= T(0)) {
                throw core::invalid_argument_error("All power, gain, frequency, and distance values must be positive");
            }

            const T wavelength = c / frequencies[i];
            const T wavelength_factor = wavelength / (four_pi * distances[i]);
            const T path_gain = wavelength_factor * wavelength_factor;
            const T polarization_factor = std::cos(mismatch_angles[i]);
            const T polarization_loss = polarization_factor * polarization_factor;

            received_powers[i] = transmit_powers[i] * transmit_gains[i] * receive_gains[i] *
                               path_gain * polarization_loss;
        }
    }

    static void polarization_loss_batch(const vector_type& mismatch_angles, vector_type& losses_db) {
        if (mismatch_angles.size() != losses_db.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = mismatch_angles.size();
        for (size_t i = 0; i < n; ++i) {
            losses_db[i] = polarization_mismatch_loss_db(mismatch_angles[i]);
        }
    }

    static constexpr T perfect_match_loss_db() { return T(0.0); }
    static constexpr T orthogonal_linear_loss_db() { return T(100.0); }
    static constexpr T linear_to_circular_loss_db() { return T(3.0); }
    static constexpr T opposite_circular_loss_db() { return T(100.0); }

    static constexpr T typical_multipath_depolarization_db() { return T(3.0); }
    static constexpr T rain_depolarization_db_per_degree() { return T(0.1); }
};

using friis_polarization_f = friis_polarization<float>;
using friis_polarization_d = friis_polarization<double>;

} 
} 
} 

#endif 