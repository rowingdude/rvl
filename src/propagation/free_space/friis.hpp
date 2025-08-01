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

#ifndef RVL_PROPAGATION_FREE_SPACE_FRIIS_HPP
#define RVL_PROPAGATION_FREE_SPACE_FRIIS_HPP

#include "../../core/constants.hpp"
#include "../../core/units.hpp"
#include "../../core/error.hpp"
#include "../../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace propagation {
namespace free_space {

template<typename T>
class friis {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static T received_power_watts(T transmit_power_w, T transmit_gain_numeric, 
                                T receive_gain_numeric, T wavelength_m, T distance_m) {
        core::check_positive(transmit_power_w, "Transmit power");
        core::check_positive(transmit_gain_numeric, "Transmit gain");
        core::check_positive(receive_gain_numeric, "Receive gain");
        core::check_positive(wavelength_m, "Wavelength");
        core::check_positive(distance_m, "Distance");

        const T wavelength_factor = wavelength_m / (T(4.0) * constants::mathematical<T>::pi * distance_m);
        const T path_gain = wavelength_factor * wavelength_factor;

        return transmit_power_w * transmit_gain_numeric * receive_gain_numeric * path_gain;
    }

    static T received_power_dbm(T transmit_power_dbm, T transmit_gain_dbi,
                              T receive_gain_dbi, T frequency_hz, T distance_m) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(distance_m, "Distance");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T path_loss_db = T(20.0) * std::log10(T(4.0) * constants::mathematical<T>::pi * distance_m / wavelength);

        return transmit_power_dbm + transmit_gain_dbi + receive_gain_dbi - path_loss_db;
    }

    static T effective_isotropic_radiated_power_dbm(T transmit_power_dbm, T transmit_gain_dbi) {
        return transmit_power_dbm + transmit_gain_dbi;
    }

    static T effective_isotropic_radiated_power_watts(T transmit_power_w, T transmit_gain_numeric) {
        core::check_positive(transmit_power_w, "Transmit power");
        core::check_positive(transmit_gain_numeric, "Transmit gain");

        return transmit_power_w * transmit_gain_numeric;
    }

    static T range_for_received_power(T transmit_power_w, T transmit_gain_numeric,
                                    T receive_gain_numeric, T wavelength_m, T required_power_w) {
        core::check_positive(transmit_power_w, "Transmit power");
        core::check_positive(transmit_gain_numeric, "Transmit gain");
        core::check_positive(receive_gain_numeric, "Receive gain");
        core::check_positive(wavelength_m, "Wavelength");
        core::check_positive(required_power_w, "Required power");

        const T gain_product = transmit_gain_numeric * receive_gain_numeric;
        const T power_ratio = transmit_power_w / required_power_w;
        const T range_factor = wavelength_m / (T(4.0) * constants::mathematical<T>::pi);

        return range_factor * std::sqrt(gain_product * power_ratio);
    }

    static T range_for_received_power_db(T transmit_power_dbm, T transmit_gain_dbi,
                                       T receive_gain_dbi, T frequency_hz, T required_power_dbm) {
        core::check_positive(frequency_hz, "Frequency");

        const T link_budget_db = transmit_power_dbm + transmit_gain_dbi + receive_gain_dbi - required_power_dbm;
        const T wavelength = constants::physical<T>::c / frequency_hz;

        const T log_lambda_4pi = std::log10(wavelength / (T(4.0) * constants::mathematical<T>::pi));
        const T log_range = (link_budget_db / T(20.0)) + log_lambda_4pi;

        return std::pow(T(10.0), log_range);
    }

    static T power_density_at_range(T eirp_watts, T distance_m) {
        core::check_positive(eirp_watts, "EIRP");
        core::check_positive(distance_m, "Distance");

        const T sphere_area = T(4.0) * constants::mathematical<T>::pi * distance_m * distance_m;
        return eirp_watts / sphere_area;
    }

    static T field_strength_at_range(T eirp_watts, T distance_m) {
        core::check_positive(eirp_watts, "EIRP");
        core::check_positive(distance_m, "Distance");

        const T power_density = power_density_at_range(eirp_watts, distance_m);
        const T impedance = constants::physical<T>::eta_0;

        return std::sqrt(power_density * impedance);
    }

    static void received_power_batch(const vector_type& transmit_powers,
                                   const vector_type& transmit_gains,
                                   const vector_type& receive_gains,
                                   const vector_type& frequencies,
                                   const vector_type& distances,
                                   vector_type& received_powers) {
        const size_t n = transmit_powers.size();
        if (transmit_gains.size() != n || receive_gains.size() != n || 
            frequencies.size() != n || distances.size() != n || received_powers.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        const T c = constants::physical<T>::c;
        const T four_pi = T(4.0) * constants::mathematical<T>::pi;

        for (size_t i = 0; i < n; ++i) {
            if (transmit_powers[i] <= T(0) || transmit_gains[i] <= T(0) || 
                receive_gains[i] <= T(0) || frequencies[i] <= T(0) || distances[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            const T wavelength = c / frequencies[i];
            const T wavelength_factor = wavelength / (four_pi * distances[i]);
            const T path_gain = wavelength_factor * wavelength_factor;

            received_powers[i] = transmit_powers[i] * transmit_gains[i] * receive_gains[i] * path_gain;
        }
    }

    static void received_power_dbm_batch(const vector_type& transmit_powers_dbm,
                                       const vector_type& transmit_gains_dbi,
                                       const vector_type& receive_gains_dbi,
                                       const vector_type& frequencies,
                                       const vector_type& distances,
                                       vector_type& received_powers_dbm) {
        const size_t n = transmit_powers_dbm.size();
        if (transmit_gains_dbi.size() != n || receive_gains_dbi.size() != n || 
            frequencies.size() != n || distances.size() != n || received_powers_dbm.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        const T c = constants::physical<T>::c;
        const T four_pi = T(4.0) * constants::mathematical<T>::pi;

        for (size_t i = 0; i < n; ++i) {
            if (frequencies[i] <= T(0) || distances[i] <= T(0)) {
                throw core::invalid_argument_error("Frequencies and distances must be positive");
            }

            const T wavelength = c / frequencies[i];
            const T path_loss_db = T(20.0) * std::log10(four_pi * distances[i] / wavelength);

            received_powers_dbm[i] = transmit_powers_dbm[i] + transmit_gains_dbi[i] + 
                                   receive_gains_dbi[i] - path_loss_db;
        }
    }
};

using friis_f = friis<float>;
using friis_d = friis<double>;

} 
} 
} 

#endif 