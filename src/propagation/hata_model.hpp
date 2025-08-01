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

#ifndef RVL_PROPAGATION_HATA_MODEL_HPP
#define RVL_PROPAGATION_HATA_MODEL_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace propagation {

template<typename T>
class hata_model {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    enum class environment_type {
        URBAN_DENSE,
        URBAN,
        SUBURBAN,
        RURAL_OPEN,
        RURAL_QUASI_OPEN
    };

    enum class antenna_type {
        SMALL_CITY,      
        LARGE_CITY       
    };

    static T receiver_height_correction_small_city(T receiver_height_m, T frequency_mhz) {
        core::check_positive(receiver_height_m, "Receiver height");
        core::check_positive(frequency_mhz, "Frequency");

        const T freq_term = T(1.1) * std::log10(frequency_mhz) - T(0.7);
        const T height_term = T(1.56) * std::log10(receiver_height_m);

        return freq_term * receiver_height_m - height_term - T(0.8);
    }

    static T receiver_height_correction_large_city(T receiver_height_m, T frequency_mhz) {
        core::check_positive(receiver_height_m, "Receiver height");
        core::check_positive(frequency_mhz, "Frequency");

        T correction;
        if (frequency_mhz <= T(300.0)) {
            const T log_freq = std::log10(frequency_mhz);
            correction = T(8.29) * std::pow(log_freq, T(2.0)) - T(1.1);
        } else {
            correction = T(3.2) * std::pow(std::log10(T(11.75) * receiver_height_m), T(2.0)) - T(4.97);
        }

        return correction;
    }

    static T path_loss_urban_db(T frequency_mhz, T transmitter_height_m, T receiver_height_m,
                               T distance_km, antenna_type ant_type = antenna_type::SMALL_CITY) {
        core::check_range(frequency_mhz, T(150.0), T(1500.0), "Frequency (MHz)");
        core::check_range(transmitter_height_m, T(30.0), T(200.0), "Transmitter height");
        core::check_range(receiver_height_m, T(1.0), T(10.0), "Receiver height");
        core::check_range(distance_km, T(1.0), T(20.0), "Distance (km)");

        const T log_freq = std::log10(frequency_mhz);
        const T log_ht = std::log10(transmitter_height_m);
        const T log_d = std::log10(distance_km);

        T a_hr;
        if (ant_type == antenna_type::SMALL_CITY) {
            a_hr = receiver_height_correction_small_city(receiver_height_m, frequency_mhz);
        } else {
            a_hr = receiver_height_correction_large_city(receiver_height_m, frequency_mhz);
        }

        const T path_loss = T(69.55) + T(26.16) * log_freq - T(13.82) * log_ht - a_hr +
                           (T(44.9) - T(6.55) * log_ht) * log_d;

        return path_loss;
    }

    static T path_loss_suburban_db(T frequency_mhz, T transmitter_height_m, T receiver_height_m,
                                  T distance_km, antenna_type ant_type = antenna_type::SMALL_CITY) {
        const T urban_loss = path_loss_urban_db(frequency_mhz, transmitter_height_m, receiver_height_m,
                                               distance_km, ant_type);

        const T log_freq = std::log10(frequency_mhz);
        const T suburban_correction = T(2.0) * std::pow(log_freq / T(28.0), T(2.0)) + T(5.4);

        return urban_loss - suburban_correction;
    }

    static T path_loss_rural_open_db(T frequency_mhz, T transmitter_height_m, T receiver_height_m,
                                    T distance_km, antenna_type ant_type = antenna_type::SMALL_CITY) {
        const T urban_loss = path_loss_urban_db(frequency_mhz, transmitter_height_m, receiver_height_m,
                                               distance_km, ant_type);

        const T log_freq = std::log10(frequency_mhz);
        const T rural_correction = T(4.78) * std::pow(log_freq, T(2.0)) - T(18.33) * log_freq + T(40.94);

        return urban_loss - rural_correction;
    }

    static T path_loss_by_environment(T frequency_mhz, T transmitter_height_m, T receiver_height_m,
                                     T distance_km, environment_type env,
                                     antenna_type ant_type = antenna_type::SMALL_CITY) {
        switch (env) {
            case environment_type::URBAN_DENSE:
                return path_loss_urban_db(frequency_mhz, transmitter_height_m, receiver_height_m,
                                        distance_km, antenna_type::LARGE_CITY);
            case environment_type::URBAN:
                return path_loss_urban_db(frequency_mhz, transmitter_height_m, receiver_height_m,
                                        distance_km, ant_type);
            case environment_type::SUBURBAN:
                return path_loss_suburban_db(frequency_mhz, transmitter_height_m, receiver_height_m,
                                           distance_km, ant_type);
            case environment_type::RURAL_OPEN:
                return path_loss_rural_open_db(frequency_mhz, transmitter_height_m, receiver_height_m,
                                              distance_km, ant_type);
            case environment_type::RURAL_QUASI_OPEN:
                return path_loss_rural_open_db(frequency_mhz, transmitter_height_m, receiver_height_m,
                                              distance_km, ant_type) + T(10.0);
            default:
                return path_loss_urban_db(frequency_mhz, transmitter_height_m, receiver_height_m,
                                        distance_km, ant_type);
        }
    }

    static T coverage_radius_km(T transmit_power_dbm, T transmit_gain_dbi, T receive_gain_dbi,
                              T receiver_sensitivity_dbm, T frequency_mhz, T transmitter_height_m,
                              T receiver_height_m, environment_type env,
                              antenna_type ant_type = antenna_type::SMALL_CITY) {
        const T link_budget = transmit_power_dbm + transmit_gain_dbi + receive_gain_dbi - receiver_sensitivity_dbm;

        T distance_km = T(1.0);
        T path_loss = T(0.0);

        for (int iter = 0; iter < 20; ++iter) {
            path_loss = path_loss_by_environment(frequency_mhz, transmitter_height_m, receiver_height_m,
                                               distance_km, env, ant_type);

            if (std::abs(path_loss - link_budget) < T(0.1)) {
                break;
            }

            const T error = path_loss - link_budget;
            distance_km *= std::pow(T(10.0), -error / T(44.9));

            if (distance_km < T(0.1)) distance_km = T(0.1);
            if (distance_km > T(100.0)) distance_km = T(100.0);
        }

        return distance_km;
    }

    static void path_loss_batch(const vector_type& frequencies_mhz,
                              const vector_type& transmitter_heights_m,
                              const vector_type& receiver_heights_m,
                              const vector_type& distances_km,
                              vector_type& path_losses_db,
                              environment_type env = environment_type::URBAN,
                              antenna_type ant_type = antenna_type::SMALL_CITY) {
        const size_t n = frequencies_mhz.size();
        if (transmitter_heights_m.size() != n || receiver_heights_m.size() != n ||
            distances_km.size() != n || path_losses_db.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            path_losses_db[i] = path_loss_by_environment(frequencies_mhz[i], transmitter_heights_m[i],
                                                       receiver_heights_m[i], distances_km[i], env, ant_type);
        }
    }

    static constexpr T typical_base_station_height_m() { return T(50.0); }
    static constexpr T typical_mobile_height_m() { return T(1.5); }
    static constexpr T typical_cell_radius_urban_km() { return T(2.0); }
    static constexpr T typical_cell_radius_suburban_km() { return T(5.0); }
    static constexpr T typical_cell_radius_rural_km() { return T(15.0); }

    static bool is_frequency_valid(T frequency_mhz) {
        return frequency_mhz >= T(150.0) && frequency_mhz <= T(1500.0);
    }

    static bool is_distance_valid(T distance_km) {
        return distance_km >= T(1.0) && distance_km <= T(20.0);
    }

    static bool is_transmitter_height_valid(T height_m) {
        return height_m >= T(30.0) && height_m <= T(200.0);
    }

    static bool is_receiver_height_valid(T height_m) {
        return height_m >= T(1.0) && height_m <= T(10.0);
    }
};

using hata_model_f = hata_model<float>;
using hata_model_d = hata_model<double>;

} 
} 

#endif 