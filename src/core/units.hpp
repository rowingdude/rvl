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

#ifndef RVL_CORE_UNITS_HPP
#define RVL_CORE_UNITS_HPP

#include <cmath>
#include "constants.hpp"

namespace rvl {
namespace units {

template<typename T>
inline constexpr T db_to_linear(T db) noexcept {
    return std::exp(db * constants::radio<T>::db_to_linear_factor);
}

template<typename T>
inline constexpr T linear_to_db(T linear) noexcept {
    return std::log(linear) * constants::radio<T>::linear_to_db_factor;
}

template<typename T>
inline constexpr T dbm_to_watts(T dbm) noexcept {
    return T(0.001) * db_to_linear(dbm);
}

template<typename T>
inline constexpr T watts_to_dbm(T watts) noexcept {
    return linear_to_db(watts * T(1000.0));
}

template<typename T>
inline constexpr T frequency_to_wavelength(T freq_hz) noexcept {
    return constants::physical<T>::c / freq_hz;
}

template<typename T>
inline constexpr T wavelength_to_frequency(T wavelength_m) noexcept {
    return constants::physical<T>::c / wavelength_m;
}

template<typename T>
inline constexpr T meters_to_feet(T meters) noexcept {
    return meters * T(3.28084);
}

template<typename T>
inline constexpr T feet_to_meters(T feet) noexcept {
    return feet * T(0.3048);
}

template<typename T>
inline constexpr T meters_to_miles(T meters) noexcept {
    return meters * T(0.000621371);
}

template<typename T>
inline constexpr T miles_to_meters(T miles) noexcept {
    return miles * T(1609.344);
}

template<typename T>
inline constexpr T km_to_meters(T km) noexcept {
    return km * T(1000.0);
}

template<typename T>
inline constexpr T meters_to_km(T meters) noexcept {
    return meters * T(0.001);
}

template<typename T>
inline constexpr T mhz_to_hz(T mhz) noexcept {
    return mhz * T(1e6);
}

template<typename T>
inline constexpr T hz_to_mhz(T hz) noexcept {
    return hz * T(1e-6);
}

template<typename T>
inline constexpr T ghz_to_hz(T ghz) noexcept {
    return ghz * T(1e9);
}

template<typename T>
inline constexpr T hz_to_ghz(T hz) noexcept {
    return hz * T(1e-9);
}

template<typename T>
inline constexpr T celsius_to_kelvin(T celsius) noexcept {
    return celsius + T(273.15);
}

template<typename T>
inline constexpr T kelvin_to_celsius(T kelvin) noexcept {
    return kelvin - T(273.15);
}

template<typename T>
inline constexpr T reflection_coeff_to_vswr(T gamma) noexcept {
    T gamma_mag = std::abs(gamma);
    return (T(1.0) + gamma_mag) / (T(1.0) - gamma_mag);
}

template<typename T>
inline constexpr T vswr_to_reflection_coeff(T vswr) noexcept {
    return (vswr - T(1.0)) / (vswr + T(1.0));
}

template<typename T>
inline constexpr T vswr_to_return_loss_db(T vswr) noexcept {
    T gamma = vswr_to_reflection_coeff(vswr);
    return -T(20.0) * std::log10(gamma);
}

template<typename T>
inline constexpr T reflection_coeff_to_return_loss_db(T gamma) noexcept {
    return -T(20.0) * std::log10(std::abs(gamma));
}

template<typename T>
inline constexpr T vswr_to_mismatch_loss_db(T vswr) noexcept {
    T gamma = vswr_to_reflection_coeff(vswr);
    return -T(10.0) * std::log10(T(1.0) - gamma * gamma);
}

} 
} 

#endif 