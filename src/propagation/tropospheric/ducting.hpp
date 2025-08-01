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

#ifndef RVL_PROPAGATION_TROPOSPHERIC_DUCTING_HPP
#define RVL_PROPAGATION_TROPOSPHERIC_DUCTING_HPP

#include "../../core/constants.hpp"
#include "../../core/error.hpp"
#include "../../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace propagation {
namespace tropospheric {

template<typename T>
class ducting {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    enum class duct_type {
        SURFACE_DUCT,        
        ELEVATED_DUCT,       
        EVAPORATION_DUCT     
    };

    static T duct_height(T wavelength_m, T refractive_index_gradient_per_m, T refractive_index = T(1.0003)) {
        core::check_positive(wavelength_m, "Wavelength");
        core::check_positive(std::abs(refractive_index_gradient_per_m), "Refractive index gradient magnitude");
        core::check_positive(refractive_index, "Refractive index");

        const T two_pi = T(2.0) * constants::mathematical<T>::pi;
        return (wavelength_m * std::abs(refractive_index_gradient_per_m)) / (two_pi * refractive_index);
    }

    static T refractive_modulus_gradient(T temperature_gradient_k_per_m, T humidity_gradient_percent_per_m,
                                       T pressure_pa = T(101325.0), T temperature_k = T(288.15)) {
        core::check_positive(pressure_pa, "Atmospheric pressure");
        core::check_positive(temperature_k, "Temperature");

        const T temp_term = T(77.6) * pressure_pa / (temperature_k * temperature_k) * temperature_gradient_k_per_m;
        const T humidity_term = T(3.73e5) * humidity_gradient_percent_per_m / (temperature_k * temperature_k);

        return -(temp_term + humidity_term);
    }

    static T standard_atmosphere_gradient() {
        return T(-7.32e-3);
    }

    static T critical_angle_mrad(T duct_height_m, T wavelength_m) {
        core::check_positive(duct_height_m, "Duct height");
        core::check_positive(wavelength_m, "Wavelength");

        return std::sqrt(T(2.0) * wavelength_m / duct_height_m) * T(1000.0);
    }

    static T maximum_range_km(T duct_height_m, T frequency_hz, T antenna_height_m) {
        core::check_positive(duct_height_m, "Duct height");
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(antenna_height_m, "Antenna height");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T effective_earth_radius = T(8500.0);

        if (antenna_height_m > duct_height_m) {
            return T(0.0);
        }

        const T range_factor = std::sqrt(T(2.0) * effective_earth_radius * duct_height_m);
        const T frequency_factor = std::sqrt(frequency_hz / T(1e9));

        return range_factor * frequency_factor / T(1000.0);
    }

    static T evaporation_duct_height_m(T sea_surface_temperature_c, T air_temperature_c,
                                     T relative_humidity_percent, T wind_speed_ms) {
        core::check_range(sea_surface_temperature_c, T(-5.0), T(40.0), "Sea surface temperature");
        core::check_range(air_temperature_c, T(-10.0), T(45.0), "Air temperature");
        core::check_range(relative_humidity_percent, T(10.0), T(100.0), "Relative humidity");
        core::check_positive(wind_speed_ms, "Wind speed");

        const T temp_diff = air_temperature_c - sea_surface_temperature_c;
        const T humidity_factor = T(100.0) - relative_humidity_percent;

        T base_height = T(2.0) + T(0.1) * humidity_factor;

        if (temp_diff > T(0.0)) {
            base_height += T(0.5) * temp_diff;
        }

        const T wind_factor = T(1.0) + T(0.03) * wind_speed_ms;
        base_height *= wind_factor;

        return std::max(T(1.0), std::min(T(40.0), base_height));
    }

    static T ducting_probability_percent(T frequency_ghz, T path_length_km, T climate_factor = T(1.0)) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_positive(path_length_km, "Path length");
        core::check_positive(climate_factor, "Climate factor");

        const T base_probability = T(0.1) * std::pow(frequency_ghz, T(0.5)) * 
                                  std::pow(path_length_km, T(0.3));

        return std::min(T(100.0), base_probability * climate_factor);
    }

    static T signal_enhancement_db(T frequency_ghz, T duct_height_m, T antenna_height_m) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_positive(duct_height_m, "Duct height");
        core::check_positive(antenna_height_m, "Antenna height");

        if (antenna_height_m > duct_height_m) {
            return T(0.0);
        }

        const T height_factor = T(1.0) - antenna_height_m / duct_height_m;
        const T freq_factor = std::pow(frequency_ghz, T(0.3));
        const T enhancement = T(15.0) * height_factor * freq_factor;

        return std::min(T(30.0), enhancement);
    }

    static void duct_height_batch(const vector_type& wavelengths,
                                const vector_type& refractive_gradients,
                                vector_type& duct_heights,
                                T refractive_index = T(1.0003)) {
        const size_t n = wavelengths.size();
        if (refractive_gradients.size() != n || duct_heights.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        const T two_pi = T(2.0) * constants::mathematical<T>::pi;

        for (size_t i = 0; i < n; ++i) {
            if (wavelengths[i] <= T(0) || std::abs(refractive_gradients[i]) <= T(0)) {
                throw core::invalid_argument_error("Wavelengths must be positive and gradients non-zero");
            }

            duct_heights[i] = (wavelengths[i] * std::abs(refractive_gradients[i])) / 
                             (two_pi * refractive_index);
        }
    }

    static void evaporation_duct_batch(const vector_type& sea_temps,
                                     const vector_type& air_temps,
                                     const vector_type& humidities,
                                     const vector_type& wind_speeds,
                                     vector_type& duct_heights) {
        const size_t n = sea_temps.size();
        if (air_temps.size() != n || humidities.size() != n || 
            wind_speeds.size() != n || duct_heights.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            duct_heights[i] = evaporation_duct_height_m(sea_temps[i], air_temps[i],
                                                       humidities[i], wind_speeds[i]);
        }
    }

    static constexpr T typical_surface_duct_height_m() { return T(50.0); }
    static constexpr T typical_elevated_duct_height_m() { return T(200.0); }
    static constexpr T typical_evaporation_duct_height_m() { return T(15.0); }

    static constexpr T strong_ducting_gradient() { return T(-0.157e-6); }
    static constexpr T weak_ducting_gradient() { return T(-0.079e-6); }

    static constexpr T maritime_climate_factor() { return T(2.0); }
    static constexpr T continental_climate_factor() { return T(0.5); }
    static constexpr T tropical_climate_factor() { return T(3.0); }

    static bool is_ducting_conditions(T refractive_gradient) {
        return refractive_gradient < standard_atmosphere_gradient();
    }

    static bool is_super_refractive(T refractive_gradient) {
        return refractive_gradient < T(-0.079e-6) && refractive_gradient >= T(-0.157e-6);
    }

    static bool is_ducting(T refractive_gradient) {
        return refractive_gradient < T(-0.157e-6);
    }
};

using ducting_f = ducting<float>;
using ducting_d = ducting<double>;

} 
} 
} 

#endif 