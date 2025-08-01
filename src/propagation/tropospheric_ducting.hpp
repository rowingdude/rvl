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

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <algorithm>

namespace rvl {
namespace propagation {

template<typename T>
class tropospheric_ducting {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    struct atmospheric_profile {
        T height_m;
        T temperature_k;
        T pressure_pa;
        T humidity_percent;
        T refractive_index;
    };

    static T calculate_refractivity_n_units(T temperature_k, T pressure_pa, T water_vapor_pressure_pa) {
        core::check_positive(temperature_k, "Temperature");
        core::check_positive(pressure_pa, "Pressure");
        core::check_non_negative(water_vapor_pressure_pa, "Water vapor pressure");

        const T dry_term = T(77.6) * pressure_pa / temperature_k;
        const T wet_term = T(373000.0) * water_vapor_pressure_pa / (temperature_k * temperature_k);

        return dry_term + wet_term;
    }

    static T calculate_modified_refractivity(T refractivity_n, T height_m) {
        core::check_non_negative(height_m, "Height");

        const T earth_radius = constants::physical<T>::earth_radius;
        const T curvature_term = height_m * T(1e6) / earth_radius;

        return refractivity_n + curvature_term;
    }

    static T calculate_m_factor_gradient(T height_lower_m, T height_upper_m, 
                                       T m_lower, T m_upper) {
        const T height_diff = height_upper_m - height_lower_m;
        if (std::abs(height_diff) < T(1e-10)) {
            return T(0);
        }

        return (m_upper - m_lower) / height_diff;
    }

    static bool is_ducting_condition(T m_gradient_per_km) {
        return m_gradient_per_km < T(-157.0);
    }

    static bool is_super_refraction(T m_gradient_per_km) {
        return m_gradient_per_km < T(0.0) && m_gradient_per_km >= T(-157.0);
    }

    static T calculate_duct_height(const std::vector<atmospheric_profile>& profile) {
        if (profile.size() < 2) {
            return T(0);
        }

        for (size_t i = 1; i < profile.size(); ++i) {
            const auto& lower = profile[i-1];
            const auto& upper = profile[i];

            const T m_lower = calculate_modified_refractivity(
                calculate_refractivity_n_units(lower.temperature_k, lower.pressure_pa, 
                                              lower.pressure_pa * lower.humidity_percent / T(100.0) * T(0.01)),
                lower.height_m);

            const T m_upper = calculate_modified_refractivity(
                calculate_refractivity_n_units(upper.temperature_k, upper.pressure_pa,
                                              upper.pressure_pa * upper.humidity_percent / T(100.0) * T(0.01)),
                upper.height_m);

            const T gradient = calculate_m_factor_gradient(lower.height_m, upper.height_m, m_lower, m_upper);
            const T gradient_per_km = gradient * T(1000.0);

            if (!is_ducting_condition(gradient_per_km)) {
                return upper.height_m;
            }
        }

        return profile.back().height_m;
    }

    static T calculate_evaporation_duct_height(T sea_surface_temp_k, T air_temp_k, 
                                             T wind_speed_ms, T humidity_percent) {
        core::check_positive(sea_surface_temp_k, "Sea surface temperature");
        core::check_positive(air_temp_k, "Air temperature");
        core::check_non_negative(wind_speed_ms, "Wind speed");
        core::check_range(humidity_percent, T(0), T(100), "Humidity");

        const T temp_diff = sea_surface_temp_k - air_temp_k;
        const T wind_factor = T(1.0) + T(0.1) * wind_speed_ms;
        const T humidity_factor = T(1.0) - humidity_percent / T(100.0);

        const T base_height = T(10.0) + T(2.0) * temp_diff;

        return base_height * wind_factor * (T(1.0) + humidity_factor);
    }

    static T calculate_duct_strength(T m_gradient_per_km) {
        if (!is_ducting_condition(m_gradient_per_km)) {
            return T(0);
        }

        return std::abs(m_gradient_per_km + T(157.0));
    }

    static T calculate_ducting_range_km(T frequency_ghz, T duct_height_m, T duct_strength) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_positive(duct_height_m, "Duct height");
        core::check_non_negative(duct_strength, "Duct strength");

        if (duct_strength < T(1.0)) {
            return T(0);
        }

        const T wavelength_m = constants::physical<T>::c / (frequency_ghz * T(1e9));
        const T duct_parameter = duct_height_m / wavelength_m;

        const T range_factor = std::sqrt(duct_strength) * duct_parameter;

        return T(50.0) * std::pow(range_factor, T(0.75)) / std::sqrt(frequency_ghz);
    }

    static T calculate_signal_enhancement_db(T frequency_ghz, T path_length_km, 
                                           T duct_height_m, T duct_strength) {
        if (duct_strength < T(1.0)) {
            return T(0);
        }

        const T max_range = calculate_ducting_range_km(frequency_ghz, duct_height_m, duct_strength);

        if (path_length_km > max_range) {
            return T(0);
        }

        const T enhancement_factor = duct_strength * (T(1.0) - path_length_km / max_range);

        return T(10.0) * std::log10(T(1.0) + enhancement_factor);
    }

    static T calculate_ducting_probability_percent(T frequency_ghz, T path_length_km, 
                                                 T latitude_deg, bool maritime_path = false) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_positive(path_length_km, "Path length");
        core::check_range(latitude_deg, T(-90), T(90), "Latitude");

        T base_probability;
        if (maritime_path) {
            base_probability = T(15.0) - T(0.1) * std::abs(latitude_deg);
        } else {
            base_probability = T(5.0) - T(0.05) * std::abs(latitude_deg);
        }

        base_probability = std::max(base_probability, T(0.1));

        const T frequency_factor = T(1.0) / (T(1.0) + frequency_ghz / T(10.0));
        const T distance_factor = std::exp(-path_length_km / T(200.0));

        return base_probability * frequency_factor * distance_factor;
    }

    static T calculate_multipath_fading_db(T frequency_ghz, T path_length_km, 
                                         T roughness_factor = T(1.0)) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_positive(path_length_km, "Path length");
        core::check_positive(roughness_factor, "Roughness factor");

        const T wavelength_m = constants::physical<T>::c / (frequency_ghz * T(1e9));
        const T path_difference_m = roughness_factor * std::sqrt(path_length_km * T(1000.0) * wavelength_m);

        const T phase_difference_rad = T(2.0) * constants::mathematical<T>::pi * path_difference_m / wavelength_m;

        const T fading_amplitude = T(2.0) * std::sin(phase_difference_rad / T(2.0));

        return T(20.0) * std::log10(std::abs(fading_amplitude) + T(1e-10));
    }

    static void calculate_duct_heights_batch(const vector_type& sea_temps,
                                           const vector_type& air_temps,
                                           const vector_type& wind_speeds,
                                           const vector_type& humidities,
                                           vector_type& duct_heights) {
        const size_t n = sea_temps.size();
        if (air_temps.size() != n || wind_speeds.size() != n || 
            humidities.size() != n || duct_heights.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (sea_temps[i] <= T(0) || air_temps[i] <= T(0) || wind_speeds[i] < T(0) ||
                humidities[i] < T(0) || humidities[i] > T(100)) {
                throw core::invalid_argument_error("Invalid parameters at index");
            }

            duct_heights[i] = calculate_evaporation_duct_height(sea_temps[i], air_temps[i],
                                                              wind_speeds[i], humidities[i]);
        }
    }

    static void calculate_ducting_ranges_batch(const vector_type& frequencies,
                                             const vector_type& duct_heights,
                                             const vector_type& duct_strengths,
                                             vector_type& ranges) {
        const size_t n = frequencies.size();
        if (duct_heights.size() != n || duct_strengths.size() != n || ranges.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (frequencies[i] <= T(0) || duct_heights[i] <= T(0) || duct_strengths[i] < T(0)) {
                throw core::invalid_argument_error("Invalid parameters at index");
            }

            ranges[i] = calculate_ducting_range_km(frequencies[i], duct_heights[i], duct_strengths[i]);
        }
    }

    static atmospheric_profile create_standard_atmosphere_profile(T height_m) {
        atmospheric_profile profile;
        profile.height_m = height_m;

        if (height_m <= T(11000.0)) {
            profile.temperature_k = T(288.15) - T(0.0065) * height_m;
            profile.pressure_pa = T(101325.0) * std::pow(profile.temperature_k / T(288.15), T(5.256));
        } else {
            profile.temperature_k = T(216.65);
            profile.pressure_pa = T(22632.0) * std::exp(-T(0.0001577) * (height_m - T(11000.0)));
        }

        profile.humidity_percent = std::max(T(1.0), T(80.0) * std::exp(-height_m / T(2000.0)));

        return profile;
    }

    static bool is_anomalous_propagation(T m_gradient_per_km) {
        return m_gradient_per_km < T(-157.0);
    }

    static constexpr T standard_atmosphere_gradient() { return T(118.0); }
    static constexpr T typical_evaporation_duct_height() { return T(15.0); }
    static constexpr T strong_ducting_threshold() { return T(50.0); }
};

using tropospheric_ducting_f = tropospheric_ducting<float>;
using tropospheric_ducting_d = tropospheric_ducting<double>;

} 
} 

#endif 