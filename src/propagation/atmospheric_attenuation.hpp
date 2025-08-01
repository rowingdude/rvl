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

#ifndef RVL_PROPAGATION_ATMOSPHERIC_ATTENUATION_HPP
#define RVL_PROPAGATION_ATMOSPHERIC_ATTENUATION_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <algorithm>

namespace rvl {
namespace propagation {

template<typename T>
class atmospheric_attenuation {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    struct atmospheric_conditions {
        T temperature_k;
        T pressure_pa;
        T water_vapor_density_gm3;
        T humidity_percent;
    };

    static T calculate_oxygen_specific_attenuation_db_km(T frequency_ghz, T pressure_pa, T temperature_k) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_positive(pressure_pa, "Pressure");
        core::check_positive(temperature_k, "Temperature");

        const T f = frequency_ghz;
        const T P = pressure_pa / T(101325.0);
        const T T_inv = T(300.0) / temperature_k;

        const T f_squared = f * f;
        const T term1 = T(7.2) * T_inv * T_inv;
        const T term2 = T(0.62) * T_inv * T_inv * T_inv;

        const T numerator1 = term1 * f_squared;
        const T denominator1 = f_squared + T(0.34) * P * P * T_inv * T_inv;

        const T numerator2 = term2 * f_squared;  
        const T denominator2 = (f - T(118.75)) * (f - T(118.75)) + T(2.3) * P * P * T_inv * T_inv;

        const T numerator3 = term2 * f_squared;
        const T denominator3 = (f + T(118.75)) * (f + T(118.75)) + T(2.3) * P * P * T_inv * T_inv;

        return P * (numerator1 / denominator1 + numerator2 / denominator2 + numerator3 / denominator3) * T(1e-3);
    }

    static T calculate_water_vapor_specific_attenuation_db_km(T frequency_ghz, T water_vapor_density_gm3, 
                                                            T pressure_pa, T temperature_k) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_non_negative(water_vapor_density_gm3, "Water vapor density");
        core::check_positive(pressure_pa, "Pressure");
        core::check_positive(temperature_k, "Temperature");

        const T f = frequency_ghz;
        const T rho = water_vapor_density_gm3;
        const T P = pressure_pa / T(101325.0);
        const T T_inv = T(300.0) / temperature_k;

        if (rho < T(1e-6)) {
            return T(0);
        }

        const T f_squared = f * f;
        const T eta1 = T(0.955) * P * T_inv * T_inv + T(0.006) * rho;
        const T eta2 = T(0.735) * P * T_inv * T_inv + T(0.5) * rho;

        const T term1_num = T(3.98) * eta1 * std::exp(T(2.23) * (T(1.0) - T_inv));
        const T term1_den = (f - T(22.235)) * (f - T(22.235)) + T(9.42) * eta1 * eta1;

        const T term2_num = T(11.96) * eta2 * std::exp(T(0.7) * (T(1.0) - T_inv));  
        const T term2_den = (f - T(183.31)) * (f - T(183.31)) + T(11.14) * eta2 * eta2;

        return rho * f_squared * T_inv * T_inv * T_inv * (term1_num / term1_den + term2_num / term2_den) * T(1e-4);
    }

    static T calculate_total_gaseous_attenuation_db_km(T frequency_ghz, const atmospheric_conditions& conditions) {
        const T oxygen_atten = calculate_oxygen_specific_attenuation_db_km(
            frequency_ghz, conditions.pressure_pa, conditions.temperature_k);

        const T water_vapor_atten = calculate_water_vapor_specific_attenuation_db_km(
            frequency_ghz, conditions.water_vapor_density_gm3, 
            conditions.pressure_pa, conditions.temperature_k);

        return oxygen_atten + water_vapor_atten;
    }

    static T calculate_rain_specific_attenuation_db_km(T frequency_ghz, T rain_rate_mmh, 
                                                      T temperature_k = T(273.15), 
                                                      bool horizontal_polarization = true) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_non_negative(rain_rate_mmh, "Rain rate");
        core::check_positive(temperature_k, "Temperature");

        if (rain_rate_mmh < T(0.1)) {
            return T(0);
        }

        T a, b;
        if (frequency_ghz < T(2.9)) {
            a = horizontal_polarization ? T(0.0000387) : T(0.0000352);
            b = horizontal_polarization ? T(0.912) : T(0.880);
        } else if (frequency_ghz < T(54.0)) {
            const T log_f = std::log10(frequency_ghz);
            if (horizontal_polarization) {
                a = std::pow(T(10.0), T(-4.33) + T(0.973) * log_f + T(0.0372) * log_f * log_f);
                b = T(1.076) + T(0.0460) * log_f - T(0.0238) * log_f * log_f;
            } else {
                a = std::pow(T(10.0), T(-4.35) + T(0.972) * log_f + T(0.0334) * log_f * log_f);
                b = T(1.065) + T(0.0456) * log_f - T(0.0220) * log_f * log_f;
            }
        } else {
            a = horizontal_polarization ? T(0.187) : T(0.167);
            b = horizontal_polarization ? T(0.735) : T(0.691);
        }

        const T temperature_factor = std::pow(T(273.15) / temperature_k, T(0.5));

        return a * std::pow(rain_rate_mmh, b) * temperature_factor;
    }

    static T calculate_cloud_fog_specific_attenuation_db_km(T frequency_ghz, T liquid_water_content_gm3, 
                                                          T temperature_k = T(273.15)) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_non_negative(liquid_water_content_gm3, "Liquid water content");
        core::check_positive(temperature_k, "Temperature");

        if (liquid_water_content_gm3 < T(1e-6)) {
            return T(0);
        }

        const T f = frequency_ghz;
        const T M = liquid_water_content_gm3;
        const T theta = T(300.0) / temperature_k;

        const T epsilon_s = T(103.3) * (theta - T(1.0)) + T(77.66);
        const T epsilon_inf = T(0.0671) * epsilon_s;
        const T lambda_s = T(0.0047) + T(2.164e-5) * temperature_k;

        const T f_p = T(20.1) * lambda_s * (epsilon_s - epsilon_inf) / (T(2.0) * constants::mathematical<T>::pi);

        const T epsilon_1 = epsilon_inf + (epsilon_s - epsilon_inf) / (T(1.0) + (f / f_p) * (f / f_p));
        const T epsilon_2 = (epsilon_s - epsilon_inf) * (f / f_p) / (T(1.0) + (f / f_p) * (f / f_p));

        const T K_l = T(0.819) * f / (epsilon_2 * (T(1.0) + (epsilon_1 - T(1.0)) * (epsilon_1 - T(1.0)) / (epsilon_2 * epsilon_2)));

        return K_l * M;
    }

    static T calculate_path_attenuation_db(T frequency_ghz, T path_length_km, 
                                         const atmospheric_conditions& conditions,
                                         T rain_rate_mmh = T(0),
                                         T liquid_water_content_gm3 = T(0)) {
        const T gaseous_atten = calculate_total_gaseous_attenuation_db_km(frequency_ghz, conditions);
        const T rain_atten = calculate_rain_specific_attenuation_db_km(frequency_ghz, rain_rate_mmh, conditions.temperature_k);
        const T cloud_atten = calculate_cloud_fog_specific_attenuation_db_km(frequency_ghz, liquid_water_content_gm3, conditions.temperature_k);

        const T total_specific_atten = gaseous_atten + rain_atten + cloud_atten;

        return total_specific_atten * path_length_km;
    }

    static T calculate_water_vapor_density_from_humidity(T temperature_k, T pressure_pa, T humidity_percent) {
        core::check_positive(temperature_k, "Temperature");
        core::check_positive(pressure_pa, "Pressure");
        core::check_range(humidity_percent, T(0), T(100), "Humidity");

        const T T_c = temperature_k - T(273.15);
        const T e_s = T(6.1078) * std::pow(T(10.0), T(7.5) * T_c / (T_c + T(237.3)));
        const T e = e_s * humidity_percent / T(100.0);

        return T(216.7) * e / temperature_k;
    }

    static atmospheric_conditions create_standard_atmosphere(T height_m, T humidity_percent = T(50.0)) {
        atmospheric_conditions conditions;

        if (height_m <= T(11000.0)) {
            conditions.temperature_k = T(288.15) - T(0.0065) * height_m;
            conditions.pressure_pa = T(101325.0) * std::pow(conditions.temperature_k / T(288.15), T(5.256));
        } else {
            conditions.temperature_k = T(216.65);
            conditions.pressure_pa = T(22632.0) * std::exp(-T(0.0001577) * (height_m - T(11000.0)));
        }

        conditions.humidity_percent = std::max(T(1.0), humidity_percent * std::exp(-height_m / T(2000.0)));
        conditions.water_vapor_density_gm3 = calculate_water_vapor_density_from_humidity(
            conditions.temperature_k, conditions.pressure_pa, conditions.humidity_percent);

        return conditions;
    }

    static void calculate_gaseous_attenuation_batch(const vector_type& frequencies,
                                                  const vector_type& temperatures,
                                                  const vector_type& pressures,
                                                  const vector_type& humidities,
                                                  vector_type& attenuations) {
        const size_t n = frequencies.size();
        if (temperatures.size() != n || pressures.size() != n || 
            humidities.size() != n || attenuations.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            atmospheric_conditions conditions;
            conditions.temperature_k = temperatures[i];
            conditions.pressure_pa = pressures[i];
            conditions.humidity_percent = humidities[i];
            conditions.water_vapor_density_gm3 = calculate_water_vapor_density_from_humidity(
                temperatures[i], pressures[i], humidities[i]);

            attenuations[i] = calculate_total_gaseous_attenuation_db_km(frequencies[i], conditions);
        }
    }

    static void calculate_rain_attenuation_batch(const vector_type& frequencies,
                                               const vector_type& rain_rates,
                                               vector_type& attenuations) {
        const size_t n = frequencies.size();
        if (rain_rates.size() != n || attenuations.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            attenuations[i] = calculate_rain_specific_attenuation_db_km(frequencies[i], rain_rates[i]);
        }
    }

    static T calculate_zenith_attenuation_db(T frequency_ghz, T height_km, 
                                           const atmospheric_conditions& surface_conditions) {
        const T scale_height = (frequency_ghz > T(57.0)) ? T(2.0) : T(6.0);
        const T surface_atten = calculate_total_gaseous_attenuation_db_km(frequency_ghz, surface_conditions);

        return surface_atten * scale_height * (T(1.0) - std::exp(-height_km / scale_height));
    }

    static T calculate_slant_path_factor(T elevation_angle_deg) {
        core::check_range(elevation_angle_deg, T(0), T(90), "Elevation angle");

        const T elevation_rad = elevation_angle_deg * constants::mathematical<T>::deg_to_rad;

        if (elevation_angle_deg > T(10.0)) {
            return T(1.0) / std::sin(elevation_rad);
        } else {
            return T(1.0) / (std::sin(elevation_rad) + T(0.15) * std::pow(elevation_angle_deg + T(3.885), T(-1.253)));
        }
    }

    static constexpr T oxygen_line_frequency_ghz() { return T(60.0); }
    static constexpr T water_vapor_line_frequency_ghz() { return T(22.235); }
    static constexpr T standard_temperature_k() { return T(288.15); }
    static constexpr T standard_pressure_pa() { return T(101325.0); }
};

using atmospheric_attenuation_f = atmospheric_attenuation<float>;
using atmospheric_attenuation_d = atmospheric_attenuation<double>;

} 
} 

#endif 