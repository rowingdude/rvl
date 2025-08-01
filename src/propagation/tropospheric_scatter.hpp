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

#ifndef RVL_PROPAGATION_TROPOSPHERIC_SCATTER_HPP
#define RVL_PROPAGATION_TROPOSPHERIC_SCATTER_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <vector>
#include <limits>
#include <algorithm>

namespace rvl {
namespace propagation {

template<typename T>
class tropospheric_scatter {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;

    struct scatter_geometry {
        T transmitter_height_m;
        T receiver_height_m; 
        T path_distance_km;
        T scattering_angle_rad;
        T common_volume_height_m;
        T transmitter_lat_deg;        
        T transmitter_lon_deg;        
        T receiver_lat_deg;           
        T receiver_lon_deg;           
    };

    struct atmospheric_conditions {
        T temperature_k;               
        T pressure_hpa;                
        T humidity_percent;            
        T turbulence_strength_cn2;     
        T outer_scale_length_m;        
        T refractive_index_gradient;   
    };

    struct scatter_analysis {
        T path_loss_db;                
        T basic_transmission_loss_db;  
        T atmospheric_absorption_db;   
        T scattering_loss_db;          
        T fade_margin_db;              
        T availability_percent;        
        T coherence_bandwidth_hz;      
        T coherence_time_ms;           
        T signal_std_deviation_db;     
    };

    struct reliability_statistics {
        T availability_50_percent_db;  
        T availability_90_percent_db;  
        T availability_99_percent_db;  
        vector_type time_percentiles;  
        vector_type signal_levels_db;  
    };

    static T calculate_scattering_cross_section_per_volume(T wavenumber, T structure_constant, 
                                                         T outer_scale_m, T scattering_angle_rad) {
        core::check_positive(wavenumber, "Wavenumber");
        core::check_positive(structure_constant, "Structure constant");
        core::check_positive(outer_scale_m, "Outer scale");
        core::check_range(scattering_angle_rad, T(0), constants::mathematical<T>::pi, "Scattering angle");

        const T k = wavenumber;
        const T L = outer_scale_m;
        const T theta = scattering_angle_rad;
        const T Cn2 = structure_constant;

        const T sin_half_theta = std::sin(theta / T(2.0));
        const T q_squared = T(4.0) * k * k * sin_half_theta * sin_half_theta;

        const T spectrum_factor = T(7.0 * M_PI / 3.0) * std::pow(L, T(7.0/3.0));
        const T frequency_factor = std::pow(T(1.0) + q_squared * L * L / T(4.0), T(-7.0/6.0));

        return spectrum_factor * Cn2 * frequency_factor / (k * k);
    }

    static T calculate_bistatic_radar_cross_section(T wavenumber, T structure_constant,
                                                  T outer_scale_m, T scattering_angle_rad) {
        return T(7.0 * M_PI / 3.0) * structure_constant * 
               std::pow(outer_scale_m, T(7.0)) * 
               std::pow(std::sin(scattering_angle_rad / T(2.0)), T(2.0)) /
               (wavenumber * wavenumber);
    }

    static T calculate_common_volume_height(T tx_height_m, T rx_height_m, T distance_km) {
        core::check_positive(tx_height_m, "Transmitter height");
        core::check_positive(rx_height_m, "Receiver height");
        core::check_positive(distance_km, "Distance");

        const T earth_radius = constants::physical<T>::earth_radius;
        const T distance_m = distance_km * T(1000.0);

        const T horizon_tx = std::sqrt(T(2.0) * earth_radius * tx_height_m);
        const T horizon_rx = std::sqrt(T(2.0) * earth_radius * rx_height_m);

        if (distance_m <= horizon_tx + horizon_rx) {
            return std::max(tx_height_m, rx_height_m) + T(1000.0);
        } else {
            const T effective_height = (tx_height_m + rx_height_m) / T(2.0);
            const T additional_height = distance_m * distance_m / (T(8.0) * earth_radius);
            return effective_height + additional_height;
        }
    }

    static T calculate_scattering_angle(T tx_height_m, T rx_height_m, T distance_km, T common_volume_height_m) {
        const T distance_m = distance_km * T(1000.0);
        const T half_distance = distance_m / T(2.0);

        const T angle_tx = std::atan((common_volume_height_m - tx_height_m) / half_distance);
        const T angle_rx = std::atan((common_volume_height_m - rx_height_m) / half_distance);

        return angle_tx + angle_rx;
    }

    static T calculate_effective_scattering_area(const scatter_geometry& geometry, 
                                               T frequency_hz, T structure_constant,
                                               T outer_scale_m = T(1000.0)) {
        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T wavenumber = T(2.0) * constants::mathematical<T>::pi / wavelength;

        const T sigma_s = calculate_scattering_cross_section_per_volume(
            wavenumber, structure_constant, outer_scale_m, geometry.scattering_angle_rad);

        const T scattering_volume = calculate_scattering_volume(geometry);

        return sigma_s * scattering_volume;
    }

    static T calculate_scattering_volume(const scatter_geometry& geometry) {
        const T distance_m = geometry.path_distance_km * T(1000.0);
        const T wavelength_factor = T(100.0);

        const T volume_length = distance_m / T(4.0);
        const T volume_radius = std::sqrt(wavelength_factor * distance_m);

        return constants::mathematical<T>::pi * volume_radius * volume_radius * volume_length;
    }

    static T calculate_troposcatter_path_loss_db(T frequency_hz, const scatter_geometry& geometry,
                                               T tx_antenna_gain_db, T rx_antenna_gain_db,
                                               T structure_constant = T(1e-14)) {
        core::check_positive(frequency_hz, "Frequency");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T distance_m = geometry.path_distance_km * T(1000.0);

        const T free_space_loss_db = T(20.0) * std::log10(T(4.0) * constants::mathematical<T>::pi * distance_m / wavelength);

        const T effective_area = calculate_effective_scattering_area(geometry, frequency_hz, structure_constant);

        const T scattering_loss_db = T(-10.0) * std::log10(effective_area / (wavelength * wavelength));

        const T antenna_gain_total = tx_antenna_gain_db + rx_antenna_gain_db;

        return free_space_loss_db + scattering_loss_db - antenna_gain_total;
    }

    static T calculate_troposcatter_availability_percent(T frequency_ghz, T path_distance_km,
                                                       T climate_factor = T(1.0)) {
        core::check_positive(frequency_ghz, "Frequency");
        core::check_positive(path_distance_km, "Path distance");
        core::check_positive(climate_factor, "Climate factor");

        T base_availability;
        if (path_distance_km <= T(200.0)) {
            base_availability = T(99.9);
        } else if (path_distance_km <= T(500.0)) {
            base_availability = T(99.5);
        } else {
            base_availability = T(99.0);
        }

        const T frequency_factor = T(1.0) - T(0.01) * (frequency_ghz - T(1.0));
        const T distance_factor = std::exp(-path_distance_km / T(300.0));

        return base_availability * frequency_factor * distance_factor * climate_factor;
    }

    static T calculate_signal_fading_db(T time_percentage, T frequency_ghz, T path_distance_km) {
        core::check_range(time_percentage, T(0.001), T(50.0), "Time percentage");
        core::check_positive(frequency_ghz, "Frequency");
        core::check_positive(path_distance_km, "Path distance");

        const T log_p = std::log10(time_percentage);
        const T frequency_factor = T(20.0) * std::log10(frequency_ghz);
        const T distance_factor = T(10.0) * std::log10(path_distance_km / T(100.0));

        return T(-5.0) * log_p + frequency_factor + distance_factor;
    }

    static T calculate_aperture_antenna_gain_db(T frequency_hz, T antenna_diameter_m, T efficiency = T(0.6)) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(antenna_diameter_m, "Antenna diameter");
        core::check_range(efficiency, T(0.1), T(1.0), "Efficiency");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T area = constants::mathematical<T>::pi * antenna_diameter_m * antenna_diameter_m / T(4.0);
        const T aperture_efficiency = T(4.0) * constants::mathematical<T>::pi * area * efficiency;

        return T(10.0) * std::log10(aperture_efficiency / (wavelength * wavelength));
    }

    static scatter_analysis calculate_enhanced_scatter_analysis(T frequency_hz,
                                                              const scatter_geometry& geometry,
                                                              const atmospheric_conditions& atmosphere,
                                                              T tx_antenna_gain_dbi = T(0),
                                                              T rx_antenna_gain_dbi = T(0)) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(geometry.path_distance_km, "Path distance");
        core::check_positive(atmosphere.turbulence_strength_cn2, "Turbulence strength");

        scatter_analysis analysis;

        analysis.path_loss_db = calculate_troposcatter_path_loss_db(
            frequency_hz, geometry, tx_antenna_gain_dbi, rx_antenna_gain_dbi, 
            atmosphere.turbulence_strength_cn2);

        const T wavelength_m = constants::physical<T>::c / frequency_hz;
        const T distance_m = geometry.path_distance_km * T(1000.0);

        analysis.basic_transmission_loss_db = T(20.0) * std::log10(
            T(4.0) * constants::mathematical<T>::pi * distance_m / wavelength_m);

        T frequency_ghz = frequency_hz / T(1e9);
        analysis.atmospheric_absorption_db = T(0.1) * frequency_ghz * geometry.path_distance_km;

        analysis.scattering_loss_db = analysis.path_loss_db - analysis.basic_transmission_loss_db - 
                                    analysis.atmospheric_absorption_db + tx_antenna_gain_dbi + rx_antenna_gain_dbi;

        T path_length_factor = std::log10(geometry.path_distance_km / T(100.0));
        T turbulence_factor = std::log10(atmosphere.turbulence_strength_cn2 / T(1e-14));
        analysis.fade_margin_db = T(15.0) + T(5.0) * path_length_factor + T(3.0) * turbulence_factor;

        T availability_factor = analysis.fade_margin_db / T(20.0);
        analysis.availability_percent = T(50.0) + T(40.0) * std::tanh(availability_factor);
        analysis.availability_percent = std::min(analysis.availability_percent, T(99.9));

        analysis.coherence_bandwidth_hz = calculate_coherence_bandwidth(
            geometry.path_distance_km, atmosphere.turbulence_strength_cn2);
        analysis.coherence_time_ms = calculate_coherence_time(
            frequency_hz, atmosphere.turbulence_strength_cn2);

        analysis.signal_std_deviation_db = calculate_signal_standard_deviation(
            geometry.path_distance_km, atmosphere.turbulence_strength_cn2);

        return analysis;
    }

    static reliability_statistics calculate_availability_statistics(
        const scatter_analysis& base_analysis,
        size_t num_percentiles = 100) {

        reliability_statistics stats;

        stats.availability_50_percent_db = base_analysis.path_loss_db;
        stats.availability_90_percent_db = base_analysis.path_loss_db + 
                                         T(1.28) * base_analysis.signal_std_deviation_db;
        stats.availability_99_percent_db = base_analysis.path_loss_db + 
                                         T(2.33) * base_analysis.signal_std_deviation_db;

        stats.time_percentiles.resize(num_percentiles);
        stats.signal_levels_db.resize(num_percentiles);

        for (size_t i = 0; i < num_percentiles; ++i) {
            T percentile = T(i + 1) / T(num_percentiles) * T(100.0);
            stats.time_percentiles[i] = percentile;

            T z_score = inverse_normal_cdf(percentile / T(100.0));
            stats.signal_levels_db[i] = base_analysis.path_loss_db + 
                                      z_score * base_analysis.signal_std_deviation_db;
        }

        return stats;
    }

    static T calculate_frequency_diversity_improvement(T frequency1_hz,
                                                     T frequency2_hz,
                                                     const scatter_geometry& geometry,
                                                     const atmospheric_conditions& atmosphere) {
        core::check_positive(frequency1_hz, "Frequency 1");
        core::check_positive(frequency2_hz, "Frequency 2");

        T coherence_bw = calculate_coherence_bandwidth(
            geometry.path_distance_km, atmosphere.turbulence_strength_cn2);

        T frequency_separation = std::abs(frequency2_hz - frequency1_hz);

        T correlation_coefficient = std::exp(-frequency_separation / coherence_bw);

        T diversity_gain_db = T(10.0) * std::log10(T(2.0) / (T(1.0) + correlation_coefficient));

        return std::min(diversity_gain_db, T(6.0)); 
    }

    static T calculate_space_diversity_improvement(T antenna_separation_m,
                                                 T frequency_hz,
                                                 const scatter_geometry& geometry) {
        core::check_positive(antenna_separation_m, "Antenna separation");
        core::check_positive(frequency_hz, "Frequency");

        const T wavelength_m = constants::physical<T>::c / frequency_hz;

        T fresnel_radius = std::sqrt(wavelength_m * geometry.path_distance_km * T(1000.0) / T(2.0));

        T spatial_correlation = std::exp(-std::pow(antenna_separation_m / fresnel_radius, T(2.0)));

        T diversity_gain_db = T(10.0) * std::log10(T(2.0) / (T(1.0) + spatial_correlation));

        return std::min(diversity_gain_db, T(5.0)); 
    }

    static T find_optimal_frequency(const vector_type& test_frequencies_hz,
                                  const scatter_geometry& geometry,
                                  const atmospheric_conditions& atmosphere,
                                  T required_availability_percent = T(90.0)) {

        T best_frequency = test_frequencies_hz[0];
        T best_margin = -std::numeric_limits<T>::max();

        for (const auto& freq : test_frequencies_hz) {
            auto analysis = calculate_enhanced_scatter_analysis(freq, geometry, atmosphere);
            auto stats = calculate_availability_statistics(analysis);

            T z_score = inverse_normal_cdf(required_availability_percent / T(100.0));
            T signal_level_at_availability = analysis.path_loss_db + 
                                            z_score * analysis.signal_std_deviation_db;

            T margin = signal_level_at_availability - T(-100.0);

            if (margin > best_margin) {
                best_margin = margin;
                best_frequency = freq;
            }
        }

        return best_frequency;
    }

    static std::vector<scatter_analysis> analyze_frequency_sweep(
        const vector_type& frequencies_hz,
        const scatter_geometry& geometry,
        const atmospheric_conditions& atmosphere,
        T tx_antenna_gain_dbi = T(0),
        T rx_antenna_gain_dbi = T(0)) {

        std::vector<scatter_analysis> results;
        results.reserve(frequencies_hz.size());

        for (const auto& freq : frequencies_hz) {
            results.push_back(calculate_enhanced_scatter_analysis(
                freq, geometry, atmosphere, tx_antenna_gain_dbi, rx_antenna_gain_dbi));
        }

        return results;
    }

    static T calculate_diversity_improvement_db(T correlation_coefficient) {
        core::check_range(correlation_coefficient, T(0), T(1), "Correlation coefficient");

        if (correlation_coefficient > T(0.99)) {
            return T(0);
        }

        return T(-10.0) * std::log10(T(2.0) * (T(1.0) + correlation_coefficient));
    }

    static void calculate_path_losses_batch(const vector_type& frequencies,
                                          const vector_type& distances,
                                          const vector_type& tx_gains,
                                          const vector_type& rx_gains,
                                          vector_type& path_losses) {
        const size_t n = frequencies.size();
        if (distances.size() != n || tx_gains.size() != n || 
            rx_gains.size() != n || path_losses.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (frequencies[i] <= T(0) || distances[i] <= T(0)) {
                throw core::invalid_argument_error("Invalid parameters at index");
            }

            scatter_geometry geom;
            geom.transmitter_height_m = T(30.0);
            geom.receiver_height_m = T(30.0);
            geom.path_distance_km = distances[i];
            geom.common_volume_height_m = calculate_common_volume_height(
                geom.transmitter_height_m, geom.receiver_height_m, distances[i]);
            geom.scattering_angle_rad = calculate_scattering_angle(
                geom.transmitter_height_m, geom.receiver_height_m, 
                distances[i], geom.common_volume_height_m);

            path_losses[i] = calculate_troposcatter_path_loss_db(
                frequencies[i], geom, tx_gains[i], rx_gains[i]);
        }
    }

    static T calculate_minimum_takeoff_angle_rad(T tx_height_m, T distance_km) {
        const T earth_radius = constants::physical<T>::earth_radius;
        const T distance_m = distance_km * T(1000.0);

        const T horizon_distance = std::sqrt(T(2.0) * earth_radius * tx_height_m);

        if (distance_m <= horizon_distance) {
            return T(0);
        }

        const T angle = std::asin((distance_m - horizon_distance) / distance_m);
        return std::max(angle, T(0.001));
    }

    static bool is_beyond_horizon(T tx_height_m, T rx_height_m, T distance_km) {
        const T earth_radius = constants::physical<T>::earth_radius;
        const T horizon_distance = std::sqrt(T(2.0) * earth_radius * tx_height_m) +
                                 std::sqrt(T(2.0) * earth_radius * rx_height_m);

        return distance_km * T(1000.0) > horizon_distance;
    }

private:

    static T calculate_coherence_bandwidth(T distance_km, T cn2) {

        T path_factor = std::pow(distance_km / T(100.0), T(-0.6));
        T turbulence_factor = std::pow(cn2 / T(1e-14), T(-0.4));

        return T(1e6) * path_factor * turbulence_factor; 
    }

    static T calculate_coherence_time(T frequency_hz, T cn2) {

        T frequency_factor = std::pow(frequency_hz / T(1e9), T(-0.3));
        T turbulence_factor = std::pow(cn2 / T(1e-14), T(-0.5));

        return T(100.0) * frequency_factor * turbulence_factor; 
    }

    static T calculate_signal_standard_deviation(T distance_km, T cn2) {

        T distance_factor = std::log10(distance_km / T(100.0));
        T turbulence_factor = std::log10(cn2 / T(1e-14));

        return T(5.0) + T(2.0) * distance_factor + T(1.5) * turbulence_factor;
    }

    static T inverse_normal_cdf(T probability) {

        if (probability <= T(0.0)) return -T(10.0);
        if (probability >= T(1.0)) return T(10.0);

        static const T a0 = T(2.50662823884);
        static const T a1 = T(-18.61500062529);
        static const T a2 = T(41.39119773534);
        static const T a3 = T(-25.44106049637);

        static const T b0 = T(-8.47351093090);
        static const T b1 = T(23.08336743743);
        static const T b2 = T(-21.06224101826);
        static const T b3 = T(3.13082909833);

        T u = probability - T(0.5);
        T r = u * u;

        if (std::abs(u) < T(0.42)) {
            return u * (a0 + a1*r + a2*r*r + a3*r*r*r) / 
                   (T(1.0) + b0*r + b1*r*r + b2*r*r*r + b3*r*r*r*r);
        }

        T sign = (u > T(0)) ? T(1.0) : T(-1.0);
        r = (u > T(0)) ? (T(1.0) - probability) : probability;
        r = std::log(-std::log(r));

        return sign * (a0 + a1*r + a2*r*r + a3*r*r*r);
    }

public:

    static atmospheric_conditions create_standard_atmosphere(T height_km = T(1.0),
                                                           T season_factor = T(1.0)) {
        atmospheric_conditions conditions;

        conditions.temperature_k = T(288.15) - T(6.5) * height_km; 
        conditions.pressure_hpa = T(1013.25) * std::pow(conditions.temperature_k / T(288.15), T(5.26));
        conditions.humidity_percent = T(60.0) * std::exp(-height_km / T(2.0)); 

        conditions.refractive_index_gradient = -T(40.0) * season_factor; 

        conditions.turbulence_strength_cn2 = T(1e-14) * season_factor * 
                                           std::exp(-height_km / T(1.5));
        conditions.outer_scale_length_m = T(100.0);

        return conditions;
    }

    static scatter_geometry create_path_geometry(T transmitter_lat_deg,
                                               T transmitter_lon_deg,
                                               T transmitter_height_m,
                                               T receiver_lat_deg,
                                               T receiver_lon_deg,
                                               T receiver_height_m) {
        scatter_geometry geometry;

        geometry.transmitter_lat_deg = transmitter_lat_deg;
        geometry.transmitter_lon_deg = transmitter_lon_deg;
        geometry.receiver_lat_deg = receiver_lat_deg;
        geometry.receiver_lon_deg = receiver_lon_deg;
        geometry.transmitter_height_m = transmitter_height_m;
        geometry.receiver_height_m = receiver_height_m;

        const T deg_to_rad = constants::mathematical<T>::pi / T(180.0);
        const T earth_radius_km = constants::physical<T>::earth_radius / T(1000.0);

        T tx_lat_rad = transmitter_lat_deg * deg_to_rad;
        T tx_lon_rad = transmitter_lon_deg * deg_to_rad;
        T rx_lat_rad = receiver_lat_deg * deg_to_rad;
        T rx_lon_rad = receiver_lon_deg * deg_to_rad;

        T delta_lat = rx_lat_rad - tx_lat_rad;
        T delta_lon = rx_lon_rad - tx_lon_rad;
        T a = std::sin(delta_lat/T(2.0)) * std::sin(delta_lat/T(2.0)) +
              std::cos(tx_lat_rad) * std::cos(rx_lat_rad) *
              std::sin(delta_lon/T(2.0)) * std::sin(delta_lon/T(2.0));
        T c = T(2.0) * std::atan2(std::sqrt(a), std::sqrt(T(1.0) - a));
        geometry.path_distance_km = earth_radius_km * c;

        geometry.common_volume_height_m = calculate_common_volume_height(
            transmitter_height_m, receiver_height_m, geometry.path_distance_km);
        geometry.scattering_angle_rad = calculate_scattering_angle(
            transmitter_height_m, receiver_height_m, geometry.path_distance_km, 
            geometry.common_volume_height_m);

        return geometry;
    }

    static constexpr T typical_structure_constant_continental() { return T(1e-14); }
    static constexpr T typical_structure_constant_maritime() { return T(5e-15); }
    static constexpr T typical_structure_constant_strong_turbulence() { return T(1e-13); }
    static constexpr T typical_structure_constant_weak_turbulence() { return T(1e-15); }
    static constexpr T typical_outer_scale_m() { return T(1000.0); }
    static constexpr T minimum_elevation_angle_deg() { return T(0.1); }
    static constexpr T maximum_practical_distance_km() { return T(800.0); }
    static constexpr T minimum_practical_distance_km() { return T(100.0); }
    static constexpr T typical_fade_margin_db() { return T(30.0); }
};

using tropospheric_scatter_f = tropospheric_scatter<float>;
using tropospheric_scatter_d = tropospheric_scatter<double>;

} 
} 

#endif 