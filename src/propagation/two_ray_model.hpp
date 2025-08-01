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

#ifndef RVL_PROPAGATION_TWO_RAY_MODEL_HPP
#define RVL_PROPAGATION_TWO_RAY_MODEL_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>

namespace rvl {
namespace propagation {

template<typename T>
class two_ray_model {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;

    struct ground_reflection_parameters {
        T relative_permittivity;
        T conductivity_s_per_m;
        T surface_roughness_factor;
    };

    struct two_ray_geometry {
        T transmitter_height_m;
        T receiver_height_m;
        T distance_m;
        T operating_frequency_hz;
    };

    static complex_type calculate_ground_reflection_coefficient(T grazing_angle_rad, 
                                                              T frequency_hz,
                                                              const ground_reflection_parameters& ground,
                                                              bool horizontal_polarization = true) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_range(grazing_angle_rad, T(0), constants::mathematical<T>::pi / T(2.0), "Grazing angle");

        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;
        const T epsilon_0 = constants::physical<T>::epsilon_0;

        const complex_type epsilon_r(ground.relative_permittivity, 
                                   -ground.conductivity_s_per_m / (omega * epsilon_0));

        const T cos_theta = std::cos(grazing_angle_rad);
        const T sin_theta = std::sin(grazing_angle_rad);

        const complex_type sqrt_term = std::sqrt(epsilon_r - complex_type(sin_theta * sin_theta, T(0)));

        complex_type reflection_coeff;
        if (horizontal_polarization) {

            const complex_type numerator = complex_type(cos_theta, T(0)) - sqrt_term;
            const complex_type denominator = complex_type(cos_theta, T(0)) + sqrt_term;
            reflection_coeff = numerator / denominator;
        } else {

            const complex_type numerator = epsilon_r * complex_type(cos_theta, T(0)) - sqrt_term;
            const complex_type denominator = epsilon_r * complex_type(cos_theta, T(0)) + sqrt_term;
            reflection_coeff = numerator / denominator;
        }

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T roughness_loss = std::exp(-T(2.0) * std::pow(T(2.0) * constants::mathematical<T>::pi * 
                                                           ground.surface_roughness_factor * cos_theta / wavelength, T(2)));

        return reflection_coeff * complex_type(roughness_loss, T(0));
    }

    static T calculate_grazing_angle_rad(const two_ray_geometry& geometry) {
        const T height_diff = std::abs(geometry.transmitter_height_m - geometry.receiver_height_m);
        const T total_height = geometry.transmitter_height_m + geometry.receiver_height_m;

        return std::atan(total_height / geometry.distance_m);
    }

    static complex_type calculate_two_ray_field(const two_ray_geometry& geometry,
                                              const ground_reflection_parameters& ground,
                                              bool horizontal_polarization = true) {
        const T wavelength = constants::physical<T>::c / geometry.operating_frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;

        const T direct_path = std::sqrt(geometry.distance_m * geometry.distance_m + 
                                      std::pow(geometry.transmitter_height_m - geometry.receiver_height_m, T(2)));

        const T reflected_path = std::sqrt(geometry.distance_m * geometry.distance_m + 
                                         std::pow(geometry.transmitter_height_m + geometry.receiver_height_m, T(2)));

        const complex_type direct_component = std::exp(complex_type(T(0), -k * direct_path)) / direct_path;

        const T grazing_angle = calculate_grazing_angle_rad(geometry);
        const complex_type reflection_coeff = calculate_ground_reflection_coefficient(
            grazing_angle, geometry.operating_frequency_hz, ground, horizontal_polarization);

        const complex_type reflected_component = reflection_coeff * 
                                               std::exp(complex_type(T(0), -k * reflected_path)) / reflected_path;

        return direct_component + reflected_component;
    }

    static T calculate_two_ray_path_loss_db(const two_ray_geometry& geometry,
                                          const ground_reflection_parameters& ground,
                                          bool horizontal_polarization = true) {
        const complex_type total_field = calculate_two_ray_field(geometry, ground, horizontal_polarization);
        const T field_magnitude = std::abs(total_field);

        if (field_magnitude <= T(0)) {
            return T(300.0); 
        }

        const T wavelength = constants::physical<T>::c / geometry.operating_frequency_hz;
        const T reference_field = wavelength / (T(4.0) * constants::mathematical<T>::pi);

        return -T(20.0) * std::log10(field_magnitude / reference_field);
    }

    static T calculate_two_ray_path_loss_simplified_db(const two_ray_geometry& geometry) {

        if (geometry.distance_m <= T(4.0) * geometry.transmitter_height_m * geometry.receiver_height_m / 
                                   (constants::physical<T>::c / geometry.operating_frequency_hz)) {

            const T wavelength = constants::physical<T>::c / geometry.operating_frequency_hz;
            return T(20.0) * std::log10(T(4.0) * constants::mathematical<T>::pi * geometry.distance_m / wavelength);
        }

        return T(40.0) * std::log10(geometry.distance_m) - 
               T(20.0) * std::log10(geometry.transmitter_height_m * geometry.receiver_height_m);
    }

    static T calculate_breakpoint_distance(const two_ray_geometry& geometry) {
        const T wavelength = constants::physical<T>::c / geometry.operating_frequency_hz;
        return T(4.0) * geometry.transmitter_height_m * geometry.receiver_height_m / wavelength;
    }

    static T calculate_interference_pattern_period_m(const two_ray_geometry& geometry) {
        const T wavelength = constants::physical<T>::c / geometry.operating_frequency_hz;
        const T total_height = geometry.transmitter_height_m + geometry.receiver_height_m;

        return wavelength * geometry.distance_m / (T(2.0) * total_height);
    }

    static void calculate_path_losses_batch(const vector_type& distances,
                                          const vector_type& tx_heights,
                                          const vector_type& rx_heights,
                                          const vector_type& frequencies,
                                          T ground_permittivity,
                                          T ground_conductivity,
                                          vector_type& path_losses) {
        const size_t n = distances.size();
        if (tx_heights.size() != n || rx_heights.size() != n || 
            frequencies.size() != n || path_losses.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        ground_reflection_parameters ground;
        ground.relative_permittivity = ground_permittivity;
        ground.conductivity_s_per_m = ground_conductivity;
        ground.surface_roughness_factor = T(0.01); 

        for (size_t i = 0; i < n; ++i) {
            two_ray_geometry geometry;
            geometry.distance_m = distances[i];
            geometry.transmitter_height_m = tx_heights[i];
            geometry.receiver_height_m = rx_heights[i];
            geometry.operating_frequency_hz = frequencies[i];

            path_losses[i] = calculate_two_ray_path_loss_db(geometry, ground);
        }
    }

    static T calculate_field_strength_variation_db(const two_ray_geometry& geometry,
                                                 const ground_reflection_parameters& ground,
                                                 T height_variation_m) {

        const complex_type nominal_field = calculate_two_ray_field(geometry, ground);
        const T nominal_magnitude = std::abs(nominal_field);

        two_ray_geometry varied_geometry = geometry;
        varied_geometry.receiver_height_m += height_variation_m;

        const complex_type varied_field = calculate_two_ray_field(varied_geometry, ground);
        const T varied_magnitude = std::abs(varied_field);

        if (nominal_magnitude <= T(0) || varied_magnitude <= T(0)) {
            return T(0);
        }

        return T(20.0) * std::log10(varied_magnitude / nominal_magnitude);
    }

    static ground_reflection_parameters create_ground_parameters(const std::string& surface_type) {
        ground_reflection_parameters params;

        if (surface_type == "concrete" || surface_type == "urban") {
            params.relative_permittivity = T(5.0);
            params.conductivity_s_per_m = T(0.001);
            params.surface_roughness_factor = T(0.02);
        } else if (surface_type == "water" || surface_type == "sea") {
            params.relative_permittivity = T(81.0);
            params.conductivity_s_per_m = T(5.0);
            params.surface_roughness_factor = T(0.1);
        } else if (surface_type == "dry_ground" || surface_type == "rural") {
            params.relative_permittivity = T(4.0);
            params.conductivity_s_per_m = T(0.001);
            params.surface_roughness_factor = T(0.05);
        } else if (surface_type == "wet_ground") {
            params.relative_permittivity = T(30.0);
            params.conductivity_s_per_m = T(0.02);
            params.surface_roughness_factor = T(0.03);
        } else {

            params.relative_permittivity = T(15.0);
            params.conductivity_s_per_m = T(0.005);
            params.surface_roughness_factor = T(0.01);
        }

        return params;
    }

    static bool is_two_ray_model_applicable(const two_ray_geometry& geometry) {

        const T min_height = std::min(geometry.transmitter_height_m, geometry.receiver_height_m);
        const T wavelength = constants::physical<T>::c / geometry.operating_frequency_hz;

        if (min_height < wavelength / T(10.0)) {
            return false;
        }

        const T near_field_boundary = T(2.0) * std::max(geometry.transmitter_height_m, geometry.receiver_height_m) *
                                     std::max(geometry.transmitter_height_m, geometry.receiver_height_m) / wavelength;

        return geometry.distance_m > near_field_boundary;
    }

    static T calculate_fade_margin_db_for_height_variation(const two_ray_geometry& geometry,
                                                         const ground_reflection_parameters& ground,
                                                         T height_std_deviation_m) {

        const T variation_1_sigma = calculate_field_strength_variation_db(geometry, ground, height_std_deviation_m);
        const T variation_neg_1_sigma = calculate_field_strength_variation_db(geometry, ground, -height_std_deviation_m);

        const T rms_variation = std::sqrt((variation_1_sigma * variation_1_sigma + 
                                         variation_neg_1_sigma * variation_neg_1_sigma) / T(2.0));

        return T(2.33) * rms_variation;
    }

    static constexpr T typical_urban_permittivity() { return T(5.0); }
    static constexpr T typical_rural_permittivity() { return T(4.0); }
    static constexpr T typical_urban_conductivity() { return T(0.001); }
    static constexpr T typical_rural_conductivity() { return T(0.001); }
    static constexpr T typical_surface_roughness_m() { return T(0.01); }
};

using two_ray_model_f = two_ray_model<float>;
using two_ray_model_d = two_ray_model<double>;

} 
} 

#endif 