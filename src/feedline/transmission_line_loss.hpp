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

#ifndef RVL_FEEDLINE_TRANSMISSION_LINE_LOSS_HPP
#define RVL_FEEDLINE_TRANSMISSION_LINE_LOSS_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace feedline {

template<typename T>
class transmission_line_loss {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static T calculate_db(T attenuation_db_per_meter, T length_meters) {
        core::check_non_negative(attenuation_db_per_meter, "Attenuation constant");
        core::check_positive(length_meters, "Length");

        return attenuation_db_per_meter * length_meters;
    }

    static T calculate_linear(T attenuation_db_per_meter, T length_meters) {
        const T loss_db = calculate_db(attenuation_db_per_meter, length_meters);
        return std::pow(T(10.0), loss_db / T(10.0));
    }

    static T attenuation_from_loss(T total_loss_db, T length_meters) {
        core::check_positive(total_loss_db, "Total loss");
        core::check_positive(length_meters, "Length");

        return total_loss_db / length_meters;
    }

    static T length_from_loss(T total_loss_db, T attenuation_db_per_meter) {
        core::check_positive(total_loss_db, "Total loss");
        core::check_positive(attenuation_db_per_meter, "Attenuation constant");

        return total_loss_db / attenuation_db_per_meter;
    }

    static T coaxial_conductor_loss_db_per_meter(T frequency_hz, T inner_diameter_m, 
                                                T outer_diameter_m, T conductivity_s_per_m) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(inner_diameter_m, "Inner diameter");
        core::check_positive(outer_diameter_m, "Outer diameter");
        core::check_positive(conductivity_s_per_m, "Conductivity");

        const T mu_0 = constants::physical<T>::mu_0;
        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;

        const T skin_depth = T(1.0) / std::sqrt(constants::mathematical<T>::pi * frequency_hz * mu_0 * conductivity_s_per_m);

        const T rs = T(1.0) / (conductivity_s_per_m * skin_depth);

        const T a = inner_diameter_m / T(2.0);
        const T b = outer_diameter_m / T(2.0);

        const T z0 = T(60.0) * std::log(b / a);

        const T loss_np_per_m = rs / (T(2.0) * z0) * (T(1.0) / a + T(1.0) / b);

        const T loss_db_per_m = loss_np_per_m * T(20.0) * std::log10(constants::mathematical<T>::euler);

        return loss_db_per_m;
    }

    static T dielectric_loss_db_per_meter(T frequency_hz, T relative_permittivity, T loss_tangent) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(relative_permittivity, "Relative permittivity");
        core::check_positive(loss_tangent, "Loss tangent");

        const T c = constants::physical<T>::c;
        const T wavelength_0 = c / frequency_hz;

        const T loss_np_per_m = constants::mathematical<T>::pi * std::sqrt(relative_permittivity) * 
                               loss_tangent / wavelength_0;

        const T loss_db_per_m = loss_np_per_m * T(20.0) * std::log10(constants::mathematical<T>::euler);

        return loss_db_per_m;
    }

    static T total_coaxial_loss_db_per_meter(T frequency_hz, T inner_diameter_m, T outer_diameter_m,
                                           T conductivity_s_per_m, T relative_permittivity, T loss_tangent) {
        const T conductor_loss = coaxial_conductor_loss_db_per_meter(frequency_hz, inner_diameter_m,
                                                                   outer_diameter_m, conductivity_s_per_m);
        const T dielectric_loss = dielectric_loss_db_per_meter(frequency_hz, relative_permittivity, loss_tangent);

        return conductor_loss + dielectric_loss;
    }

    static T matched_line_loss_with_vswr(T matched_loss_db, T vswr) {
        core::check_positive(matched_loss_db, "Matched loss");
        core::check_range(vswr, T(1.0), std::numeric_limits<T>::max(), "VSWR");

        const T gamma = (vswr - T(1.0)) / (vswr + T(1.0));
        const T matched_loss_linear = std::pow(T(10.0), matched_loss_db / T(10.0));

        const T numerator = matched_loss_linear * (T(1.0) - gamma * gamma);
        const T denominator = T(1.0) - gamma * gamma * matched_loss_linear * matched_loss_linear;

        if (denominator <= T(0)) {
            return std::numeric_limits<T>::infinity();
        }

        const T actual_loss_linear = numerator / denominator;
        return T(10.0) * std::log10(actual_loss_linear);
    }

    static void calculate_batch(const vector_type& attenuations, const vector_type& lengths,
                              vector_type& losses) {
        const size_t n = attenuations.size();
        if (lengths.size() != n || losses.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (attenuations[i] < T(0) || lengths[i] <= T(0)) {
                throw core::invalid_argument_error("Attenuations must be non-negative and lengths must be positive");
            }

            losses[i] = attenuations[i] * lengths[i];
        }
    }

    static void coaxial_conductor_loss_batch(const vector_type& frequencies,
                                           const vector_type& inner_diameters,
                                           const vector_type& outer_diameters,
                                           const vector_type& conductivities,
                                           vector_type& losses_db_per_m) {
        const size_t n = frequencies.size();
        if (inner_diameters.size() != n || outer_diameters.size() != n || 
            conductivities.size() != n || losses_db_per_m.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (frequencies[i] <= T(0) || inner_diameters[i] <= T(0) || 
                outer_diameters[i] <= T(0) || conductivities[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            losses_db_per_m[i] = coaxial_conductor_loss_db_per_meter(frequencies[i], inner_diameters[i],
                                                                   outer_diameters[i], conductivities[i]);
        }
    }

    static constexpr T rg58_loss_100mhz_db_per_m() { return constants::cable_specs<T>::rg58_loss_100mhz_db_per_m; }
    static constexpr T rg58_loss_1ghz_db_per_m() { return constants::cable_specs<T>::rg58_loss_1ghz_db_per_m; }
    static constexpr T rg213_loss_100mhz_db_per_m() { return constants::cable_specs<T>::rg213_loss_100mhz_db_per_m; }
    static constexpr T rg213_loss_1ghz_db_per_m() { return constants::cable_specs<T>::rg213_loss_1ghz_db_per_m; }

    static constexpr T copper_conductivity() { return constants::materials<T>::copper_conductivity; }
    static constexpr T silver_conductivity() { return constants::materials<T>::silver_conductivity; }
    static constexpr T aluminum_conductivity() { return constants::materials<T>::aluminum_conductivity; }
};

using transmission_line_loss_f = transmission_line_loss<float>;
using transmission_line_loss_d = transmission_line_loss<double>;

} 
} 

#endif 