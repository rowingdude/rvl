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

#ifndef RVL_FEEDLINE_TWIN_LEAD_HPP
#define RVL_FEEDLINE_TWIN_LEAD_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace feedline {

template<typename T>
class twin_lead {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static T characteristic_impedance(T center_to_center_spacing_m, T wire_diameter_m, 
                                    T relative_permittivity = T(1.0)) {
        core::check_positive(center_to_center_spacing_m, "Center-to-center spacing");
        core::check_positive(wire_diameter_m, "Wire diameter");
        core::check_positive(relative_permittivity, "Relative permittivity");

        if (center_to_center_spacing_m <= wire_diameter_m) {
            throw core::invalid_argument_error("Center spacing must be greater than wire diameter");
        }

        const T spacing_ratio = center_to_center_spacing_m / wire_diameter_m;
        const T log_term = std::log10(T(2.0) * spacing_ratio);

        return (T(276.0) / std::sqrt(relative_permittivity)) * log_term;
    }

    static T characteristic_impedance_precise(T center_to_center_spacing_m, T wire_diameter_m,
                                            T relative_permittivity = T(1.0)) {
        core::check_positive(center_to_center_spacing_m, "Center-to-center spacing");
        core::check_positive(wire_diameter_m, "Wire diameter");
        core::check_positive(relative_permittivity, "Relative permittivity");

        if (center_to_center_spacing_m <= wire_diameter_m) {
            throw core::invalid_argument_error("Center spacing must be greater than wire diameter");
        }

        const T D = center_to_center_spacing_m;
        const T d = wire_diameter_m;
        const T ratio = D / d;

        T log_term;
        if (ratio > T(2.0)) {
            log_term = std::log(T(2.0) * ratio - T(1.0));
        } else {
            log_term = std::acosh(ratio);
        }

        const T impedance_0 = constants::physical<T>::eta_0;
        return (impedance_0 / (constants::mathematical<T>::pi * std::sqrt(relative_permittivity))) * log_term;
    }

    static T velocity_factor(T relative_permittivity) {
        core::check_positive(relative_permittivity, "Relative permittivity");
        return T(1.0) / std::sqrt(relative_permittivity);
    }

    static T propagation_delay_ns_per_m(T relative_permittivity) {
        const T vf = velocity_factor(relative_permittivity);
        const T c_ns_per_m = constants::physical<T>::c / T(1e9);
        return T(1.0) / (vf * c_ns_per_m);
    }

    static T capacitance_pf_per_m(T center_to_center_spacing_m, T wire_diameter_m,
                                 T relative_permittivity = T(1.0)) {
        core::check_positive(center_to_center_spacing_m, "Center-to-center spacing");
        core::check_positive(wire_diameter_m, "Wire diameter");
        core::check_positive(relative_permittivity, "Relative permittivity");

        const T ratio = center_to_center_spacing_m / wire_diameter_m;
        const T epsilon_0_pf_per_m = constants::physical<T>::epsilon_0 * T(1e12);

        const T capacitance = (constants::mathematical<T>::pi * epsilon_0_pf_per_m * relative_permittivity) /
                             std::log(T(2.0) * ratio);

        return capacitance;
    }

    static T inductance_nh_per_m(T center_to_center_spacing_m, T wire_diameter_m) {
        core::check_positive(center_to_center_spacing_m, "Center-to-center spacing");
        core::check_positive(wire_diameter_m, "Wire diameter");

        const T ratio = center_to_center_spacing_m / wire_diameter_m;
        const T mu_0_nh_per_m = constants::physical<T>::mu_0 * T(1e9);

        const T inductance = (mu_0_nh_per_m / constants::mathematical<T>::pi) * 
                            std::log(T(2.0) * ratio);

        return inductance;
    }

    static T loss_db_per_m(T frequency_hz, T center_to_center_spacing_m, T wire_diameter_m,
                          T conductivity_s_per_m = T(5.8e7), T relative_permittivity = T(1.0),
                          T loss_tangent = T(0.0)) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(center_to_center_spacing_m, "Center-to-center spacing");
        core::check_positive(wire_diameter_m, "Wire diameter");
        core::check_positive(conductivity_s_per_m, "Conductivity");

        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;
        const T mu_0 = constants::physical<T>::mu_0;

        const T skin_depth = T(1.0) / std::sqrt(constants::mathematical<T>::pi * frequency_hz * mu_0 * conductivity_s_per_m);
        const T surface_resistance = T(1.0) / (conductivity_s_per_m * skin_depth);

        const T z0 = characteristic_impedance(center_to_center_spacing_m, wire_diameter_m, relative_permittivity);

        const T conductor_loss = surface_resistance / (constants::mathematical<T>::pi * wire_diameter_m * z0);

        T dielectric_loss = T(0.0);
        if (loss_tangent > T(0.0)) {
            const T wavelength_0 = constants::physical<T>::c / frequency_hz;
            dielectric_loss = constants::mathematical<T>::pi * std::sqrt(relative_permittivity) * 
                             loss_tangent / wavelength_0;
        }

        const T total_loss_np_per_m = conductor_loss + dielectric_loss;
        return total_loss_np_per_m * T(20.0) * std::log10(constants::mathematical<T>::euler);
    }

    static T radiation_loss_db_per_m(T frequency_hz, T center_to_center_spacing_m) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(center_to_center_spacing_m, "Center-to-center spacing");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T spacing_in_wavelengths = center_to_center_spacing_m / wavelength;

        if (spacing_in_wavelengths < T(0.01)) {
            return T(0.0);
        }

        const T radiation_resistance = T(790.0) * std::pow(spacing_in_wavelengths, T(4.0));
        const T line_impedance = T(300.0);

        const T loss_per_wavelength = radiation_resistance / line_impedance;
        const T loss_per_meter = loss_per_wavelength / wavelength;

        return T(20.0) * std::log10(constants::mathematical<T>::euler) * loss_per_meter;
    }

    static void characteristic_impedance_batch(const vector_type& spacings,
                                             const vector_type& diameters,
                                             const vector_type& permittivities,
                                             vector_type& impedances) {
        const size_t n = spacings.size();
        if (diameters.size() != n || permittivities.size() != n || impedances.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (spacings[i] <= T(0) || diameters[i] <= T(0) || permittivities[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            if (spacings[i] <= diameters[i]) {
                throw core::invalid_argument_error("Spacing must be greater than diameter");
            }

            impedances[i] = characteristic_impedance(spacings[i], diameters[i], permittivities[i]);
        }
    }

    static void loss_batch(const vector_type& frequencies,
                         const vector_type& spacings,
                         const vector_type& diameters,
                         vector_type& losses_db_per_m,
                         T conductivity = T(5.8e7),
                         T relative_permittivity = T(1.0),
                         T loss_tangent = T(0.0)) {
        const size_t n = frequencies.size();
        if (spacings.size() != n || diameters.size() != n || losses_db_per_m.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (frequencies[i] <= T(0) || spacings[i] <= T(0) || diameters[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            losses_db_per_m[i] = loss_db_per_m(frequencies[i], spacings[i], diameters[i],
                                              conductivity, relative_permittivity, loss_tangent);
        }
    }

    static constexpr T standard_300_ohm_spacing_mm() { return constants::transmission_lines<T>::twin_lead_300_spacing_mm; }
    static constexpr T standard_300_ohm_diameter_mm() { return constants::transmission_lines<T>::twin_lead_300_diameter_mm; }
    static constexpr T standard_450_ohm_spacing_mm() { return constants::transmission_lines<T>::twin_lead_450_spacing_mm; }
    static constexpr T standard_450_ohm_diameter_mm() { return constants::transmission_lines<T>::twin_lead_450_diameter_mm; }
    static constexpr T standard_600_ohm_spacing_mm() { return constants::transmission_lines<T>::twin_lead_600_spacing_mm; }
    static constexpr T standard_600_ohm_diameter_mm() { return constants::transmission_lines<T>::twin_lead_600_diameter_mm; }

    static constexpr T polyethylene_permittivity() { return constants::materials<T>::polyethylene_dielectric; }
    static constexpr T teflon_permittivity() { return constants::materials<T>::teflon_dielectric; }
    static constexpr T air_permittivity() { return constants::atmospheric<T>::air_permittivity; }

    static constexpr T copper_conductivity() { return constants::materials<T>::copper_conductivity; }
    static constexpr T aluminum_conductivity() { return constants::materials<T>::aluminum_conductivity; }
};

using twin_lead_f = twin_lead<float>;
using twin_lead_d = twin_lead<double>;

} 
} 

#endif 