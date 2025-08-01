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

#ifndef RVL_FEEDLINE_TRANSMISSION_LINE_IMPEDANCE_HPP
#define RVL_FEEDLINE_TRANSMISSION_LINE_IMPEDANCE_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace feedline {

template<typename T>
class transmission_line_impedance {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static T calculate_from_lc(T inductance_per_meter, T capacitance_per_meter) {
        core::check_positive(inductance_per_meter, "Inductance per meter");
        core::check_positive(capacitance_per_meter, "Capacitance per meter");

        return std::sqrt(inductance_per_meter / capacitance_per_meter);
    }

    static T coaxial_cable_impedance(T inner_diameter, T outer_diameter, T relative_permittivity = T(1.0)) {
        core::check_positive(inner_diameter, "Inner diameter");
        core::check_positive(outer_diameter, "Outer diameter");
        core::check_positive(relative_permittivity, "Relative permittivity");

        if (outer_diameter <= inner_diameter) {
            throw core::invalid_argument_error("Outer diameter must be greater than inner diameter");
        }

        const T impedance_0 = constants::physical<T>::eta_0;
        const T log_ratio = std::log(outer_diameter / inner_diameter);

        return (impedance_0 / (T(2.0) * constants::mathematical<T>::pi * std::sqrt(relative_permittivity))) * log_ratio;
    }

    static T parallel_wire_impedance(T wire_diameter, T spacing, T relative_permittivity = T(1.0)) {
        core::check_positive(wire_diameter, "Wire diameter");
        core::check_positive(spacing, "Wire spacing");
        core::check_positive(relative_permittivity, "Relative permittivity");

        if (spacing <= wire_diameter) {
            throw core::invalid_argument_error("Wire spacing must be greater than wire diameter");
        }

        const T impedance_0 = constants::physical<T>::eta_0;
        const T D_over_d = spacing / wire_diameter;

        T log_term;
        if (D_over_d > T(2.0)) {
            log_term = std::log(T(2.0) * D_over_d);
        } else {
            log_term = std::acosh(D_over_d);
        }

        return (impedance_0 / (constants::mathematical<T>::pi * std::sqrt(relative_permittivity))) * log_term;
    }

    static T microstrip_impedance(T trace_width, T substrate_thickness, T relative_permittivity) {
        core::check_positive(trace_width, "Trace width");
        core::check_positive(substrate_thickness, "Substrate thickness");
        core::check_positive(relative_permittivity, "Relative permittivity");

        const T w_over_h = trace_width / substrate_thickness;
        const T sqrt_er = std::sqrt(relative_permittivity);

        T z0;
        if (w_over_h <= T(1.0)) {
            const T term1 = std::log(T(8.0) / w_over_h + w_over_h / T(4.0));
            z0 = (T(60.0) / sqrt_er) * term1;
        } else {
            const T term1 = w_over_h + T(1.393) + T(0.667) * std::log(w_over_h + T(1.444));
            z0 = (T(120.0) * constants::mathematical<T>::pi) / (sqrt_er * term1);
        }

        return z0;
    }

    static T stripline_impedance(T trace_width, T substrate_thickness, T relative_permittivity) {
        core::check_positive(trace_width, "Trace width");
        core::check_positive(substrate_thickness, "Substrate thickness");
        core::check_positive(relative_permittivity, "Relative permittivity");

        const T w_over_b = trace_width / substrate_thickness;
        const T sqrt_er = std::sqrt(relative_permittivity);

        T cf;  
        if (w_over_b <= T(0.35)) {
            cf = T(30.0) - constants::physical<T>::eta_0 * constants::mathematical<T>::pi / (T(2.0) * sqrt_er);
        } else {
            cf = T(5.98) * substrate_thickness / (T(0.8) * trace_width + substrate_thickness);
        }

        const T z0 = (T(60.0) / sqrt_er) * std::log(T(4.0) / (constants::mathematical<T>::pi * w_over_b)) + cf;

        return z0;
    }

    static void coaxial_impedance_batch(const vector_type& inner_diameters,
                                      const vector_type& outer_diameters,
                                      const vector_type& relative_permittivities,
                                      vector_type& impedances) {
        const size_t n = inner_diameters.size();
        if (outer_diameters.size() != n || relative_permittivities.size() != n || impedances.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        const T impedance_0 = constants::physical<T>::eta_0;
        const T factor = impedance_0 / (T(2.0) * constants::mathematical<T>::pi);

        for (size_t i = 0; i < n; ++i) {
            if (inner_diameters[i] <= T(0) || outer_diameters[i] <= T(0) || 
                relative_permittivities[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            if (outer_diameters[i] <= inner_diameters[i]) {
                throw core::invalid_argument_error("Outer diameter must be greater than inner diameter");
            }

            const T log_ratio = std::log(outer_diameters[i] / inner_diameters[i]);
            impedances[i] = factor * log_ratio / std::sqrt(relative_permittivities[i]);
        }
    }

    static void microstrip_impedance_batch(const vector_type& trace_widths,
                                         const vector_type& substrate_thicknesses,
                                         const vector_type& relative_permittivities,
                                         vector_type& impedances) {
        const size_t n = trace_widths.size();
        if (substrate_thicknesses.size() != n || relative_permittivities.size() != n || impedances.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (trace_widths[i] <= T(0) || substrate_thicknesses[i] <= T(0) || 
                relative_permittivities[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            impedances[i] = microstrip_impedance(trace_widths[i], substrate_thicknesses[i], 
                                               relative_permittivities[i]);
        }
    }

    static constexpr T standard_50_ohm() { return constants::transmission_lines<T>::standard_50_ohm; }
    static constexpr T standard_75_ohm() { return constants::transmission_lines<T>::standard_75_ohm; }
    static constexpr T standard_300_ohm() { return constants::transmission_lines<T>::standard_300_ohm; }
    static constexpr T standard_600_ohm() { return constants::transmission_lines<T>::standard_600_ohm; }
};

using transmission_line_impedance_f = transmission_line_impedance<float>;
using transmission_line_impedance_d = transmission_line_impedance<double>;

} 
} 

#endif 