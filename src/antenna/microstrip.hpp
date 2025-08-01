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

#ifndef RVL_ANTENNA_MICROSTRIP_HPP
#define RVL_ANTENNA_MICROSTRIP_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace antenna {

template<typename T>
class microstrip {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    static T effective_dielectric_constant(T relative_permittivity, T width, T height) {
        core::check_positive(relative_permittivity, "Relative permittivity");
        core::check_positive(width, "Width");
        core::check_positive(height, "Height");

        const T w_over_h = width / height;

        if (w_over_h >= T(1.0)) {
            const T term1 = (relative_permittivity + T(1.0)) / T(2.0);
            const T term2 = (relative_permittivity - T(1.0)) / T(2.0);
            const T term3 = T(1.0) / std::sqrt(T(1.0) + T(12.0) / w_over_h);

            return term1 + term2 * term3;
        } else {
            const T term1 = (relative_permittivity + T(1.0)) / T(2.0);
            const T term2 = (relative_permittivity - T(1.0)) / T(2.0);
            const T term3 = T(1.0) / std::sqrt(T(1.0) + T(12.0) * height / width);
            const T term4 = T(0.04) * std::pow(T(1.0) - width / height, T(2.0));

            return term1 + term2 * (term3 + term4);
        }
    }

    static T resonant_frequency(T patch_length, T effective_permittivity) {
        core::check_positive(patch_length, "Patch length");
        core::check_positive(effective_permittivity, "Effective permittivity");

        return constants::physical<T>::c / (T(2.0) * patch_length * std::sqrt(effective_permittivity));
    }

    static T resonant_frequency_with_substrate(T patch_length, T relative_permittivity, 
                                             T width, T height) {
        const T eff_permittivity = effective_dielectric_constant(relative_permittivity, width, height);
        return resonant_frequency(patch_length, eff_permittivity);
    }

    static T patch_length_for_frequency(T frequency, T effective_permittivity) {
        core::check_positive(frequency, "Frequency");
        core::check_positive(effective_permittivity, "Effective permittivity");

        return constants::physical<T>::c / (T(2.0) * frequency * std::sqrt(effective_permittivity));
    }

    static T patch_width(T frequency, T relative_permittivity, T height) {
        core::check_positive(frequency, "Frequency");
        core::check_positive(relative_permittivity, "Relative permittivity");
        core::check_positive(height, "Height");

        const T c = constants::physical<T>::c;
        const T width = c / (T(2.0) * frequency) * std::sqrt(T(2.0) / (relative_permittivity + T(1.0)));

        return width;
    }

    static T effective_length(T physical_length, T width, T height, T relative_permittivity) {
        core::check_positive(physical_length, "Physical length");
        core::check_positive(width, "Width");
        core::check_positive(height, "Height");
        core::check_positive(relative_permittivity, "Relative permittivity");

        const T eff_permittivity = effective_dielectric_constant(relative_permittivity, width, height);
        const T w_over_h = width / height;

        const T delta_l = T(0.412) * height * 
                         (eff_permittivity + T(0.3)) * (w_over_h + T(0.264)) /
                         ((eff_permittivity - T(0.258)) * (w_over_h + T(0.8)));

        return physical_length + T(2.0) * delta_l;
    }

    static T directivity_db() {
        return T(6.0);
    }

    static T efficiency_typical() {
        return T(0.9);
    }

    static T gain_db(T efficiency = T(0.9)) {
        core::check_range(efficiency, T(0.0), T(1.0), "Efficiency");
        return directivity_db() + T(10.0) * std::log10(efficiency);
    }

    static T bandwidth_typical_percent() {
        return T(2.0);
    }

    static T input_impedance_center_feed() {
        return T(200.0);
    }

    static T input_impedance_edge_feed() {
        return T(300.0);
    }

    static T inset_feed_position(T patch_length, T desired_impedance = T(50.0)) {
        core::check_positive(patch_length, "Patch length");
        core::check_positive(desired_impedance, "Desired impedance");

        const T impedance_at_center = input_impedance_center_feed();
        const T ratio = desired_impedance / impedance_at_center;

        if (ratio > T(1.0)) {
            return T(0.0);
        }

        return patch_length * T(0.5) * (T(1.0) - std::sqrt(ratio));
    }

    static void resonant_frequency_batch(const vector_type& patch_lengths,
                                       const vector_type& relative_permittivities,
                                       const vector_type& widths,
                                       const vector_type& heights,
                                       vector_type& frequencies) {
        const size_t n = patch_lengths.size();
        if (relative_permittivities.size() != n || widths.size() != n || 
            heights.size() != n || frequencies.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (patch_lengths[i] <= T(0) || relative_permittivities[i] <= T(0) ||
                widths[i] <= T(0) || heights[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            frequencies[i] = resonant_frequency_with_substrate(patch_lengths[i], 
                                                             relative_permittivities[i],
                                                             widths[i], heights[i]);
        }
    }

    static void effective_dielectric_batch(const vector_type& relative_permittivities,
                                         const vector_type& widths,
                                         const vector_type& heights,
                                         vector_type& effective_permittivities) {
        const size_t n = relative_permittivities.size();
        if (widths.size() != n || heights.size() != n || effective_permittivities.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (relative_permittivities[i] <= T(0) || widths[i] <= T(0) || heights[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            effective_permittivities[i] = effective_dielectric_constant(relative_permittivities[i],
                                                                      widths[i], heights[i]);
        }
    }

    static constexpr T typical_substrate_height_mm() { return constants::antenna_design<T>::typical_substrate_height_mm; }
    static constexpr T fr4_permittivity() { return constants::materials<T>::pcb_fr4_dielectric; }
    static constexpr T rogers_permittivity() { return constants::materials<T>::rogers_dielectric; }
    static constexpr T alumina_permittivity() { return constants::materials<T>::alumina_dielectric; }
};

using microstrip_f = microstrip<float>;
using microstrip_d = microstrip<double>;

} 
} 

#endif 