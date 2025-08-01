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

#ifndef RVL_FEEDLINE_WAVEGUIDE_HPP
#define RVL_FEEDLINE_WAVEGUIDE_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace feedline {

template<typename T>
class waveguide {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    enum class mode_type {
        TE10,    
        TE20,
        TE01,
        TE11,
        TM11
    };

    static T cutoff_frequency_te10(T broad_dimension_a, T narrow_dimension_b = T(0.0)) {
        core::check_positive(broad_dimension_a, "Broad dimension");

        return constants::physical<T>::c / (T(2.0) * broad_dimension_a);
    }

    static T cutoff_frequency_te_mn(T broad_dimension_a, T narrow_dimension_b, int m, int n) {
        core::check_positive(broad_dimension_a, "Broad dimension");
        core::check_positive(narrow_dimension_b, "Narrow dimension");

        if (m < 0 || n < 0 || (m == 0 && n == 0)) {
            throw core::invalid_argument_error("Invalid mode indices");
        }

        const T term_m = (m == 0) ? T(0.0) : std::pow(T(m) / (T(2.0) * broad_dimension_a), T(2.0));
        const T term_n = (n == 0) ? T(0.0) : std::pow(T(n) / (T(2.0) * narrow_dimension_b), T(2.0));

        return constants::physical<T>::c * std::sqrt(term_m + term_n);
    }

    static T characteristic_impedance_te10(T frequency_hz, T broad_dimension_a) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(broad_dimension_a, "Broad dimension");

        const T fc = cutoff_frequency_te10(broad_dimension_a);

        if (frequency_hz <= fc) {
            throw core::invalid_argument_error("Frequency must be above cutoff frequency");
        }

        const T impedance_0 = constants::physical<T>::eta_0;
        const T freq_ratio = frequency_hz / fc;

        return impedance_0 / std::sqrt(T(1.0) - T(1.0) / (freq_ratio * freq_ratio));
    }

    static T guide_wavelength_te10(T frequency_hz, T broad_dimension_a) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(broad_dimension_a, "Broad dimension");

        const T fc = cutoff_frequency_te10(broad_dimension_a);

        if (frequency_hz <= fc) {
            throw core::invalid_argument_error("Frequency must be above cutoff frequency");
        }

        const T free_space_wavelength = constants::physical<T>::c / frequency_hz;
        const T freq_ratio = frequency_hz / fc;

        return free_space_wavelength / std::sqrt(T(1.0) - T(1.0) / (freq_ratio * freq_ratio));
    }

    static T attenuation_te10_np_per_m(T frequency_hz, T broad_dimension_a, T narrow_dimension_b,
                                     T surface_resistance_ohms, T intrinsic_impedance = constants::physical<T>::eta_0) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(broad_dimension_a, "Broad dimension");
        core::check_positive(narrow_dimension_b, "Narrow dimension");
        core::check_positive(surface_resistance_ohms, "Surface resistance");
        core::check_positive(intrinsic_impedance, "Intrinsic impedance");

        const T fc = cutoff_frequency_te10(broad_dimension_a);

        if (frequency_hz <= fc) {
            throw core::invalid_argument_error("Frequency must be above cutoff frequency");
        }

        const T k = T(2.0) * constants::mathematical<T>::pi * frequency_hz / constants::physical<T>::c;
        const T kc = constants::mathematical<T>::pi / broad_dimension_a;
        const T beta = std::sqrt(k * k - kc * kc);

        const T term1 = T(2.0) * surface_resistance_ohms / (intrinsic_impedance * broad_dimension_a * narrow_dimension_b);
        const T term2 = broad_dimension_a / narrow_dimension_b + narrow_dimension_b / broad_dimension_a;
        const T term3 = (k * k * broad_dimension_a * narrow_dimension_b) / (beta * beta);
        const T term4 = T(1.0) / (k * k - kc * kc);

        return term1 * k / beta * (term2 + term3 * term4);
    }

    static T attenuation_te10_db_per_m(T frequency_hz, T broad_dimension_a, T narrow_dimension_b,
                                     T surface_resistance_ohms, T intrinsic_impedance = constants::physical<T>::eta_0) {
        const T alpha_np = attenuation_te10_np_per_m(frequency_hz, broad_dimension_a, narrow_dimension_b,
                                                   surface_resistance_ohms, intrinsic_impedance);

        return alpha_np * T(20.0) * std::log10(constants::mathematical<T>::euler);
    }

    static T surface_resistance(T frequency_hz, T conductivity_s_per_m) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(conductivity_s_per_m, "Conductivity");

        const T mu_0 = constants::physical<T>::mu_0;
        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;

        return std::sqrt(omega * mu_0 / (T(2.0) * conductivity_s_per_m));
    }

    static T power_handling_kw(T broad_dimension_a, T narrow_dimension_b, T frequency_hz,
                             T breakdown_field_mv_per_m = T(30.0)) {
        core::check_positive(broad_dimension_a, "Broad dimension");
        core::check_positive(narrow_dimension_b, "Narrow dimension");
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(breakdown_field_mv_per_m, "Breakdown field");

        const T fc = cutoff_frequency_te10(broad_dimension_a);

        if (frequency_hz <= fc) {
            throw core::invalid_argument_error("Frequency must be above cutoff frequency");
        }

        const T impedance = characteristic_impedance_te10(frequency_hz, broad_dimension_a);
        const T breakdown_field_v_per_m = breakdown_field_mv_per_m * T(1e6);
        const T max_voltage = breakdown_field_v_per_m * narrow_dimension_b / T(2.0);

        const T power_watts = max_voltage * max_voltage / impedance;
        return power_watts / T(1000.0);
    }

    static bool is_single_mode(T frequency_hz, T broad_dimension_a, T narrow_dimension_b) {
        const T fc_te10 = cutoff_frequency_te10(broad_dimension_a);
        const T fc_te20 = cutoff_frequency_te_mn(broad_dimension_a, narrow_dimension_b, 2, 0);
        const T fc_te01 = cutoff_frequency_te_mn(broad_dimension_a, narrow_dimension_b, 0, 1);

        const T next_mode_fc = std::min(fc_te20, fc_te01);

        return frequency_hz > fc_te10 && frequency_hz < next_mode_fc;
    }

    static void cutoff_frequency_batch(const vector_type& broad_dimensions,
                                     vector_type& cutoff_frequencies) {
        if (broad_dimensions.size() != cutoff_frequencies.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const T c_over_2 = constants::physical<T>::c / T(2.0);
        const size_t n = broad_dimensions.size();

        for (size_t i = 0; i < n; ++i) {
            if (broad_dimensions[i] <= T(0)) {
                throw core::invalid_argument_error("All dimensions must be positive");
            }

            cutoff_frequencies[i] = c_over_2 / broad_dimensions[i];
        }
    }

    static void attenuation_batch(const vector_type& frequencies,
                                const vector_type& broad_dimensions,
                                const vector_type& narrow_dimensions,
                                const vector_type& surface_resistances,
                                vector_type& attenuations_db_per_m) {
        const size_t n = frequencies.size();
        if (broad_dimensions.size() != n || narrow_dimensions.size() != n ||
            surface_resistances.size() != n || attenuations_db_per_m.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (frequencies[i] <= T(0) || broad_dimensions[i] <= T(0) || 
                narrow_dimensions[i] <= T(0) || surface_resistances[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            attenuations_db_per_m[i] = attenuation_te10_db_per_m(frequencies[i], broad_dimensions[i],
                                                                narrow_dimensions[i], surface_resistances[i]);
        }
    }

    static constexpr T wr90_broad_mm() { return T(22.86); }      
    static constexpr T wr90_narrow_mm() { return T(10.16); }

    static constexpr T wr62_broad_mm() { return T(15.80); }      
    static constexpr T wr62_narrow_mm() { return T(7.90); }

    static constexpr T wr42_broad_mm() { return T(10.67); }      
    static constexpr T wr42_narrow_mm() { return T(4.32); }

    static constexpr T wr28_broad_mm() { return T(7.11); }       
    static constexpr T wr28_narrow_mm() { return T(3.56); }

    static constexpr T wr22_broad_mm() { return T(5.69); }       
    static constexpr T wr22_narrow_mm() { return T(2.84); }

    static constexpr T wr15_broad_mm() { return T(3.76); }       
    static constexpr T wr15_narrow_mm() { return T(1.88); }

    static constexpr T copper_conductivity() { return T(5.8e7); }
    static constexpr T aluminum_conductivity() { return T(3.5e7); }
    static constexpr T silver_conductivity() { return T(6.1e7); }

    static constexpr T air_breakdown_mv_per_m() { return T(30.0); }
    static constexpr T dry_nitrogen_breakdown_mv_per_m() { return T(35.0); }
};

using waveguide_f = waveguide<float>;
using waveguide_d = waveguide<double>;

} 
} 

#endif 