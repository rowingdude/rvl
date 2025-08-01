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

#ifndef RVL_ANTENNA_HELICAL_ANTENNA_HPP
#define RVL_ANTENNA_HELICAL_ANTENNA_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace antenna {

template<typename T>
class helical_antenna {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    enum class mode {
        NORMAL,     
        AXIAL       
    };

    static T gain_axial_mode(T num_turns, T turn_spacing, T circumference, T wavelength) {
        core::check_positive(num_turns, "Number of turns");
        core::check_positive(turn_spacing, "Turn spacing");
        core::check_positive(circumference, "Circumference");
        core::check_positive(wavelength, "Wavelength");

        const T c_over_lambda = circumference / wavelength;
        const T s_over_lambda = turn_spacing / wavelength;

        const T gain = T(15.0) * num_turns * s_over_lambda * 
                       c_over_lambda * c_over_lambda;

        return gain;
    }

    static T gain_axial_mode_db(T num_turns, T turn_spacing, T circumference, T wavelength) {
        const T gain_numeric = gain_axial_mode(num_turns, turn_spacing, circumference, wavelength);
        return T(10.0) * std::log10(gain_numeric);
    }

    static T optimal_circumference(T wavelength) {
        core::check_positive(wavelength, "Wavelength");
        return wavelength;
    }

    static T optimal_turn_spacing(T wavelength) {
        core::check_positive(wavelength, "Wavelength");
        return wavelength / T(4.0);
    }

    static T pitch_angle_rad(T turn_spacing, T circumference) {
        core::check_positive(turn_spacing, "Turn spacing");
        core::check_positive(circumference, "Circumference");

        return std::atan(turn_spacing / circumference);
    }

    static T pitch_angle_deg(T turn_spacing, T circumference) {
        return pitch_angle_rad(turn_spacing, circumference) * constants::mathematical<T>::rad_to_deg;
    }

    static T axial_length(T num_turns, T turn_spacing) {
        core::check_positive(num_turns, "Number of turns");
        core::check_positive(turn_spacing, "Turn spacing");

        return num_turns * turn_spacing;
    }

    static T wire_length(T num_turns, T circumference, T turn_spacing) {
        core::check_positive(num_turns, "Number of turns");
        core::check_positive(circumference, "Circumference");
        core::check_positive(turn_spacing, "Turn spacing");

        const T turn_length = std::sqrt(circumference * circumference + 
                                       turn_spacing * turn_spacing);
        return num_turns * turn_length;
    }

    static T beamwidth_deg(T num_turns, T turn_spacing, T wavelength) {
        core::check_positive(num_turns, "Number of turns");
        core::check_positive(turn_spacing, "Turn spacing");
        core::check_positive(wavelength, "Wavelength");

        const T antenna_length = num_turns * turn_spacing;
        const T length_in_wavelengths = antenna_length / wavelength;

        if (length_in_wavelengths < T(0.5)) {
            return T(90.0);
        }

        return T(52.0) / std::sqrt(length_in_wavelengths);
    }

    static T input_impedance_axial(T circumference, T wavelength) {
        core::check_positive(circumference, "Circumference");
        core::check_positive(wavelength, "Wavelength");

        const T c_over_lambda = circumference / wavelength;

        if (c_over_lambda >= T(0.75) && c_over_lambda <= T(1.33)) {
            return T(140.0);
        } else {
            return T(140.0) * c_over_lambda;
        }
    }

    static void gain_axial_batch(const vector_type& num_turns,
                               const vector_type& turn_spacings,
                               const vector_type& circumferences,
                               const vector_type& wavelengths,
                               vector_type& gains) {
        const size_t n = num_turns.size();
        if (turn_spacings.size() != n || circumferences.size() != n || 
            wavelengths.size() != n || gains.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (num_turns[i] <= T(0) || turn_spacings[i] <= T(0) || 
                circumferences[i] <= T(0) || wavelengths[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            gains[i] = gain_axial_mode(num_turns[i], turn_spacings[i], 
                                     circumferences[i], wavelengths[i]);
        }
    }

    static void gain_axial_db_batch(const vector_type& num_turns,
                                  const vector_type& turn_spacings,
                                  const vector_type& circumferences,
                                  const vector_type& wavelengths,
                                  vector_type& gains_db) {
        gain_axial_batch(num_turns, turn_spacings, circumferences, wavelengths, gains_db);

        const size_t n = gains_db.size();
        for (size_t i = 0; i < n; ++i) {
            gains_db[i] = T(10.0) * std::log10(gains_db[i]);
        }
    }

    static bool is_axial_mode_range(T circumference, T wavelength) {
        const T ratio = circumference / wavelength;
        return (ratio >= T(0.75) && ratio <= T(1.33));
    }

    static bool is_optimal_pitch_angle(T turn_spacing, T circumference) {
        const T pitch_rad = pitch_angle_rad(turn_spacing, circumference);
        const T pitch_deg = pitch_rad * constants::mathematical<T>::rad_to_deg;
        return (pitch_deg >= T(12.0) && pitch_deg <= T(16.0));
    }

    static constexpr T typical_num_turns() { return T(10.0); }
    static constexpr T typical_gain_axial_db() { return T(12.0); }
};

using helical_antenna_f = helical_antenna<float>;
using helical_antenna_d = helical_antenna<double>;

} 
} 

#endif 