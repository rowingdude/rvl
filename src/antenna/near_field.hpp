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

#ifndef RVL_ANTENNA_NEAR_FIELD_HPP
#define RVL_ANTENNA_NEAR_FIELD_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>

namespace rvl {
namespace antenna {

template<typename T>
class near_field {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    enum class field_region {
        REACTIVE_NEAR_FIELD,    
        RADIATING_NEAR_FIELD,   
        FAR_FIELD              
    };

    static T fraunhofer_distance(T max_antenna_dimension, T wavelength) {
        core::check_positive(max_antenna_dimension, "Maximum antenna dimension");
        core::check_positive(wavelength, "Wavelength");

        return T(2.0) * max_antenna_dimension * max_antenna_dimension / wavelength;
    }

    static T fresnel_distance(T max_antenna_dimension, T wavelength) {
        core::check_positive(max_antenna_dimension, "Maximum antenna dimension");
        core::check_positive(wavelength, "Wavelength");

        const T d_cubed = max_antenna_dimension * max_antenna_dimension * max_antenna_dimension;
        return T(0.62) * std::sqrt(d_cubed / wavelength);
    }

    static T reactive_near_field_boundary(T max_antenna_dimension) {
        core::check_positive(max_antenna_dimension, "Maximum antenna dimension");

        return T(0.159) * max_antenna_dimension;
    }

    static field_region classify_field_region(T distance, T max_antenna_dimension, T wavelength) {
        core::check_positive(distance, "Distance");
        core::check_positive(max_antenna_dimension, "Maximum antenna dimension");
        core::check_positive(wavelength, "Wavelength");

        const T reactive_boundary = reactive_near_field_boundary(max_antenna_dimension);
        const T fresnel_boundary = fresnel_distance(max_antenna_dimension, wavelength);
        const T fraunhofer_boundary = fraunhofer_distance(max_antenna_dimension, wavelength);

        if (distance < reactive_boundary || distance < fresnel_boundary) {
            return field_region::REACTIVE_NEAR_FIELD;
        } else if (distance < fraunhofer_boundary) {
            return field_region::RADIATING_NEAR_FIELD;
        } else {
            return field_region::FAR_FIELD;
        }
    }

    static T minimum_measurement_distance(T max_antenna_dimension, T wavelength, 
                                        T accuracy_requirement_db = T(1.0)) {
        core::check_positive(max_antenna_dimension, "Maximum antenna dimension");
        core::check_positive(wavelength, "Wavelength");
        core::check_positive(accuracy_requirement_db, "Accuracy requirement");

        T multiplier;
        if (accuracy_requirement_db <= T(0.1)) {
            multiplier = T(10.0);
        } else if (accuracy_requirement_db <= T(0.5)) {
            multiplier = T(5.0);
        } else if (accuracy_requirement_db <= T(1.0)) {
            multiplier = T(2.0);
        } else {
            multiplier = T(1.0);
        }

        return multiplier * fraunhofer_distance(max_antenna_dimension, wavelength);
    }

    static T phase_center_location_uncertainty(T distance, T max_antenna_dimension, T wavelength) {
        core::check_positive(distance, "Distance");
        core::check_positive(max_antenna_dimension, "Maximum antenna dimension");
        core::check_positive(wavelength, "Wavelength");

        const T far_field_distance = fraunhofer_distance(max_antenna_dimension, wavelength);

        if (distance >= far_field_distance) {
            return T(0.0);
        }

        const T ratio = distance / far_field_distance;
        return max_antenna_dimension * (T(1.0) - ratio) / T(8.0);
    }

    static T aperture_efficiency_from_near_field(T measured_gain, T physical_area, T wavelength) {
        core::check_positive(measured_gain, "Measured gain");
        core::check_positive(physical_area, "Physical area");
        core::check_positive(wavelength, "Wavelength");

        const T theoretical_gain = T(4.0) * constants::mathematical<T>::pi * physical_area / 
                                  (wavelength * wavelength);

        return measured_gain / theoretical_gain;
    }

    static T near_field_amplitude_correction_db(T distance, T max_antenna_dimension, T wavelength) {
        core::check_positive(distance, "Distance");
        core::check_positive(max_antenna_dimension, "Maximum antenna dimension");
        core::check_positive(wavelength, "Wavelength");

        const T far_field_distance = fraunhofer_distance(max_antenna_dimension, wavelength);

        if (distance >= far_field_distance) {
            return T(0.0);
        }

        const T ratio = distance / far_field_distance;
        const T correction_factor = T(1.0) - T(0.25) * (T(1.0) - ratio);

        return T(20.0) * std::log10(correction_factor);
    }

    static void fraunhofer_distance_batch(const vector_type& max_dimensions,
                                        const vector_type& wavelengths,
                                        vector_type& distances) {
        const size_t n = max_dimensions.size();
        if (wavelengths.size() != n || distances.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (max_dimensions[i] <= T(0) || wavelengths[i] <= T(0)) {
                throw core::invalid_argument_error("All dimensions and wavelengths must be positive");
            }

            distances[i] = T(2.0) * max_dimensions[i] * max_dimensions[i] / wavelengths[i];
        }
    }

    static void field_region_classification_batch(const vector_type& distances,
                                                const vector_type& max_dimensions,
                                                const vector_type& wavelengths,
                                                std::vector<field_region>& regions) {
        const size_t n = distances.size();
        if (max_dimensions.size() != n || wavelengths.size() != n || regions.size() != n) {
            throw core::dimension_mismatch_error("All vectors must have same size");
        }

        for (size_t i = 0; i < n; ++i) {
            if (distances[i] <= T(0) || max_dimensions[i] <= T(0) || wavelengths[i] <= T(0)) {
                throw core::invalid_argument_error("All parameters must be positive");
            }

            regions[i] = classify_field_region(distances[i], max_dimensions[i], wavelengths[i]);
        }
    }

    static constexpr T typical_antenna_range_distance_multiplier() { return T(2.0); }
    static constexpr T compact_range_distance_multiplier() { return T(0.5); }
    static constexpr T planar_near_field_scan_distance_multiplier() { return T(0.1); }

    static bool is_electrically_small(T max_dimension, T wavelength) {
        return (max_dimension / wavelength) < T(0.5);
    }

    static bool is_electrically_large(T max_dimension, T wavelength) {
        return (max_dimension / wavelength) > T(10.0);
    }
};

using near_field_f = near_field<float>;
using near_field_d = near_field<double>;

} 
} 

#endif 