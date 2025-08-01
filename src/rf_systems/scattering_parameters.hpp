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

#ifndef RVL_RF_SYSTEMS_SCATTERING_PARAMETERS_HPP
#define RVL_RF_SYSTEMS_SCATTERING_PARAMETERS_HPP

#include "../core/constants.hpp"
#include "../core/units.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <complex>
#include <cmath>

namespace rvl {
namespace rf_systems {

template<typename T>
class scattering_parameters {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    static constexpr T default_impedance = T(50.0);

    static complex_type s11_from_impedance(const complex_type& z_in, T z0 = default_impedance) {
        core::check_positive(z0, "Characteristic impedance");

        const complex_type numerator = z_in - complex_type(z0, T(0.0));
        const complex_type denominator = z_in + complex_type(z0, T(0.0));

        if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
            throw core::numerical_error("Denominator too small in S11 calculation");
        }

        return numerator / denominator;
    }

    static complex_type impedance_from_s11(const complex_type& s11, T z0 = default_impedance) {
        core::check_positive(z0, "Characteristic impedance");

        const complex_type one(T(1.0), T(0.0));
        const complex_type numerator = (one + s11) * z0;
        const complex_type denominator = one - s11;

        if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
            throw core::numerical_error("Denominator too small in impedance calculation");
        }

        return numerator / denominator;
    }

    static T reflection_coefficient_magnitude(const complex_type& s11) {
        return std::abs(s11);
    }

    static T reflection_coefficient_phase_rad(const complex_type& s11) {
        return std::arg(s11);
    }

    static T reflection_coefficient_phase_deg(const complex_type& s11) {
        return std::arg(s11) * constants::mathematical<T>::rad_to_deg;
    }

    static T return_loss_db(const complex_type& s11) {
        const T magnitude = std::abs(s11);
        if (magnitude <= T(0)) {
            return T(100.0);
        }
        return -T(20.0) * std::log10(magnitude);
    }

    static T vswr_from_s11(const complex_type& s11) {
        const T gamma = std::abs(s11);
        if (gamma >= T(1.0)) {
            return std::numeric_limits<T>::infinity();
        }
        return (T(1.0) + gamma) / (T(1.0) - gamma);
    }

    static complex_type s11_from_vswr(T vswr, T phase_rad = T(0.0)) {
        core::check_range(vswr, T(1.0), std::numeric_limits<T>::max(), "VSWR");

        const T gamma = (vswr - T(1.0)) / (vswr + T(1.0));
        return gamma * std::exp(complex_type(T(0.0), phase_rad));
    }

    static T mismatch_loss_db(const complex_type& s11) {
        const T gamma_squared = std::norm(s11);
        return -T(10.0) * std::log10(T(1.0) - gamma_squared);
    }

    static void s11_from_impedance_batch(const complex_vector_type& z_in, 
                                       complex_vector_type& s11,
                                       T z0 = default_impedance) {
        if (z_in.size() != s11.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = z_in.size();
        const complex_type z0_complex(z0, T(0.0));

        for (size_t i = 0; i < n; ++i) {
            const complex_type numerator = z_in[i] - z0_complex;
            const complex_type denominator = z_in[i] + z0_complex;

            if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
                throw core::numerical_error("Denominator too small in batch S11 calculation");
            }

            s11[i] = numerator / denominator;
        }
    }

    static void return_loss_batch(const complex_vector_type& s11, vector_type& return_loss_db) {
        if (s11.size() != return_loss_db.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = s11.size();
        for (size_t i = 0; i < n; ++i) {
            const T magnitude = std::abs(s11[i]);
            if (magnitude <= T(0)) {
                return_loss_db[i] = T(100.0);
            } else {
                return_loss_db[i] = -T(20.0) * std::log10(magnitude);
            }
        }
    }

    static void vswr_batch(const complex_vector_type& s11, vector_type& vswr) {
        if (s11.size() != vswr.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = s11.size();
        for (size_t i = 0; i < n; ++i) {
            const T gamma = std::abs(s11[i]);
            if (gamma >= T(1.0)) {
                vswr[i] = std::numeric_limits<T>::infinity();
            } else {
                vswr[i] = (T(1.0) + gamma) / (T(1.0) - gamma);
            }
        }
    }
};

using scattering_parameters_f = scattering_parameters<float>;
using scattering_parameters_d = scattering_parameters<double>;

} 
} 

#endif 