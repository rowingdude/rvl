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

#ifndef RVL_FEEDLINE_SMITH_CHART_HPP
#define RVL_FEEDLINE_SMITH_CHART_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <vector>

namespace rvl {
namespace feedline {

template<typename T>
class smith_chart {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    struct smith_point {
        complex_type impedance;         
        complex_type reflection_coeff;  
        T resistance;                   
        T reactance;                    
        T vswr;                         
        T return_loss_db;               
        bool is_normalized;             
    };

    struct transmission_line_params {
        T characteristic_impedance;     
        T electrical_length_degrees;   
        T velocity_factor;              
        T loss_db_per_100m;            
        T frequency_hz;                 
    };

    struct matching_network {
        std::string network_type;       
        std::vector<complex_type> component_values; 
        std::vector<std::string> component_types;   
        T bandwidth_hz;                 
        T q_factor;                     
    };

    static complex_type impedance_to_reflection_coefficient(const complex_type& impedance,
                                                          T characteristic_impedance = T(50.0)) {
        core::check_positive(characteristic_impedance, "Characteristic impedance");

        const complex_type z0(characteristic_impedance, T(0.0));
        const complex_type numerator = impedance - z0;
        const complex_type denominator = impedance + z0;

        if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
            throw core::numerical_error("Denominator too small in reflection coefficient calculation");
        }

        return numerator / denominator;
    }

    static complex_type reflection_coefficient_to_impedance(const complex_type& gamma,
                                                          T characteristic_impedance = T(50.0)) {
        core::check_positive(characteristic_impedance, "Characteristic impedance");

        const complex_type one(T(1.0), T(0.0));
        const complex_type numerator = (one + gamma) * characteristic_impedance;
        const complex_type denominator = one - gamma;

        if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
            throw core::numerical_error("Denominator too small in impedance calculation");
        }

        return numerator / denominator;
    }

    static complex_type transform_impedance_through_line(const complex_type& load_impedance,
                                                       const transmission_line_params& line) {
        core::check_positive(line.characteristic_impedance, "Characteristic impedance");
        core::check_range(line.velocity_factor, T(0.1), T(1.0), "Velocity factor");

        const T z0 = line.characteristic_impedance;
        const complex_type zl = load_impedance;
        const complex_type z0_complex(z0, T(0.0));

        const T beta_l = line.electrical_length_degrees * constants::mathematical<T>::pi / T(180.0);

        const T alpha = line.loss_db_per_100m * constants::mathematical<T>::ln_10 / (T(20.0) * T(100.0));
        const complex_type gamma_line = complex_type(alpha, beta_l);

        const complex_type tanh_gamma_l = std::tanh(gamma_line);
        const complex_type numerator = z0_complex * (zl + z0_complex * tanh_gamma_l);
        const complex_type denominator = z0_complex + zl * tanh_gamma_l;

        if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
            throw core::numerical_error("Denominator too small in line transformation");
        }

        return numerator / denominator;
    }

    static complex_type transform_impedance_lossless(const complex_type& load_impedance,
                                                   T electrical_length_degrees,
                                                   T characteristic_impedance = T(50.0)) {
        const T beta_l = electrical_length_degrees * constants::mathematical<T>::pi / T(180.0);
        const T tan_beta_l = std::tan(beta_l);
        const complex_type z0(characteristic_impedance, T(0.0));
        const complex_type j_z0_tan = complex_type(T(0.0), characteristic_impedance * tan_beta_l);

        const complex_type numerator = load_impedance + j_z0_tan;
        const complex_type denominator = z0 + complex_type(T(0.0), load_impedance.imag() * tan_beta_l) +
                                       complex_type(load_impedance.real() * tan_beta_l, T(0.0));

        if (std::abs(denominator) < std::numeric_limits<T>::epsilon()) {
            return complex_type(std::numeric_limits<T>::infinity(), T(0.0));
        }

        return characteristic_impedance * characteristic_impedance * denominator / numerator;
    }

    static complex_type normalize_impedance(const complex_type& impedance,
                                          T characteristic_impedance = T(50.0)) {
        core::check_positive(characteristic_impedance, "Characteristic impedance");
        return impedance / characteristic_impedance;
    }

    static complex_type denormalize_impedance(const complex_type& normalized_impedance,
                                            T characteristic_impedance = T(50.0)) {
        core::check_positive(characteristic_impedance, "Characteristic impedance");
        return normalized_impedance * characteristic_impedance;
    }

    static smith_point calculate_smith_point(const complex_type& impedance,
                                           T characteristic_impedance = T(50.0)) {
        smith_point point;

        point.impedance = impedance;
        point.resistance = impedance.real();
        point.reactance = impedance.imag();
        point.is_normalized = false;

        point.reflection_coeff = impedance_to_reflection_coefficient(impedance, characteristic_impedance);

        const T gamma_magnitude = std::abs(point.reflection_coeff);
        if (gamma_magnitude >= T(1.0)) {
            point.vswr = std::numeric_limits<T>::infinity();
            point.return_loss_db = T(0.0);
        } else {
            point.vswr = (T(1.0) + gamma_magnitude) / (T(1.0) - gamma_magnitude);
            point.return_loss_db = -T(20.0) * std::log10(gamma_magnitude);
        }

        return point;
    }

    static matching_network design_l_network(const complex_type& source_impedance,
                                           const complex_type& load_impedance,
                                           T frequency_hz) {
        core::check_positive(frequency_hz, "Frequency");

        matching_network network;
        network.network_type = "L";

        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;

        const T Rs = source_impedance.real();
        const T Xs = source_impedance.imag();
        const T Rl = load_impedance.real();
        const T Xl = load_impedance.imag();

        if (Rs <= T(0) || Rl <= T(0)) {
            throw core::invalid_argument_error("Both source and load must have positive resistance");
        }

        bool source_higher = Rs > Rl;

        T R_high = source_higher ? Rs : Rl;
        T R_low = source_higher ? Rl : Rs;
        T X_high = source_higher ? Xs : Xl;
        T X_low = source_higher ? Xl : Xs;

        T Q_required = std::sqrt(R_high / R_low - T(1.0));
        network.q_factor = Q_required;

        if (source_higher) {

            T X_series = -X_low + std::sqrt(R_low * (R_high - R_low));

            T X_parallel = (R_high * R_high + (X_high + X_series) * (X_high + X_series)) / 
                          (X_high + X_series);

            network.component_values.push_back(complex_type(T(0.0), X_series));
            network.component_values.push_back(complex_type(T(0.0), X_parallel));

            if (X_series > T(0)) {
                network.component_types.push_back("series_L");
                network.component_values[0] = complex_type(X_series / omega, T(0.0));
            } else {
                network.component_types.push_back("series_C");
                network.component_values[0] = complex_type(-T(1.0) / (omega * X_series), T(0.0));
            }

            if (X_parallel > T(0)) {
                network.component_types.push_back("shunt_L");
                network.component_values[1] = complex_type(X_parallel / omega, T(0.0));
            } else {
                network.component_types.push_back("shunt_C");
                network.component_values[1] = complex_type(-T(1.0) / (omega * X_parallel), T(0.0));
            }
        } else {

            T X_series = -X_high + std::sqrt(R_high * (R_low - R_high));
            T X_parallel = (R_low * R_low + (X_low + X_series) * (X_low + X_series)) / 
                          (X_low + X_series);

            network.component_values.push_back(complex_type(T(0.0), X_series));
            network.component_values.push_back(complex_type(T(0.0), X_parallel));

            if (X_series > T(0)) {
                network.component_types.push_back("series_L");
                network.component_values[0] = complex_type(X_series / omega, T(0.0));
            } else {
                network.component_types.push_back("series_C");
                network.component_values[0] = complex_type(-T(1.0) / (omega * X_series), T(0.0));
            }

            if (X_parallel > T(0)) {
                network.component_types.push_back("shunt_L");
                network.component_values[1] = complex_type(X_parallel / omega, T(0.0));
            } else {
                network.component_types.push_back("shunt_C");
                network.component_values[1] = complex_type(-T(1.0) / (omega * X_parallel), T(0.0));
            }
        }

        network.bandwidth_hz = frequency_hz / (T(2.0) * Q_required);

        return network;
    }

    static T design_quarter_wave_transformer(const complex_type& source_impedance,
                                           const complex_type& load_impedance) {
        const T Rs = source_impedance.real();
        const T Rl = load_impedance.real();

        if (Rs <= T(0) || Rl <= T(0)) {
            throw core::invalid_argument_error("Both impedances must have positive real parts");
        }

        return std::sqrt(Rs * Rl);
    }

    static matching_network design_stub_match(const complex_type& load_impedance,
                                            T line_impedance = T(50.0),
                                            T frequency_hz = T(1e9)) {
        matching_network network;
        network.network_type = "stub";

        const complex_type normalized_load = normalize_impedance(load_impedance, line_impedance);
        const T Rl = normalized_load.real();
        const T Xl = normalized_load.imag();

        if (Rl <= T(0)) {
            throw core::invalid_argument_error("Load resistance must be positive");
        }

        T y = T(1.0) / Rl; 
        T b = -Xl / (Rl * Rl + Xl * Xl); 

        T B_target = std::sqrt(y * (y - T(1.0)));

        T B1 = B_target;
        T B2 = -B_target;

        T d1 = std::atan(B1) * T(180.0) / constants::mathematical<T>::pi;
        T d2 = std::atan(B2) * T(180.0) / constants::mathematical<T>::pi;

        T stub_length = std::atan(-B1) * T(180.0) / constants::mathematical<T>::pi;
        T line_length = d1;

        if (stub_length < T(0)) stub_length += T(180.0);
        if (line_length < T(0)) line_length += T(180.0);

        network.component_values.push_back(complex_type(line_length, T(0.0)));
        network.component_values.push_back(complex_type(stub_length, T(0.0)));
        network.component_types.push_back("transmission_line");
        network.component_types.push_back("short_circuit_stub");

        network.q_factor = std::abs(B1);
        network.bandwidth_hz = frequency_hz / network.q_factor;

        return network;
    }

    static std::vector<complex_type> calculate_vswr_circle(T vswr, size_t num_points = 360) {
        core::check_range(vswr, T(1.0), std::numeric_limits<T>::max(), "VSWR");
        core::check_positive(num_points, "Number of points");

        std::vector<complex_type> circle_points;
        circle_points.reserve(num_points);

        const T gamma_magnitude = (vswr - T(1.0)) / (vswr + T(1.0));
        const T angle_step = T(2.0) * constants::mathematical<T>::pi / T(num_points);

        for (size_t i = 0; i < num_points; ++i) {
            const T angle = T(i) * angle_step;
            const complex_type gamma = gamma_magnitude * std::exp(complex_type(T(0.0), angle));

            const complex_type z_norm = reflection_coefficient_to_impedance(gamma, T(1.0));
            circle_points.push_back(z_norm);
        }

        return circle_points;
    }

    static std::vector<complex_type> calculate_resistance_circle(T resistance, 
                                                               size_t num_points = 180) {
        core::check_non_negative(resistance, "Resistance");
        core::check_positive(num_points, "Number of points");

        std::vector<complex_type> circle_points;
        circle_points.reserve(num_points);

        const T center_r = resistance / (T(1.0) + resistance);
        const T radius = T(1.0) / (T(1.0) + resistance);

        const T angle_step = constants::mathematical<T>::pi / T(num_points - 1);

        for (size_t i = 0; i < num_points; ++i) {
            const T angle = T(i) * angle_step - constants::mathematical<T>::half_pi;
            const T x_offset = radius * std::cos(angle);
            const T y_offset = radius * std::sin(angle);

            circle_points.push_back(complex_type(center_r + x_offset, y_offset));
        }

        return circle_points;
    }

    static std::vector<complex_type> calculate_reactance_circle(T reactance,
                                                              size_t num_points = 180) {
        core::check_positive(num_points, "Number of points");

        std::vector<complex_type> circle_points;
        circle_points.reserve(num_points);

        if (std::abs(reactance) < std::numeric_limits<T>::epsilon()) {

            const T r_step = T(10.0) / T(num_points - 1);
            for (size_t i = 0; i < num_points; ++i) {
                circle_points.push_back(complex_type(T(i) * r_step, T(0.0)));
            }
            return circle_points;
        }

        const T center_x = T(1.0);
        const T center_y = T(1.0) / reactance;
        const T radius = T(1.0) / std::abs(reactance);

        const T angle_step = constants::mathematical<T>::pi / T(num_points - 1);
        const T start_angle = (reactance > T(0)) ? T(0.0) : constants::mathematical<T>::pi;

        for (size_t i = 0; i < num_points; ++i) {
            const T angle = start_angle + T(i) * angle_step;
            const T x_offset = radius * std::cos(angle);
            const T y_offset = radius * std::sin(angle);

            circle_points.push_back(complex_type(center_x + x_offset, center_y + y_offset));
        }

        return circle_points;
    }

    static void transform_impedances_batch(const complex_vector_type& load_impedances,
                                         const transmission_line_params& line,
                                         complex_vector_type& input_impedances) {
        if (load_impedances.size() != input_impedances.size()) {
            throw core::dimension_mismatch_error("Input and output vectors must have same size");
        }

        const size_t n = load_impedances.size();
        for (size_t i = 0; i < n; ++i) {
            input_impedances[i] = transform_impedance_through_line(load_impedances[i], line);
        }
    }

    static constexpr T standard_impedance_50_ohm() { return constants::transmission_lines<T>::standard_50_ohm; }
    static constexpr T standard_impedance_75_ohm() { return constants::transmission_lines<T>::standard_75_ohm; }
    static constexpr T standard_impedance_300_ohm() { return constants::transmission_lines<T>::standard_300_ohm; }
    static constexpr T standard_impedance_600_ohm() { return constants::transmission_lines<T>::standard_600_ohm; }
    static constexpr T maximum_practical_vswr() { return constants::vswr_thresholds<T>::maximum_practical; }
    static constexpr T excellent_match_vswr() { return constants::vswr_thresholds<T>::excellent_match; }
    static constexpr T good_match_vswr() { return constants::vswr_thresholds<T>::good_match; }
    static constexpr T acceptable_match_vswr() { return constants::vswr_thresholds<T>::acceptable_amateur; }
};

using smith_chart_f = smith_chart<float>;
using smith_chart_d = smith_chart<double>;

} 
} 

#endif 