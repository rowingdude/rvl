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

#ifndef RVL_ANTENNA_IMPEDANCE_FREQUENCY_HPP
#define RVL_ANTENNA_IMPEDANCE_FREQUENCY_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

namespace rvl {
namespace antenna {

template<typename T>
class impedance_frequency {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    struct antenna_geometry {
        T length_m;              
        T radius_m;              
        T height_above_ground_m; 
        T ground_conductivity;   
        T ground_permittivity;   
    };

    struct impedance_point {
        T frequency_hz;          
        complex_type impedance;  
        T resistance;            
        T reactance;             
        T magnitude;             
        T phase_rad;             
        T vswr;                  
        T return_loss_db;        
    };

    struct bandwidth_analysis {
        T center_frequency_hz;   
        T lower_frequency_hz;    
        T upper_frequency_hz;    
        T bandwidth_hz;          
        T fractional_bandwidth;  
        T q_factor;              
        T minimum_vswr;          
        std::vector<T> resonant_frequencies; 
    };

    static complex_type calculate_king_dipole_impedance(T frequency_hz, 
                                                       const antenna_geometry& geometry) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(geometry.length_m, "Antenna length");
        core::check_positive(geometry.radius_m, "Wire radius");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T half_length = geometry.length_m / T(2.0);
        const T kh = k * half_length;
        const T ln_term = std::log(geometry.length_m / geometry.radius_m);

        const T C1 = T(73.1);  
        const T C2 = T(42.5);  
        const T C3 = T(-21.25); 

        const T sin_kh = std::sin(kh);
        const T cos_kh = std::cos(kh);
        const T si_2kh = sine_integral(T(2.0) * kh);
        const T ci_2kh = cosine_integral(T(2.0) * kh);

        const T R1 = C1 * (T(1.0) - cos_kh) / (sin_kh * sin_kh);
        const T R2 = C2 * (si_2kh - T(2.0) * sine_integral(kh)) / (sin_kh * sin_kh);
        const T R3 = C3 * (T(1.0) - cos_kh) / (sin_kh * sin_kh);

        const T X1 = C1 * (sine_integral(kh) - si_2kh/T(2.0)) / (sin_kh * sin_kh);
        const T X2 = C2 * (T(1.0) - cos_kh - ci_2kh + T(2.0) * cosine_integral(kh)) / (sin_kh * sin_kh);
        const T X3 = C3 * (sine_integral(kh) - si_2kh/T(2.0)) / (sin_kh * sin_kh);

        const T resistance = R1 + R2 + R3;
        const T reactance = (X1 + X2 + X3) - T(120.0) * (ln_term - T(2.0));

        return complex_type(resistance, reactance);
    }

    static complex_type calculate_arbitrary_length_dipole_impedance(T frequency_hz,
                                                                  const antenna_geometry& geometry) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(geometry.length_m, "Antenna length");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T electrical_length = geometry.length_m / wavelength;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T half_length = geometry.length_m / T(2.0);
        const T kh = k * half_length;

        const T R_half_wave = T(73.1);
        const T X_half_wave = constants::antenna<T>::dipole_reactance;

        if (std::abs(sin(kh)) < T(0.01)) {

            const T length_deviation = electrical_length - T(0.5);
            const T resistance = R_half_wave * (T(1.0) + T(0.2) * length_deviation * length_deviation);
            const T reactance = T(200.0) * length_deviation; 

            return complex_type(resistance, reactance);
        }

        const T sin_kh_squared = std::sin(kh) * std::sin(kh);
        const T resistance = T(60.0) * (T(1.0) - std::cos(kh)) / sin_kh_squared;

        const T ln_ratio = std::log(geometry.length_m / geometry.radius_m);
        const T reactance = T(60.0) * (sine_integral(kh) / sin_kh_squared - ln_ratio);

        return complex_type(resistance, reactance);
    }

    static complex_type calculate_monopole_impedance(T frequency_hz,
                                                   const antenna_geometry& geometry) {

        const complex_type dipole_impedance = calculate_arbitrary_length_dipole_impedance(frequency_hz, geometry);
        return dipole_impedance / T(2.0);
    }

    static complex_type calculate_loop_impedance(T frequency_hz, const antenna_geometry& geometry) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(geometry.length_m, "Loop circumference");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T circumference = geometry.length_m;
        const T radius = circumference / (T(2.0) * constants::mathematical<T>::pi);
        const T area = constants::mathematical<T>::pi * radius * radius;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T ka = k * radius;

        if (ka < T(0.1)) {

            const T area_over_lambda_squared = area / (wavelength * wavelength);
            const T radiation_resistance = T(20.0) * constants::mathematical<T>::pi * constants::mathematical<T>::pi * 
                                         area_over_lambda_squared * area_over_lambda_squared;

            const T ln_ratio = std::log(T(8.0) * radius / geometry.radius_m) - T(2.0);
            const T reactance = T(120.0) * ka * ln_ratio;

            return complex_type(radiation_resistance, reactance);
        } else {

            const T resistance = T(120.0) * ka * ka / T(2.0);
            const T reactance = T(120.0) * ka * (T(1.0) - T(2.0)/constants::mathematical<T>::pi);

            return complex_type(resistance, reactance);
        }
    }

    static std::vector<impedance_point> frequency_sweep(const antenna_geometry& geometry,
                                                      T start_frequency_hz,
                                                      T stop_frequency_hz,
                                                      size_t num_points = 100,
                                                      const std::string& antenna_type = "dipole") {
        core::check_positive(num_points, "Number of frequency points");
        core::check_range(start_frequency_hz, T(0), stop_frequency_hz, "Start frequency");

        std::vector<impedance_point> sweep_data;
        sweep_data.reserve(num_points);

        const T frequency_step = (stop_frequency_hz - start_frequency_hz) / T(num_points - 1);

        for (size_t i = 0; i < num_points; ++i) {
            impedance_point point;
            point.frequency_hz = start_frequency_hz + T(i) * frequency_step;

            if (antenna_type == "dipole") {
                point.impedance = calculate_arbitrary_length_dipole_impedance(point.frequency_hz, geometry);
            } else if (antenna_type == "monopole") {
                point.impedance = calculate_monopole_impedance(point.frequency_hz, geometry);
            } else if (antenna_type == "loop") {
                point.impedance = calculate_loop_impedance(point.frequency_hz, geometry);
            } else {
                point.impedance = calculate_arbitrary_length_dipole_impedance(point.frequency_hz, geometry);
            }

            point.resistance = point.impedance.real();
            point.reactance = point.impedance.imag();
            point.magnitude = std::abs(point.impedance);
            point.phase_rad = std::arg(point.impedance);

            const complex_type z0(T(50.0), T(0.0));
            const complex_type gamma = (point.impedance - z0) / (point.impedance + z0);
            const T gamma_magnitude = std::abs(gamma);

            if (gamma_magnitude >= T(1.0)) {
                point.vswr = std::numeric_limits<T>::infinity();
                point.return_loss_db = T(0.0);
            } else {
                point.vswr = (T(1.0) + gamma_magnitude) / (T(1.0) - gamma_magnitude);
                point.return_loss_db = -T(20.0) * std::log10(gamma_magnitude);
            }

            sweep_data.push_back(point);
        }

        return sweep_data;
    }

    static bandwidth_analysis analyze_bandwidth(const std::vector<impedance_point>& sweep_data,
                                              T vswr_threshold = T(2.0),
                                              T reference_impedance = T(50.0)) {
        if (sweep_data.empty()) {
            throw core::invalid_argument_error("Sweep data cannot be empty");
        }

        bandwidth_analysis analysis;

        analysis.resonant_frequencies = find_resonant_frequencies(sweep_data);

        if (analysis.resonant_frequencies.empty()) {

            auto min_vswr_it = std::min_element(sweep_data.begin(), sweep_data.end(),
                [](const impedance_point& a, const impedance_point& b) {
                    return a.vswr < b.vswr;
                });
            analysis.center_frequency_hz = min_vswr_it->frequency_hz;
            analysis.minimum_vswr = min_vswr_it->vswr;
        } else {

            analysis.center_frequency_hz = analysis.resonant_frequencies[0];

            auto center_it = std::min_element(sweep_data.begin(), sweep_data.end(),
                [&](const impedance_point& a, const impedance_point& b) {
                    return std::abs(a.frequency_hz - analysis.center_frequency_hz) < 
                           std::abs(b.frequency_hz - analysis.center_frequency_hz);
                });
            analysis.minimum_vswr = center_it->vswr;
        }

        std::vector<T> valid_frequencies;
        for (const auto& point : sweep_data) {
            if (point.vswr <= vswr_threshold) {
                valid_frequencies.push_back(point.frequency_hz);
            }
        }

        if (valid_frequencies.empty()) {
            analysis.bandwidth_hz = T(0.0);
            analysis.fractional_bandwidth = T(0.0);
            analysis.q_factor = std::numeric_limits<T>::infinity();
        } else {
            analysis.lower_frequency_hz = *std::min_element(valid_frequencies.begin(), valid_frequencies.end());
            analysis.upper_frequency_hz = *std::max_element(valid_frequencies.begin(), valid_frequencies.end());
            analysis.bandwidth_hz = analysis.upper_frequency_hz - analysis.lower_frequency_hz;
            analysis.fractional_bandwidth = analysis.bandwidth_hz / analysis.center_frequency_hz;

            if (analysis.bandwidth_hz > T(0)) {
                analysis.q_factor = analysis.center_frequency_hz / analysis.bandwidth_hz;
            } else {
                analysis.q_factor = std::numeric_limits<T>::infinity();
            }
        }

        return analysis;
    }

    static std::vector<T> find_resonant_frequencies(const std::vector<impedance_point>& sweep_data) {
        std::vector<T> resonances;

        for (size_t i = 1; i < sweep_data.size(); ++i) {
            const T prev_reactance = sweep_data[i-1].reactance;
            const T curr_reactance = sweep_data[i].reactance;

            if ((prev_reactance < T(0) && curr_reactance > T(0)) || 
                (prev_reactance > T(0) && curr_reactance < T(0))) {

                const T prev_freq = sweep_data[i-1].frequency_hz;
                const T curr_freq = sweep_data[i].frequency_hz;

                const T crossing_freq = prev_freq - prev_reactance * (curr_freq - prev_freq) / 
                                      (curr_reactance - prev_reactance);

                resonances.push_back(crossing_freq);
            }
        }

        return resonances;
    }

    static T calculate_antenna_q_factor(const std::vector<impedance_point>& sweep_data,
                                       T resonant_frequency_hz) {

        auto resonant_it = std::min_element(sweep_data.begin(), sweep_data.end(),
            [&](const impedance_point& a, const impedance_point& b) {
                return std::abs(a.frequency_hz - resonant_frequency_hz) < 
                       std::abs(b.frequency_hz - resonant_frequency_hz);
            });

        if (resonant_it == sweep_data.end()) {
            throw core::invalid_argument_error("Resonant frequency not found in sweep data");
        }

        const T R_resonant = resonant_it->resistance;

        std::vector<T> q_frequencies;

        for (const auto& point : sweep_data) {
            if (std::abs(std::abs(point.reactance) - R_resonant) < R_resonant * T(0.1)) {
                q_frequencies.push_back(point.frequency_hz);
            }
        }

        if (q_frequencies.size() >= 2) {
            std::sort(q_frequencies.begin(), q_frequencies.end());
            const T bandwidth_q = q_frequencies.back() - q_frequencies.front();
            return resonant_frequency_hz / bandwidth_q;
        }

        const T delta_f = T(0.01) * resonant_frequency_hz;
        const T slope = estimate_reactance_slope(sweep_data, resonant_frequency_hz, delta_f);

        if (slope > T(0)) {
            return slope * resonant_frequency_hz / (T(2.0) * R_resonant);
        }

        return std::numeric_limits<T>::infinity();
    }

    static T sine_integral(T x) {
        if (std::abs(x) < T(1e-6)) {
            return x * (T(1.0) - x*x/T(18.0));
        }

        if (x < T(0)) {
            return -sine_integral(-x);
        }

        if (x < T(4.0)) {

            T sum = T(0.0);
            T term = x;
            T x_squared = x * x;

            for (int n = 0; n < 20; ++n) {
                sum += term / T(2*n + 1);
                term *= -x_squared / T((2*n + 2) * (2*n + 3));

                if (std::abs(term) < T(1e-12)) break;
            }

            return sum;
        } else {

            return constants::mathematical<T>::half_pi - std::cos(x)/x - std::sin(x)/(x*x);
        }
    }

    static T cosine_integral(T x) {
        const T euler_gamma = T(0.5772156649);

        if (x <= T(0)) {
            throw core::invalid_argument_error("Cosine integral requires positive argument");
        }

        if (x < T(4.0)) {

            T sum = euler_gamma + std::log(x);
            T term = -x*x/T(4.0);

            for (int n = 1; n < 20; ++n) {
                sum += term / T(2*n * (2*n));
                term *= -x*x / T((2*n + 1) * (2*n + 2));

                if (std::abs(term) < T(1e-12)) break;
            }

            return sum;
        } else {

            return std::sin(x)/x - std::cos(x)/(x*x);
        }
    }

private:

    static T estimate_reactance_slope(const std::vector<impedance_point>& sweep_data,
                                    T frequency_hz, T delta_f) {

        auto it = std::min_element(sweep_data.begin(), sweep_data.end(),
            [&](const impedance_point& a, const impedance_point& b) {
                return std::abs(a.frequency_hz - frequency_hz) < 
                       std::abs(b.frequency_hz - frequency_hz);
            });

        if (it == sweep_data.begin() || it == sweep_data.end() - 1) {
            return T(0.0); 
        }

        const T df = (it+1)->frequency_hz - (it-1)->frequency_hz;
        const T dx = (it+1)->reactance - (it-1)->reactance;

        return dx / df;
    }

public:
    static constexpr T typical_wire_radius_mm() { return constants::antenna_design<T>::typical_wire_radius_mm; }
    static constexpr T standard_vswr_threshold() { return constants::antenna_design<T>::standard_vswr_threshold; }
    static constexpr T good_vswr_threshold() { return constants::antenna_design<T>::good_vswr_threshold; }
    static constexpr T excellent_vswr_threshold() { return constants::antenna_design<T>::excellent_vswr_threshold; }
    static constexpr size_t default_sweep_points() { return 200; }
    static constexpr T minimum_q_factor() { return constants::antenna_design<T>::minimum_q_factor; }
    static constexpr T typical_dipole_q() { return constants::antenna_design<T>::typical_dipole_q; }
};

using impedance_frequency_f = impedance_frequency<float>;
using impedance_frequency_d = impedance_frequency<double>;

} 
} 

#endif 