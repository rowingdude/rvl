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

#ifndef RVL_ANTENNA_RADIATION_EFFICIENCY_HPP
#define RVL_ANTENNA_RADIATION_EFFICIENCY_HPP

#include "../core/constants.hpp"
#include "../core/error.hpp"
#include "../core/memory_allocator.hpp"
#include <cmath>
#include <complex>
#include <vector>

namespace rvl {
namespace antenna {

template<typename T>
class radiation_efficiency {
public:
    using value_type = T;
    using complex_type = std::complex<T>;
    using vector_type = core::memory::simd_vector<T>;
    using complex_vector_type = core::memory::simd_vector<complex_type>;

    struct conductor_properties {
        T conductivity_s_per_m;     
        T relative_permeability;    
        T diameter_m;               
        T surface_roughness_factor; 
    };

    struct dielectric_properties {
        T relative_permittivity;    
        T loss_tangent;             
        T thickness_m;              
        T volume_m3;                
    };

    struct efficiency_analysis {
        T radiation_efficiency;     
        T mismatch_efficiency;      
        T total_efficiency;         
        T radiation_resistance;     
        T loss_resistance;          
        T ohmic_loss_db;           
        T dielectric_loss_db;      
        T mismatch_loss_db;        
        T gain_reduction_db;       
    };

    struct antenna_parameters {
        T frequency_hz;             
        T physical_length_m;        
        T effective_length_m;       
        T radiation_resistance;     
        T input_impedance_magnitude; 
        complex_type input_impedance; 
    };

    static T calculate_conductor_loss_resistance(T frequency_hz,
                                               T wire_length_m,
                                               const conductor_properties& conductor) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(wire_length_m, "Wire length");
        core::check_positive(conductor.conductivity_s_per_m, "Conductivity");
        core::check_positive(conductor.diameter_m, "Wire diameter");

        const T mu_0 = constants::physical<T>::mu_0;
        const T mu_r = conductor.relative_permeability;
        const T sigma = conductor.conductivity_s_per_m;
        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;

        const T skin_depth = std::sqrt(T(2.0) / (omega * mu_0 * mu_r * sigma));

        const T surface_resistance = T(1.0) / (sigma * skin_depth);

        const T radius = conductor.diameter_m / T(2.0);

        T ac_resistance;
        if (radius > T(2.0) * skin_depth) {

            ac_resistance = surface_resistance / (T(2.0) * constants::mathematical<T>::pi * radius);
        } else {

            const T ratio = radius / skin_depth;
            const T correction_factor = ratio / T(2.0) * (T(1.0) + ratio * ratio / T(12.0));
            ac_resistance = surface_resistance * correction_factor / (constants::mathematical<T>::pi * radius * radius);
        }

        ac_resistance *= conductor.surface_roughness_factor;

        return ac_resistance * wire_length_m;
    }

    static T calculate_dielectric_loss_resistance(T frequency_hz,
                                                const dielectric_properties& dielectric,
                                                const antenna_parameters& antenna) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(dielectric.relative_permittivity, "Relative permittivity");
        core::check_non_negative(dielectric.loss_tangent, "Loss tangent");

        if (dielectric.loss_tangent <= T(0)) {
            return T(0.0); 
        }

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
        const T epsilon_0 = constants::physical<T>::epsilon_0;

        const T field_energy_density = T(0.5) * epsilon_0 * dielectric.relative_permittivity;

        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;
        const T power_loss_density = omega * epsilon_0 * dielectric.relative_permittivity * 
                                   dielectric.loss_tangent;

        const T average_field_squared = std::pow(constants::physical<T>::eta_0 * 
                                               antenna.input_impedance_magnitude, T(2.0));

        const T total_power_loss = power_loss_density * average_field_squared * dielectric.volume_m3;

        const T input_current_squared = T(1.0); 
        return total_power_loss / input_current_squared;
    }

    static T calculate_ground_loss_resistance(T frequency_hz,
                                            T antenna_height_m,
                                            const conductor_properties& ground) {
        core::check_positive(frequency_hz, "Frequency");
        core::check_positive(antenna_height_m, "Antenna height");

        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T mu_0 = constants::physical<T>::mu_0;
        const T omega = T(2.0) * constants::mathematical<T>::pi * frequency_hz;

        const T skin_depth = std::sqrt(T(2.0) / (omega * mu_0 * ground.conductivity_s_per_m));

        const complex_type z_surface = complex_type(T(1.0), T(1.0)) / 
                                     (ground.conductivity_s_per_m * skin_depth);

        const T height_factor = std::exp(-T(2.0) * antenna_height_m / wavelength);
        const T ground_resistance = std::real(z_surface) * height_factor;

        return ground_resistance;
    }

    static T calculate_theoretical_radiation_resistance(const std::string& antenna_type,
                                                      T frequency_hz,
                                                      const antenna_parameters& params) {
        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T electrical_length = params.physical_length_m / wavelength;

        if (antenna_type == "dipole") {
            if (std::abs(electrical_length - T(0.5)) < T(0.1)) {

                return T(73.1);
            } else {

                const T k = T(2.0) * constants::mathematical<T>::pi / wavelength;
                const T kh = k * params.physical_length_m / T(2.0);
                const T sin_kh = std::sin(kh);

                if (std::abs(sin_kh) < T(0.01)) {
                    return T(1000.0); 
                }

                return T(60.0) * (T(1.0) - std::cos(kh)) / (sin_kh * sin_kh);
            }
        } else if (antenna_type == "monopole") {

            return calculate_theoretical_radiation_resistance("dipole", frequency_hz, params) / T(2.0);
        } else if (antenna_type == "loop") {
            const T circumference = params.physical_length_m;
            const T area = std::pow(circumference / (T(2.0) * constants::mathematical<T>::pi), T(2.0)) * 
                          constants::mathematical<T>::pi;
            const T area_over_lambda_squared = area / (wavelength * wavelength);

            return T(20.0) * std::pow(constants::mathematical<T>::pi * area_over_lambda_squared, T(2.0));
        } else if (antenna_type == "helical") {

            const T turns = params.physical_length_m / wavelength; 
            return T(140.0) * turns; 
        }

        return T(73.1);
    }

    static T calculate_mismatch_efficiency(const complex_type& reflection_coefficient) {
        const T gamma_magnitude_squared = std::norm(reflection_coefficient);
        return T(1.0) - gamma_magnitude_squared;
    }

    static T calculate_mismatch_efficiency_from_vswr(T vswr) {
        core::check_range(vswr, T(1.0), std::numeric_limits<T>::max(), "VSWR");

        const T gamma = (vswr - T(1.0)) / (vswr + T(1.0));
        return T(1.0) - gamma * gamma;
    }

    static efficiency_analysis analyze_antenna_efficiency(const antenna_parameters& antenna,
                                                        const conductor_properties& conductor,
                                                        const dielectric_properties& dielectric,
                                                        const std::string& antenna_type = "dipole",
                                                        T reference_impedance = T(50.0)) {
        efficiency_analysis analysis;

        analysis.radiation_resistance = antenna.radiation_resistance > T(0) ? 
            antenna.radiation_resistance : 
            calculate_theoretical_radiation_resistance(antenna_type, antenna.frequency_hz, antenna);

        const T conductor_loss = calculate_conductor_loss_resistance(
            antenna.frequency_hz, antenna.physical_length_m, conductor);

        const T dielectric_loss = calculate_dielectric_loss_resistance(
            antenna.frequency_hz, dielectric, antenna);

        analysis.loss_resistance = conductor_loss + dielectric_loss;

        const T total_resistance = analysis.radiation_resistance + analysis.loss_resistance;
        if (total_resistance > T(0)) {
            analysis.radiation_efficiency = analysis.radiation_resistance / total_resistance;
        } else {
            analysis.radiation_efficiency = T(1.0);
        }

        const complex_type gamma = (antenna.input_impedance - complex_type(reference_impedance, T(0.0))) /
                                 (antenna.input_impedance + complex_type(reference_impedance, T(0.0)));
        analysis.mismatch_efficiency = calculate_mismatch_efficiency(gamma);

        analysis.total_efficiency = analysis.radiation_efficiency * analysis.mismatch_efficiency;

        analysis.ohmic_loss_db = -T(10.0) * std::log10(analysis.radiation_efficiency);
        analysis.dielectric_loss_db = dielectric_loss > T(0) ? 
            -T(10.0) * std::log10(T(1.0) - dielectric_loss / total_resistance) : T(0.0);
        analysis.mismatch_loss_db = -T(10.0) * std::log10(analysis.mismatch_efficiency);
        analysis.gain_reduction_db = -T(10.0) * std::log10(analysis.total_efficiency);

        return analysis;
    }

    static std::vector<efficiency_analysis> efficiency_frequency_sweep(
        const antenna_parameters& base_antenna,
        const conductor_properties& conductor,
        const dielectric_properties& dielectric,
        const vector_type& frequencies,
        const std::string& antenna_type = "dipole") {

        std::vector<efficiency_analysis> results;
        results.reserve(frequencies.size());

        for (const auto& freq : frequencies) {
            antenna_parameters antenna = base_antenna;
            antenna.frequency_hz = freq;

            const T wavelength = constants::physical<T>::c / freq;
            antenna.effective_length_m = antenna.physical_length_m; 

            results.push_back(analyze_antenna_efficiency(antenna, conductor, dielectric, antenna_type));
        }

        return results;
    }

    static std::vector<efficiency_analysis> compare_conductor_materials(
        const antenna_parameters& antenna,
        const std::vector<conductor_properties>& materials,
        const dielectric_properties& dielectric,
        const std::string& antenna_type = "dipole") {

        std::vector<efficiency_analysis> results;
        results.reserve(materials.size());

        for (const auto& material : materials) {
            results.push_back(analyze_antenna_efficiency(antenna, material, dielectric, antenna_type));
        }

        return results;
    }

    static conductor_properties create_conductor_material(const std::string& material_name) {
        conductor_properties props;
        props.relative_permeability = T(1.0);
        props.surface_roughness_factor = T(1.0);

        if (material_name == "copper") {
            props.conductivity_s_per_m = constants::materials<T>::copper_conductivity;
            props.diameter_m = T(0.001); 
        } else if (material_name == "aluminum") {
            props.conductivity_s_per_m = constants::materials<T>::aluminum_conductivity;
            props.diameter_m = T(0.001);
        } else if (material_name == "silver") {
            props.conductivity_s_per_m = constants::materials<T>::silver_conductivity;
            props.diameter_m = T(0.001);
        } else if (material_name == "brass") {
            props.conductivity_s_per_m = T(1.5e7); 
            props.diameter_m = T(0.001);
        } else if (material_name == "steel") {
            props.conductivity_s_per_m = T(1.0e7); 
            props.relative_permeability = T(100.0); 
            props.diameter_m = T(0.001);
        } else {

            props.conductivity_s_per_m = constants::materials<T>::copper_conductivity;
            props.diameter_m = T(0.001);
        }

        return props;
    }

    static dielectric_properties create_dielectric_material(const std::string& material_name) {
        dielectric_properties props;
        props.thickness_m = T(0.001); 
        props.volume_m3 = T(0.001);   

        if (material_name == "air") {
            props.relative_permittivity = T(1.0);
            props.loss_tangent = T(0.0);
        } else if (material_name == "fr4") {
            props.relative_permittivity = T(4.3);
            props.loss_tangent = T(0.02);
        } else if (material_name == "teflon") {
            props.relative_permittivity = constants::materials<T>::teflon_dielectric;
            props.loss_tangent = T(0.0001);
        } else if (material_name == "polyethylene") {
            props.relative_permittivity = constants::materials<T>::polyethylene_dielectric;
            props.loss_tangent = T(0.0005);
        } else if (material_name == "ceramic") {
            props.relative_permittivity = T(9.8);
            props.loss_tangent = T(0.001);
        } else {

            props.relative_permittivity = T(1.0);
            props.loss_tangent = T(0.0);
        }

        return props;
    }

    static T calculate_aperture_efficiency(T physical_aperture_area,
                                         T effective_aperture_area,
                                         T illumination_efficiency = T(0.55),
                                         T spillover_efficiency = T(0.9),
                                         T phase_efficiency = T(0.95)) {
        core::check_positive(physical_aperture_area, "Physical aperture area");
        core::check_positive(effective_aperture_area, "Effective aperture area");

        const T geometric_efficiency = effective_aperture_area / physical_aperture_area;

        return geometric_efficiency * illumination_efficiency * spillover_efficiency * phase_efficiency;
    }

    static T calculate_finite_ground_efficiency(T antenna_height_m,
                                              T ground_plane_radius_m,
                                              T frequency_hz) {
        const T wavelength = constants::physical<T>::c / frequency_hz;
        const T normalized_radius = ground_plane_radius_m / wavelength;
        const T normalized_height = antenna_height_m / wavelength;

        const T efficiency_factor = T(1.0) - std::exp(-T(2.0) * normalized_radius) * 
                                  (T(1.0) + normalized_height);

        return std::max(efficiency_factor, T(0.1)); 
    }

    static void calculate_efficiency_batch(const std::vector<antenna_parameters>& antennas,
                                         const conductor_properties& conductor,
                                         const dielectric_properties& dielectric,
                                         std::vector<efficiency_analysis>& results,
                                         const std::string& antenna_type = "dipole") {
        results.clear();
        results.reserve(antennas.size());

        for (const auto& antenna : antennas) {
            results.push_back(analyze_antenna_efficiency(antenna, conductor, dielectric, antenna_type));
        }
    }

    static constexpr T excellent_efficiency_threshold() { return constants::antenna_design<T>::excellent_efficiency_threshold; }
    static constexpr T good_efficiency_threshold() { return constants::antenna_design<T>::good_efficiency_threshold; }
    static constexpr T acceptable_efficiency_threshold() { return constants::antenna_design<T>::acceptable_efficiency_threshold; }
    static constexpr T poor_efficiency_threshold() { return constants::antenna_design<T>::poor_efficiency_threshold; }

    static constexpr T copper_conductivity() { return constants::materials<T>::copper_conductivity; }
    static constexpr T aluminum_conductivity() { return constants::materials<T>::aluminum_conductivity; }
    static constexpr T silver_conductivity() { return constants::materials<T>::silver_conductivity; }
    static constexpr T brass_conductivity() { return constants::materials<T>::brass_conductivity; }

    static constexpr T air_loss_tangent() { return constants::loss_tangents<T>::air; }
    static constexpr T fr4_loss_tangent() { return T(0.02); }
    static constexpr T teflon_loss_tangent() { return T(0.0001); }
    static constexpr T ceramic_loss_tangent() { return T(0.001); }
};

using radiation_efficiency_f = radiation_efficiency<float>;
using radiation_efficiency_d = radiation_efficiency<double>;

} 
} 

#endif 