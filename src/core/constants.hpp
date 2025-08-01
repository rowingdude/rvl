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

#ifndef RVL_CORE_CONSTANTS_HPP
#define RVL_CORE_CONSTANTS_HPP

#include <cmath>

namespace rvl {
namespace constants {

template<typename T = double>
struct physical {
    static constexpr T c = T(299792458.0);                    
    static constexpr T mu_0 = T(4.0e-7) * M_PI;              
    static constexpr T epsilon_0 = T(8.854187817e-12);        
    static constexpr T eta_0 = T(376.730313668);              
    static constexpr T k_B = T(1.380649e-23);                 
    static constexpr T q_e = T(1.602176634e-19);              
    static constexpr T m_e = T(9.1093837015e-31);             
    static constexpr T earth_radius = T(6371000.0);           
    static constexpr T standard_temp = T(290.0);              
    static constexpr T standard_pressure = T(101325.0);       
    static constexpr T G = T(6.67430e-11);                    
    static constexpr T M_earth = T(5.972e24);                 
    static constexpr T g_0 = T(9.80665);                      
    static constexpr T N_A = T(6.02214076e23);                
    static constexpr T h = T(6.62607015e-34);                 
    static constexpr T h_bar = h / (T(2.0) * M_PI);           
    static constexpr T sigma_sb = T(5.670374419e-8);          
    static constexpr T R = T(8.314462618);                    
    static constexpr T k_B_per_q = k_B / q_e;                 
};

template<typename T = double>
struct mathematical {
    static constexpr T pi = M_PI;
    static constexpr T two_pi = T(2.0) * M_PI;
    static constexpr T half_pi = M_PI / T(2.0);
    static constexpr T quarter_pi = M_PI / T(4.0);
    static constexpr T inv_pi = T(1.0) / M_PI;
    static constexpr T inv_two_pi = T(1.0) / (T(2.0) * M_PI);
    static constexpr T sqrt_two = M_SQRT2;
    static constexpr T inv_sqrt_two = T(1.0) / M_SQRT2;
    static constexpr T euler = M_E;
    static constexpr T ln_2 = M_LN2;
    static constexpr T ln_10 = M_LN10;
    static constexpr T rad_to_deg = T(180.0) / M_PI;
    static constexpr T deg_to_rad = M_PI / T(180.0);    
    static constexpr T sqrt_pi = T(1.772453850905516);        
    static constexpr T inv_sqrt_pi = T(1.0) / sqrt_pi;        
    static constexpr T e_inv = T(1.0) / M_E;                  
    static constexpr T golden_ratio = T(1.618033988749895);   
};

template<typename T = double>
struct radio {
    static constexpr T db_to_linear_factor = T(0.1) * M_LN10;
    static constexpr T linear_to_db_factor = T(10.0) / M_LN10;
    static constexpr T vswr_perfect = T(1.0);
    static constexpr T reflection_coeff_perfect = T(0.0);
    static constexpr T mismatch_loss_perfect = T(0.0);    
    static constexpr T dbm_to_watts = T(0.001);               
    static constexpr T watts_to_dbm = T(30.0);                
    static constexpr T max_vswr_practical = T(20.0);          
    static constexpr T min_gamma_stable = T(0.01);            
};

template<typename T = double>
struct antenna {
    static constexpr T dipole_half_wave_length_factor = T(0.49);  
    static constexpr T monopole_length_factor = T(0.245);         
    static constexpr T dipole_input_impedance = T(73.0);          
    static constexpr T dipole_reactance = T(42.5);                
    static constexpr T ground_plane_reflection_coeff = T(-1.0);   
    static constexpr T yagi_director_factor = T(0.95);            
    static constexpr T yagi_reflector_factor = T(1.05);           
};

template<typename T = double>
struct propagation {
    static constexpr T atmospheric_refraction = T(0.25);          
    static constexpr T earth_curvature_factor = T(4.0/3.0);       
    static constexpr T standard_atmosphere_temp_lapse = T(-0.0065); 
    static constexpr T dry_air_molecular_weight = T(28.9647e-3);  
    static constexpr T water_vapor_molecular_weight = T(18.0153e-3); 
    static constexpr T max_ionospheric_frequency = T(30e6);       
    static constexpr T min_ionospheric_frequency = T(1.6e6);      
};

template<typename T = double>
struct materials {
    static constexpr T copper_conductivity = T(5.96e7);           
    static constexpr T aluminum_conductivity = T(3.77e7);         
    static constexpr T silver_conductivity = T(6.30e7);           
    static constexpr T gold_conductivity = T(4.10e7);             
    static constexpr T copper_relative_permeability = T(1.0);     
    static constexpr T iron_relative_permeability = T(5000.0);    
    static constexpr T teflon_dielectric = T(2.1);                
    static constexpr T polyethylene_dielectric = T(2.25);         
    static constexpr T plexiglass_dielectric = T(2.6);            
    static constexpr T pcb_fr4_dielectric = T(4.4);               
    static constexpr T water_dielectric = T(80.1);                
    static constexpr T brass_conductivity = T(1.5e7);
    static constexpr T rogers_dielectric = T(2.2);
    static constexpr T alumina_dielectric = T(9.6);
};

using physical_f = physical<float>;
using physical_d = physical<double>;
using math_f = mathematical<float>;
using math_d = mathematical<double>;
using radio_f = radio<float>;
using radio_d = radio<double>;
using antenna_f = antenna<float>;
using antenna_d = antenna<double>;
using propagation_f = propagation<float>;
using propagation_d = propagation<double>;
using materials_f = materials<float>;
using materials_d = materials<double>;

template<typename T = double>
struct transmission_lines {
    static constexpr T standard_50_ohm = T(50.0);
    static constexpr T standard_75_ohm = T(75.0);
    static constexpr T standard_300_ohm = T(300.0);
    static constexpr T standard_600_ohm = T(600.0);
    static constexpr T twin_lead_300_spacing_mm = T(7.62);
    static constexpr T twin_lead_300_diameter_mm = T(1.63);
    static constexpr T twin_lead_450_spacing_mm = T(12.7);
    static constexpr T twin_lead_450_diameter_mm = T(1.63);
    static constexpr T twin_lead_600_spacing_mm = T(25.4);
    static constexpr T twin_lead_600_diameter_mm = T(1.63);
};

template<typename T = double>
struct cable_specs {
    static constexpr T rg58_loss_100mhz_db_per_m = T(0.066);
    static constexpr T rg58_loss_1ghz_db_per_m = T(0.216);
    static constexpr T rg213_loss_100mhz_db_per_m = T(0.033);
    static constexpr T rg213_loss_1ghz_db_per_m = T(0.108);
};

template<typename T = double>
struct waveguides {
    static constexpr T wr90_broad_mm = T(22.86);
    static constexpr T wr90_narrow_mm = T(10.16);
    static constexpr T wr62_broad_mm = T(15.80);
    static constexpr T wr62_narrow_mm = T(7.90);
    static constexpr T wr42_broad_mm = T(10.67);
    static constexpr T wr42_narrow_mm = T(4.32);
    static constexpr T wr28_broad_mm = T(7.11);
    static constexpr T wr28_narrow_mm = T(3.56);
    static constexpr T wr22_broad_mm = T(5.69);
    static constexpr T wr22_narrow_mm = T(2.84);
    static constexpr T wr15_broad_mm = T(3.76);
    static constexpr T wr15_narrow_mm = T(1.88);
    static constexpr T air_breakdown_mv_per_m = T(30.0);
    static constexpr T dry_nitrogen_breakdown_mv_per_m = T(35.0);
};

template<typename T = double>
struct vswr_thresholds {
    static constexpr T perfect_match = T(1.0);
    static constexpr T excellent_match = T(1.2);
    static constexpr T good_match = T(1.5);
    static constexpr T acceptable_amateur = T(2.0);
    static constexpr T acceptable_commercial = T(1.5);
    static constexpr T maximum_practical = T(10.0);
};

template<typename T = double>
struct antenna_design {
    static constexpr T excellent_efficiency_threshold = T(0.9);
    static constexpr T good_efficiency_threshold = T(0.7);
    static constexpr T acceptable_efficiency_threshold = T(0.5);
    static constexpr T poor_efficiency_threshold = T(0.3);
    static constexpr T typical_wire_radius_awg12 = T(0.00103);
    static constexpr T typical_wire_radius_awg14 = T(0.000815);
    static constexpr T typical_wire_radius_awg16 = T(0.000648);
    static constexpr T typical_wire_radius_mm = T(1.0);
    static constexpr T minimum_samples_per_wavelength = T(20.0);
    static constexpr T maximum_useful_length_wavelengths = T(5.0);
    static constexpr T maximum_element_spacing_wavelengths = T(1.0);
    static constexpr T typical_element_spacing_wavelengths = T(0.5);
    static constexpr T minimum_element_spacing_wavelengths = T(0.1);
    static constexpr int maximum_practical_elements = 100;
    static constexpr T typical_sidelobe_level_db = T(-20.0);
    static constexpr T low_sidelobe_level_db = T(-40.0);
    static constexpr T typical_substrate_height_mm = T(1.6);
    static constexpr T standard_vswr_threshold = T(2.0);
    static constexpr T good_vswr_threshold = T(1.5);
    static constexpr T excellent_vswr_threshold = T(1.2);
    static constexpr T minimum_q_factor = T(1.0);
    static constexpr T typical_dipole_q = T(10.0);
    static constexpr T typical_antenna_range_distance_multiplier = T(2.0);
    static constexpr T compact_range_distance_multiplier = T(0.5);
    static constexpr T planar_near_field_scan_distance_multiplier = T(0.1);
    static constexpr T typical_num_turns_helical = T(10.0);
    static constexpr T typical_gain_axial_helical_db = T(12.0);
};

template<typename T = double>
struct rf_components {
    static constexpr T thermal_noise_floor_dbm_per_hz = T(-174.0);
    static constexpr T typical_fade_margin_db = T(10.0);
    static constexpr T typical_implementation_loss_db = T(3.0);
    static constexpr T typical_mixer_ip3_dbm = T(10.0);
    static constexpr T typical_lna_ip3_dbm = T(-5.0);
    static constexpr T typical_power_amp_ip3_dbm = T(30.0);
    static constexpr T typical_vga_ip3_dbm = T(20.0);
    static constexpr T excellent_ip3_dbm = T(40.0);
    static constexpr T good_ip3_dbm = T(20.0);
    static constexpr T fair_ip3_dbm = T(10.0);
    static constexpr T poor_ip3_dbm = T(0.0);
    static constexpr T default_impedance = T(50.0);
};

template<typename T = double>
struct noise_performance {
    static constexpr T excellent_phase_noise_1khz_dbc = T(-120.0);
    static constexpr T good_phase_noise_1khz_dbc = T(-100.0);
    static constexpr T fair_phase_noise_1khz_dbc = T(-80.0);
    static constexpr T poor_phase_noise_1khz_dbc = T(-60.0);
    static constexpr T crystal_oscillator_q = T(1e6);
    static constexpr T ceramic_resonator_q = T(1e3);
    static constexpr T lc_oscillator_q = T(100.0);
    static constexpr T rc_oscillator_q = T(10.0);
    static constexpr T thermal_noise_floor_dbc = T(-174.0);
    static constexpr T typical_flicker_corner_hz = T(1000.0);
};

template<typename T = double>
struct atmospheric {
    static constexpr T oxygen_line_frequency_ghz = T(60.0);
    static constexpr T water_vapor_line_frequency_ghz = T(22.235);
    static constexpr T standard_temperature_k = T(288.15);
    static constexpr T standard_pressure_pa = T(101325.0);
    static constexpr T standard_atmosphere_k_factor = T(1.33);
    static constexpr T minimum_clearance_factor = T(0.6);
    static constexpr T maximum_fresnel_parameter = T(2.4);
    static constexpr T air_permittivity = T(1.0);
};

template<typename T = double>
struct ionospheric {
    static constexpr T typical_foF2_mhz = T(10.0);
    static constexpr T typical_hmF2_km = T(300.0);
    static constexpr T typical_scale_height_km = T(50.0);
    static constexpr T typical_electron_density_f2 = T(1e12);
    static constexpr T typical_critical_frequency_mhz = T(10.0);
    static constexpr T typical_collision_frequency = T(1e3);
    static constexpr T typical_structure_constant = T(1e-50);
    static constexpr T typical_correlation_length_km = T(1.0);
    static constexpr T typical_s4_strong = T(0.8);
    static constexpr T typical_s4_weak = T(0.1);
};

template<typename T = double>
struct tropospheric {
    static constexpr T typical_structure_constant_continental = T(1e-14);
    static constexpr T typical_structure_constant_maritime = T(5e-15);
    static constexpr T typical_structure_constant_strong_turbulence = T(1e-13);
    static constexpr T typical_structure_constant_weak_turbulence = T(1e-15);
    static constexpr T typical_outer_scale_m = T(1000.0);
    static constexpr T minimum_elevation_angle_deg = T(0.1);
    static constexpr T maximum_practical_distance_km = T(800.0);
    static constexpr T minimum_practical_distance_km = T(100.0);
    static constexpr T typical_fade_margin_db = T(30.0);
    static constexpr T typical_surface_duct_height_m = T(50.0);
    static constexpr T typical_elevated_duct_height_m = T(200.0);
    static constexpr T typical_evaporation_duct_height_m = T(15.0);
    static constexpr T strong_ducting_gradient = T(-0.157e-6);
    static constexpr T weak_ducting_gradient = T(-0.079e-6);
    static constexpr T maritime_climate_factor = T(2.0);
    static constexpr T continental_climate_factor = T(0.5);
    static constexpr T tropical_climate_factor = T(3.0);
};

template<typename T = double>
struct ground_reflection {
    static constexpr T typical_urban_permittivity = T(5.0);
    static constexpr T typical_rural_permittivity = T(4.0);
    static constexpr T typical_urban_conductivity = T(0.001);
    static constexpr T typical_rural_conductivity = T(0.001);
    static constexpr T typical_surface_roughness_m = T(0.01);
};

template<typename T = double>
struct fading_models {
    static constexpr T typical_urban_k_factor_db = T(3.0);
    static constexpr T typical_suburban_k_factor_db = T(7.0);
    static constexpr T typical_rural_k_factor_db = T(15.0);
    static constexpr T minimum_k_factor_db = T(-10.0);
    static constexpr T maximum_practical_k_factor_db = T(30.0);
    static constexpr T typical_urban_rms_delay_spread_us = T(1.0);
    static constexpr T typical_doppler_at_walking_speed_hz = T(5.0);
    static constexpr T typical_doppler_at_vehicular_speed_hz = T(100.0);
    static constexpr T rayleigh_distribution_mode = T(0);
};

template<typename T = double>
struct cellular_models {
    static constexpr T typical_base_station_height_m = T(50.0);
    static constexpr T typical_mobile_height_m = T(1.5);
    static constexpr T typical_cell_radius_urban_km = T(2.0);
    static constexpr T typical_cell_radius_suburban_km = T(5.0);
    static constexpr T typical_cell_radius_rural_km = T(15.0);
};

template<typename T = double>
struct polarization {
    static constexpr T perfect_match_loss_db = T(0.0);
    static constexpr T orthogonal_linear_loss_db = T(100.0);
    static constexpr T linear_to_circular_loss_db = T(3.0);
    static constexpr T opposite_circular_loss_db = T(100.0);
    static constexpr T typical_multipath_depolarization_db = T(3.0);
    static constexpr T rain_depolarization_db_per_degree = T(0.1);
};

template<typename T = double>
struct building_materials {
    static constexpr T typical_concrete_permittivity = T(5.5);
    static constexpr T typical_glass_permittivity = T(6.0);
    static constexpr T typical_metal_conductivity = T(10000000.0);
};

template<typename T = double>
struct computational_limits {
    static constexpr int maximum_ray_bounces = 5;
    static constexpr int cache_line_size = 64;
    static constexpr int simd_alignment = 64;
};

template<typename T = double>
struct loss_tangents {
    static constexpr T air = T(0.0);
    static constexpr T fr4 = T(0.02);
    static constexpr T teflon = T(0.0001);
    static constexpr T ceramic = T(0.001);
};

using transmission_lines_f = transmission_lines<float>;
using transmission_lines_d = transmission_lines<double>;
using cable_specs_f = cable_specs<float>;
using cable_specs_d = cable_specs<double>;
using waveguides_f = waveguides<float>;
using waveguides_d = waveguides<double>;
using vswr_thresholds_f = vswr_thresholds<float>;
using vswr_thresholds_d = vswr_thresholds<double>;
using antenna_design_f = antenna_design<float>;
using antenna_design_d = antenna_design<double>;
using rf_components_f = rf_components<float>;
using rf_components_d = rf_components<double>;
using noise_performance_f = noise_performance<float>;
using noise_performance_d = noise_performance<double>;
using atmospheric_f = atmospheric<float>;
using atmospheric_d = atmospheric<double>;
using ionospheric_f = ionospheric<float>;
using ionospheric_d = ionospheric<double>;
using tropospheric_f = tropospheric<float>;
using tropospheric_d = tropospheric<double>;
using ground_reflection_f = ground_reflection<float>;
using ground_reflection_d = ground_reflection<double>;
using fading_models_f = fading_models<float>;
using fading_models_d = fading_models<double>;
using cellular_models_f = cellular_models<float>;
using cellular_models_d = cellular_models<double>;
using polarization_f = polarization<float>;
using polarization_d = polarization<double>;
using building_materials_f = building_materials<float>;
using building_materials_d = building_materials<double>;
using computational_limits_f = computational_limits<float>;
using computational_limits_d = computational_limits<double>;
using loss_tangents_f = loss_tangents<float>;
using loss_tangents_d = loss_tangents<double>;

} 
} 

#endif 