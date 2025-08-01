#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include "../src/propagation/tropospheric_ducting.hpp"
#include "../src/propagation/tropospheric_scatter.hpp"
#include "../src/propagation/atmospheric_attenuation.hpp"
#include "../src/propagation/ground_wave.hpp"
#include "../src/propagation/terrain_diffraction.hpp"

using namespace rvl::propagation;

void test_tropospheric_ducting() {
    std::cout << "Testing tropospheric ducting calculations..." << std::endl;

    const double sea_temp = 288.15;
    const double air_temp = 283.15;
    const double wind_speed = 10.0;
    const double humidity = 80.0;

    const double duct_height = tropospheric_ducting_d::calculate_evaporation_duct_height(
        sea_temp, air_temp, wind_speed, humidity);

    assert(duct_height > 0.0 && duct_height < 100.0);

    const double m_gradient = -200.0;
    assert(tropospheric_ducting_d::is_ducting_condition(m_gradient));

    const double duct_strength = tropospheric_ducting_d::calculate_duct_strength(m_gradient);
    assert(duct_strength > 0.0);

    const double ducting_range = tropospheric_ducting_d::calculate_ducting_range_km(
        3.0, duct_height, duct_strength);

    std::cout << "  Evaporation duct height: " << duct_height << " m" << std::endl;
    std::cout << "  Duct strength: " << duct_strength << std::endl;
    std::cout << "  Ducting range: " << ducting_range << " km" << std::endl;
    std::cout << "✓ Tropospheric ducting tests passed" << std::endl;
}

void test_tropospheric_scatter() {
    std::cout << "Testing tropospheric scatter calculations..." << std::endl;

    tropospheric_scatter_d::scatter_geometry geometry;
    geometry.transmitter_height_m = 30.0;
    geometry.receiver_height_m = 30.0;
    geometry.path_distance_km = 300.0;
    geometry.common_volume_height_m = tropospheric_scatter_d::calculate_common_volume_height(
        geometry.transmitter_height_m, geometry.receiver_height_m, geometry.path_distance_km);
    geometry.scattering_angle_rad = tropospheric_scatter_d::calculate_scattering_angle(
        geometry.transmitter_height_m, geometry.receiver_height_m, 
        geometry.path_distance_km, geometry.common_volume_height_m);

    assert(geometry.common_volume_height_m > 1000.0);
    assert(geometry.scattering_angle_rad > 0.0);

    const double path_loss = tropospheric_scatter_d::calculate_troposcatter_path_loss_db(
        1e9, geometry, 40.0, 40.0);

    assert(path_loss > 100.0 && path_loss < 300.0);

    const double availability = tropospheric_scatter_d::calculate_troposcatter_availability_percent(
        1.0, 300.0);

    assert(availability > 90.0 && availability <= 100.0);

    std::cout << "  Common volume height: " << geometry.common_volume_height_m << " m" << std::endl;
    std::cout << "  Scattering angle: " << geometry.scattering_angle_rad * 180.0 / M_PI << "°" << std::endl;
    std::cout << "  Path loss: " << path_loss << " dB" << std::endl;
    std::cout << "  Availability: " << availability << "%" << std::endl;
    std::cout << "✓ Tropospheric scatter tests passed" << std::endl;
}

void test_atmospheric_attenuation() {
    std::cout << "Testing atmospheric attenuation calculations..." << std::endl;

    atmospheric_attenuation_d::atmospheric_conditions conditions;
    conditions.temperature_k = 288.15;
    conditions.pressure_pa = 101325.0;
    conditions.humidity_percent = 50.0;
    conditions.water_vapor_density_gm3 = atmospheric_attenuation_d::calculate_water_vapor_density_from_humidity(
        conditions.temperature_k, conditions.pressure_pa, conditions.humidity_percent);

    const double oxygen_atten = atmospheric_attenuation_d::calculate_oxygen_specific_attenuation_db_km(
        10.0, conditions.pressure_pa, conditions.temperature_k);

    const double water_vapor_atten = atmospheric_attenuation_d::calculate_water_vapor_specific_attenuation_db_km(
        22.235, conditions.water_vapor_density_gm3, conditions.pressure_pa, conditions.temperature_k);

    const double total_atten = atmospheric_attenuation_d::calculate_total_gaseous_attenuation_db_km(
        10.0, conditions);

    assert(oxygen_atten >= 0.0);
    assert(water_vapor_atten >= 0.0);
    assert(total_atten >= oxygen_atten);

    const double rain_atten = atmospheric_attenuation_d::calculate_rain_specific_attenuation_db_km(
        10.0, 10.0);

    assert(rain_atten > 0.0);

    std::cout << "  Oxygen attenuation (10 GHz): " << oxygen_atten << " dB/km" << std::endl;
    std::cout << "  Water vapor attenuation (22.235 GHz): " << water_vapor_atten << " dB/km" << std::endl;
    std::cout << "  Total gaseous attenuation: " << total_atten << " dB/km" << std::endl;
    std::cout << "  Rain attenuation (10mm/h): " << rain_atten << " dB/km" << std::endl;
    std::cout << "✓ Atmospheric attenuation tests passed" << std::endl;
}

void test_ground_wave() {
    std::cout << "Testing ground wave calculations..." << std::endl;

    ground_wave_d::ground_parameters ground = ground_wave_d::create_ground_parameters("medium_dry_ground");

    assert(ground.relative_permittivity > 1.0);
    assert(ground.conductivity_s_per_m > 0.0);

    const double path_loss = ground_wave_d::calculate_ground_wave_path_loss_db(
        3e6, 10000.0, ground);

    assert(path_loss > 40.0 && path_loss < 120.0);

    const double field_strength = ground_wave_d::calculate_ground_wave_field_strength_db(
        1000.0, 3e6, 10000.0, ground, 30.0, 2.0);

    const double surface_wave_range = ground_wave_d::calculate_surface_wave_range_km(
        3e6, 40.0, -100.0, ground);

    assert(surface_wave_range > 0.0);

    std::cout << "  Ground conductivity: " << ground.conductivity_s_per_m << " S/m" << std::endl;
    std::cout << "  Ground permittivity: " << ground.relative_permittivity << std::endl;
    std::cout << "  Path loss (3 MHz, 10 km): " << path_loss << " dB" << std::endl;
    std::cout << "  Field strength: " << field_strength << " dBμV/m" << std::endl;
    std::cout << "  Surface wave range: " << surface_wave_range << " km" << std::endl;
    std::cout << "✓ Ground wave tests passed" << std::endl;
}

void test_terrain_diffraction() {
    std::cout << "Testing terrain diffraction calculations..." << std::endl;

    const double obstacle_height = 100.0;
    const double tx_distance = 5000.0;
    const double rx_distance = 5000.0;
    const double frequency = 100e6;
    const double tx_height = 30.0;
    const double rx_height = 10.0;

    const double fresnel_param = terrain_diffraction_d::calculate_fresnel_parameter(
        obstacle_height, tx_distance, rx_distance, frequency, tx_height, rx_height);

    const double diffraction_loss = terrain_diffraction_d::calculate_knife_edge_diffraction_loss_db(fresnel_param);

    assert(diffraction_loss >= 0.0);

    const double fresnel_radius = terrain_diffraction_d::calculate_fresnel_zone_radius(
        tx_distance, rx_distance, frequency);

    assert(fresnel_radius > 0.0);

    auto result = terrain_diffraction_d::analyze_single_knife_edge(
        obstacle_height, tx_distance, rx_distance, frequency, tx_height, rx_height);

    assert(result.diffraction_loss_db >= 0.0);

    const double smooth_earth_loss = terrain_diffraction_d::calculate_smooth_earth_diffraction_db(
        50000.0, frequency, tx_height, rx_height);

    std::cout << "  Fresnel parameter: " << fresnel_param << std::endl;
    std::cout << "  Knife-edge loss: " << diffraction_loss << " dB" << std::endl;
    std::cout << "  Fresnel zone radius: " << fresnel_radius << " m" << std::endl;
    std::cout << "  Line of sight: " << (result.line_of_sight ? "Yes" : "No") << std::endl;
    std::cout << "  Smooth earth loss: " << smooth_earth_loss << " dB" << std::endl;
    std::cout << "✓ Terrain diffraction tests passed" << std::endl;
}

void test_batch_operations() {
    std::cout << "Testing batch operations..." << std::endl;

    const size_t n = 5;

    tropospheric_ducting_d::vector_type sea_temps = {288.15, 290.15, 285.15, 292.15, 287.15};
    tropospheric_ducting_d::vector_type air_temps = {283.15, 285.15, 280.15, 287.15, 282.15};
    tropospheric_ducting_d::vector_type wind_speeds = {5.0, 10.0, 15.0, 8.0, 12.0};
    tropospheric_ducting_d::vector_type humidities = {70.0, 80.0, 60.0, 85.0, 75.0};
    tropospheric_ducting_d::vector_type duct_heights(n);

    tropospheric_ducting_d::calculate_duct_heights_batch(
        sea_temps, air_temps, wind_speeds, humidities, duct_heights);

    for (size_t i = 0; i < n; ++i) {
        assert(duct_heights[i] > 0.0 && duct_heights[i] < 100.0);
    }

    atmospheric_attenuation_d::vector_type frequencies = {1.0, 5.0, 10.0, 20.0, 60.0};
    atmospheric_attenuation_d::vector_type temperatures = {288.15, 290.15, 285.15, 292.15, 280.15};
    atmospheric_attenuation_d::vector_type pressures = {101325.0, 95000.0, 102000.0, 98000.0, 105000.0};
    atmospheric_attenuation_d::vector_type humidity_vals = {50.0, 60.0, 40.0, 70.0, 30.0};
    atmospheric_attenuation_d::vector_type attenuations(n);

    atmospheric_attenuation_d::calculate_gaseous_attenuation_batch(
        frequencies, temperatures, pressures, humidity_vals, attenuations);

    for (size_t i = 0; i < n; ++i) {
        assert(attenuations[i] >= 0.0);
    }

    std::cout << "  Processed " << n << " duct height calculations" << std::endl;
    std::cout << "  Processed " << n << " atmospheric attenuation calculations" << std::endl;
    std::cout << "✓ Batch operation tests passed" << std::endl;
}

int main() {
    std::cout << "=== RadioVectorLib Tropospheric/Ground Wave Equations Test ===" << std::endl;

    try {
        test_tropospheric_ducting();
        test_tropospheric_scatter();
        test_atmospheric_attenuation();
        test_ground_wave();
        test_terrain_diffraction();
        test_batch_operations();

        std::cout << std::endl;
        std::cout << "🎉 All tropospheric/ground wave equation tests passed!" << std::endl;
        std::cout << "Implemented equations:" << std::endl;
        std::cout << "  • Tropospheric ducting and evaporation ducts" << std::endl;
        std::cout << "  • Tropospheric scatter propagation" << std::endl;
        std::cout << "  • Atmospheric gas attenuation (O2, H2O)" << std::endl;
        std::cout << "  • Rain and cloud attenuation" << std::endl;
        std::cout << "  • Ground wave propagation (Sommerfeld solution)" << std::endl;
        std::cout << "  • Terrain diffraction (knife-edge, multiple obstacles)" << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "❌ Test failed: " << e.what() << std::endl;
        return 1;
    }
}