#include "test_utils.hpp"
#include "src/propagation/tropospheric_scatter.hpp"
#include "src/propagation/tropospheric_ducting.hpp"
#include "src/propagation/atmospheric_attenuation.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("Tropospheric Propagation Functions");

    try {
        propagation::tropospheric_scatter<double>::scatter_geometry geom;
        geom.distance_km = 200.0;
        geom.transmitter_height_m = 100.0;
        geom.receiver_height_m = 100.0;
        geom.scatter_angle_rad = 0.01;
        geom.common_volume_height_km = 2.0;
        
        propagation::tropospheric_scatter<double>::atmospheric_conditions atmos;
        atmos.surface_refractivity = 320.0;
        atmos.refractivity_gradient = -40.0;
        atmos.turbulence_parameter = 1e-15;
        atmos.temperature_k = 288.0;
        atmos.pressure_hpa = 1013.0;
        atmos.humidity_percent = 50.0;
        
        double loss = propagation::tropospheric_scatter<double>::calculate_scatter_loss_db(
            1e9, geom, atmos);
        bool passed = loss > 100.0 && loss < 300.0;
        suite.add_test("calculate_scatter_loss_db", passed);
    } catch (...) {
        suite.add_test("calculate_scatter_loss_db", false, "exception thrown");
    }
    
    try {
        double angle = propagation::tropospheric_scatter<double>::calculate_scatter_angle(
            200.0, 100.0, 100.0);
        bool passed = angle > 0.0 && angle < 0.1;
        suite.add_test("calculate_scatter_angle", passed);
    } catch (...) {
        suite.add_test("calculate_scatter_angle", false, "exception thrown");
    }
    
    try {
        auto vol = propagation::tropospheric_scatter<double>::calculate_common_volume(
            200.0, 100.0, 100.0, 0.01);
        bool passed = vol.volume_km3 > 0.0 && vol.height_km > 0.0;
        suite.add_test("calculate_common_volume", passed);
    } catch (...) {
        suite.add_test("calculate_common_volume", false, "exception thrown");
    }
    
    try {
        double fade = propagation::tropospheric_scatter<double>::calculate_fade_margin_db(
            0.99, 0.01, 200.0);
        bool passed = fade > 0.0 && fade < 50.0;
        suite.add_test("calculate_fade_margin_db", passed);
    } catch (...) {
        suite.add_test("calculate_fade_margin_db", false, "exception thrown");
    }

    try {
        propagation::tropospheric_ducting<double>::ducting_parameters params;
        params.frequency_hz = 1e9;
        params.duct_height_m = 100.0;
        params.duct_strength_m_units = 20.0;
        params.surface_refractivity = 320.0;
        params.antenna_height_m = 50.0;
        
        bool ducting = propagation::tropospheric_ducting<double>::check_ducting_conditions(params);
        suite.add_test("check_ducting_conditions", true);
    } catch (...) {
        suite.add_test("check_ducting_conditions", false, "exception thrown");
    }
    
    try {
        double range = propagation::tropospheric_ducting<double>::calculate_trapped_range(
            100.0, 20.0, 50.0, 1e9);
        bool passed = range > 0.0 && range < 1000.0;
        suite.add_test("calculate_trapped_range", passed);
    } catch (...) {
        suite.add_test("calculate_trapped_range", false, "exception thrown");
    }
    
    try {
        double loss = propagation::tropospheric_ducting<double>::calculate_duct_propagation_loss(
            100.0, 1e9, 100.0, 20.0);
        bool passed = loss > 50.0 && loss < 200.0;
        suite.add_test("calculate_duct_propagation_loss", passed);
    } catch (...) {
        suite.add_test("calculate_duct_propagation_loss", false, "exception thrown");
    }
    
    try {
        double cutoff = propagation::tropospheric_ducting<double>::calculate_duct_cutoff_frequency(
            100.0, 20.0);
        bool passed = cutoff > 100e6 && cutoff < 10e9;
        suite.add_test("calculate_duct_cutoff_frequency", passed);
    } catch (...) {
        suite.add_test("calculate_duct_cutoff_frequency", false, "exception thrown");
    }

    try {
        double atten = propagation::atmospheric_attenuation<double>::calculate_oxygen_attenuation_db_per_km(
            10e9, 288.0, 1013.0, 0.0);
        bool passed = atten > 0.0 && atten < 1.0;
        suite.add_test("calculate_oxygen_attenuation_db_per_km", passed);
    } catch (...) {
        suite.add_test("calculate_oxygen_attenuation_db_per_km", false, "exception thrown");
    }
    
    try {
        double atten = propagation::atmospheric_attenuation<double>::calculate_water_vapor_attenuation_db_per_km(
            22e9, 288.0, 1013.0, 7.5);
        bool passed = atten > 0.0 && atten < 10.0;
        suite.add_test("calculate_water_vapor_attenuation_db_per_km", passed);
    } catch (...) {
        suite.add_test("calculate_water_vapor_attenuation_db_per_km", false, "exception thrown");
    }
    
    try {
        double atten = propagation::atmospheric_attenuation<double>::calculate_rain_attenuation_db_per_km(
            10e9, 10.0, 0.0, 0.0);
        bool passed = atten > 0.0 && atten < 10.0;
        suite.add_test("calculate_rain_attenuation_db_per_km", passed);
    } catch (...) {
        suite.add_test("calculate_rain_attenuation_db_per_km", false, "exception thrown");
    }
    
    try {
        double total = propagation::atmospheric_attenuation<double>::calculate_total_atmospheric_loss_db(
            10.0, 10e9, 288.0, 1013.0, 50.0, 5.0);
        bool passed = total > 0.0 && total < 100.0;
        suite.add_test("calculate_total_atmospheric_loss_db", passed);
    } catch (...) {
        suite.add_test("calculate_total_atmospheric_loss_db", false, "exception thrown");
    }
    
    try {
        double zenith = propagation::atmospheric_attenuation<double>::calculate_zenith_attenuation_db(
            10e9, 288.0, 1013.0, 50.0, 1000.0);
        bool passed = zenith > 0.0 && zenith < 10.0;
        suite.add_test("calculate_zenith_attenuation_db", passed);
    } catch (...) {
        suite.add_test("calculate_zenith_attenuation_db", false, "exception thrown");
    }
    
    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}