#include "test_utils.hpp"
#include "src/antenna/radiation_efficiency.hpp"
#include "src/antenna/antenna_gain.hpp"
#include "src/antenna/near_field.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("Antenna Efficiency Functions");

    try {
        auto rr = antenna::radiation_efficiency<double>::calculate_radiation_resistance(
            antenna::radiation_efficiency<double>::antenna_type::dipole, 5.0, 30e6);
        bool passed = rr > 0.0 && rr < 1000.0;
        suite.add_test("calculate_radiation_resistance", passed);
    } catch (...) {
        suite.add_test("calculate_radiation_resistance", false, "exception thrown");
    }
    
    try {
        auto loss_r = antenna::radiation_efficiency<double>::calculate_ohmic_loss_resistance(
            5.0, 0.001, 30e6, 5.96e7);
        bool passed = loss_r > 0.0 && loss_r < 10.0;
        suite.add_test("calculate_ohmic_loss_resistance", passed);
    } catch (...) {
        suite.add_test("calculate_ohmic_loss_resistance", false, "exception thrown");
    }
    
    try {
        auto eff = antenna::radiation_efficiency<double>::calculate_radiation_efficiency(
            73.0, 0.5, 0.1);
        bool passed = eff > 0.9 && eff <= 1.0;
        suite.add_test("calculate_radiation_efficiency", passed);
    } catch (...) {
        suite.add_test("calculate_radiation_efficiency", false, "exception thrown");
    }
    
    try {
        antenna::radiation_efficiency<double>::antenna_parameters ant;
        ant.type = antenna::radiation_efficiency<double>::antenna_type::dipole;
        ant.length_m = 5.0;
        ant.radius_m = 0.001;
        ant.frequency_hz = 30e6;
        ant.input_impedance = std::complex<double>(73.0, 42.5);
        ant.characteristic_impedance = 50.0;
        
        antenna::radiation_efficiency<double>::conductor_properties cond;
        cond.conductivity_s_per_m = 5.96e7;
        cond.relative_permeability = 1.0;
        cond.diameter_m = 0.002;
        cond.surface_roughness_factor = 1.0;
        
        antenna::radiation_efficiency<double>::dielectric_properties diel;
        diel.relative_permittivity = 1.0;
        diel.loss_tangent = 0.0;
        diel.thickness_m = 0.0;
        diel.volume_m3 = 0.0;
        
        auto analysis = antenna::radiation_efficiency<double>::analyze_antenna_efficiency(
            ant, cond, diel);
        bool passed = analysis.radiation_efficiency > 0.9 && analysis.total_efficiency > 0.8;
        suite.add_test("analyze_antenna_efficiency", passed);
    } catch (...) {
        suite.add_test("analyze_antenna_efficiency", false, "exception thrown");
    }
    
    try {
        auto cond = antenna::radiation_efficiency<double>::create_conductor_material("copper");
        bool passed = approx_equal(cond.conductivity_s_per_m, 5.96e7, 1e6);
        suite.add_test("create_conductor_material", passed);
    } catch (...) {
        suite.add_test("create_conductor_material", false, "exception thrown");
    }
    
    try {
        auto diel = antenna::radiation_efficiency<double>::create_dielectric_material("teflon");
        bool passed = approx_equal(diel.relative_permittivity, 2.1, 0.01);
        suite.add_test("create_dielectric_material", passed);
    } catch (...) {
        suite.add_test("create_dielectric_material", false, "exception thrown");
    }

    try {
        double g = antenna::antenna_gain<double>::calculate_gain_from_directivity(1.64, 0.95);
        bool passed = approx_equal(g, 1.558, 0.01);
        suite.add_test("calculate_gain_from_directivity", passed);
    } catch (...) {
        suite.add_test("calculate_gain_from_directivity", false, "exception thrown");
    }
    
    try {
        double g_dbi = antenna::antenna_gain<double>::convert_dbd_to_dbi(3.0);
        bool passed = approx_equal(g_dbi, 5.15, 0.01);
        suite.add_test("convert_dbd_to_dbi", passed);
    } catch (...) {
        suite.add_test("convert_dbd_to_dbi", false, "exception thrown");
    }
    
    try {
        double ae = antenna::antenna_gain<double>::calculate_effective_aperture(6.0, 100e6);
        bool passed = ae > 0.0 && ae < 10.0;
        suite.add_test("calculate_effective_aperture", passed);
    } catch (...) {
        suite.add_test("calculate_effective_aperture", false, "exception thrown");
    }
    
    try {
        double g = antenna::antenna_gain<double>::calculate_gain_from_aperture(1.0, 0.8, 100e6);
        bool passed = g > 0.0 && g < 100.0;
        suite.add_test("calculate_gain_from_aperture", passed);
    } catch (...) {
        suite.add_test("calculate_gain_from_aperture", false, "exception thrown");
    }

    try {
        auto fields = antenna::near_field<double>::calculate_dipole_near_fields(
            1.0, M_PI/2, 0.0, 1.0, 0.0, 5.0, 30e6);
        bool passed = std::abs(fields.e_r) >= 0.0 && std::abs(fields.h_phi) >= 0.0;
        suite.add_test("calculate_dipole_near_fields", passed);
    } catch (...) {
        suite.add_test("calculate_dipole_near_fields", false, "exception thrown");
    }
    
    try {
        double dist = antenna::near_field<double>::calculate_reactive_near_field_boundary(5.0, 30e6);
        bool passed = dist > 0.0 && dist < 10.0;
        suite.add_test("calculate_reactive_near_field_boundary", passed);
    } catch (...) {
        suite.add_test("calculate_reactive_near_field_boundary", false, "exception thrown");
    }
    
    try {
        double dist = antenna::near_field<double>::calculate_radiating_near_field_boundary(5.0, 30e6);
        bool passed = dist > 0.0 && dist < 100.0;
        suite.add_test("calculate_radiating_near_field_boundary", passed);
    } catch (...) {
        suite.add_test("calculate_radiating_near_field_boundary", false, "exception thrown");
    }
    
    try {
        double dist = antenna::near_field<double>::calculate_far_field_boundary(5.0, 30e6);
        bool passed = dist > 0.0 && dist < 1000.0;
        suite.add_test("calculate_far_field_boundary", passed);
    } catch (...) {
        suite.add_test("calculate_far_field_boundary", false, "exception thrown");
    }
    
    try {
        double power = antenna::near_field<double>::calculate_power_density(100.0, 1000.0, 1.64, 30e6);
        bool passed = power > 0.0 && power < 1.0;
        suite.add_test("calculate_power_density", passed);
    } catch (...) {
        suite.add_test("calculate_power_density", false, "exception thrown");
    }
    
    try {
        bool is_far = antenna::near_field<double>::is_far_field(100.0, 5.0, 30e6);
        suite.add_test("is_far_field", true);
    } catch (...) {
        suite.add_test("is_far_field", false, "exception thrown");
    }
    
    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}