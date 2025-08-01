#include "test_utils.hpp"
#include "src/feedline/transmission_line_impedance.hpp"
#include "src/feedline/transmission_line_loss.hpp"
#include "src/feedline/waveguide.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("Transmission Line Functions");

    try {
        double z0 = feedline::transmission_line_impedance<double>::calculate_coaxial_impedance(
            0.002, 0.0066, 2.25);
        bool passed = z0 > 40.0 && z0 < 60.0;
        suite.add_test("calculate_coaxial_impedance", passed);
    } catch (...) {
        suite.add_test("calculate_coaxial_impedance", false, "exception thrown");
    }
    
    try {
        double z0 = feedline::transmission_line_impedance<double>::calculate_microstrip_impedance(
            0.002, 1.6e-3, 4.4);
        bool passed = z0 > 30.0 && z0 < 100.0;
        suite.add_test("calculate_microstrip_impedance", passed);
    } catch (...) {
        suite.add_test("calculate_microstrip_impedance", false, "exception thrown");
    }
    
    try {
        double z0 = feedline::transmission_line_impedance<double>::calculate_stripline_impedance(
            0.002, 1.6e-3, 4.4);
        bool passed = z0 > 30.0 && z0 < 100.0;
        suite.add_test("calculate_stripline_impedance", passed);
    } catch (...) {
        suite.add_test("calculate_stripline_impedance", false, "exception thrown");
    }
    
    try {
        double z0 = feedline::transmission_line_impedance<double>::calculate_cpw_impedance(
            0.001, 0.002, 4.4);
        bool passed = z0 > 30.0 && z0 < 100.0;
        suite.add_test("calculate_cpw_impedance", passed);
    } catch (...) {
        suite.add_test("calculate_cpw_impedance", false, "exception thrown");
    }
    
    try {
        double width = feedline::transmission_line_impedance<double>::calculate_microstrip_width_for_impedance(
            50.0, 1.6e-3, 4.4);
        bool passed = width > 0.001 && width < 0.01;
        suite.add_test("calculate_microstrip_width_for_impedance", passed);
    } catch (...) {
        suite.add_test("calculate_microstrip_width_for_impedance", false, "exception thrown");
    }

    try {
        double loss = feedline::transmission_line_loss<double>::calculate_coaxial_loss_db_per_m(
            1e9, 0.002, 0.0066, 5.96e7, 2.25, 0.0001);
        bool passed = loss > 0.0 && loss < 1.0;
        suite.add_test("calculate_coaxial_loss_db_per_m", passed);
    } catch (...) {
        suite.add_test("calculate_coaxial_loss_db_per_m", false, "exception thrown");
    }
    
    try {
        double loss = feedline::transmission_line_loss<double>::calculate_microstrip_loss_db_per_m(
            1e9, 0.002, 1.6e-3, 5.96e7, 4.4, 0.02, 0.0);
        bool passed = loss > 0.0 && loss < 1.0;
        suite.add_test("calculate_microstrip_loss_db_per_m", passed);
    } catch (...) {
        suite.add_test("calculate_microstrip_loss_db_per_m", false, "exception thrown");
    }
    
    try {
        double loss = feedline::transmission_line_loss<double>::calculate_skin_depth_loss(
            1e9, 50.0, 5.96e7, 1.0, 0.001);
        bool passed = loss > 0.0 && loss < 1.0;
        suite.add_test("calculate_skin_depth_loss", passed);
    } catch (...) {
        suite.add_test("calculate_skin_depth_loss", false, "exception thrown");
    }
    
    try {
        double loss = feedline::transmission_line_loss<double>::calculate_dielectric_loss(
            1e9, 4.4, 0.02, 50.0, 0.005);
        bool passed = loss > 0.0 && loss < 1.0;
        suite.add_test("calculate_dielectric_loss", passed);
    } catch (...) {
        suite.add_test("calculate_dielectric_loss", false, "exception thrown");
    }
    
    try {
        double loss = feedline::transmission_line_loss<double>::calculate_rg58_loss_db_per_m(100e6);
        bool passed = approx_equal(loss, 0.066, 0.001);
        suite.add_test("calculate_rg58_loss_db_per_m", passed);
    } catch (...) {
        suite.add_test("calculate_rg58_loss_db_per_m", false, "exception thrown");
    }

    try {
        double fc = feedline::waveguide<double>::calculate_cutoff_frequency_te10(0.0229, 0.0102);
        bool passed = fc > 6e9 && fc < 7e9;
        suite.add_test("calculate_cutoff_frequency_te10", passed);
    } catch (...) {
        suite.add_test("calculate_cutoff_frequency_te10", false, "exception thrown");
    }
    
    try {
        double z0 = feedline::waveguide<double>::calculate_wave_impedance_te10(10e9, 0.0229, 0.0102);
        bool passed = z0 > 300.0 && z0 < 600.0;
        suite.add_test("calculate_wave_impedance_te10", passed);
    } catch (...) {
        suite.add_test("calculate_wave_impedance_te10", false, "exception thrown");
    }
    
    try {
        double vg = feedline::waveguide<double>::calculate_group_velocity(10e9, 6.5e9);
        bool passed = vg > 0.0 && vg < 3e8;
        suite.add_test("calculate_group_velocity", passed);
    } catch (...) {
        suite.add_test("calculate_group_velocity", false, "exception thrown");
    }
    
    try {
        double loss = feedline::waveguide<double>::calculate_attenuation_te10(
            10e9, 0.0229, 0.0102, 5.96e7);
        bool passed = loss > 0.0 && loss < 1.0;
        suite.add_test("calculate_attenuation_te10", passed);
    } catch (...) {
        suite.add_test("calculate_attenuation_te10", false, "exception thrown");
    }
    
    try {
        double power = feedline::waveguide<double>::calculate_power_handling_te10(
            0.0229, 0.0102, 30e3, 1.0);
        bool passed = power > 1e6 && power < 1e9;
        suite.add_test("calculate_power_handling_te10", passed);
    } catch (...) {
        suite.add_test("calculate_power_handling_te10", false, "exception thrown");
    }

    try {
        core::memory::simd_vector<double> frequencies(3);
        core::memory::simd_vector<double> losses(3);
        frequencies[0] = 100e6; frequencies[1] = 1e9; frequencies[2] = 10e9;
        
        feedline::transmission_line_loss<double>::calculate_coaxial_loss_batch(
            frequencies, 0.002, 0.0066, 5.96e7, 2.25, 0.0001, losses);
        bool passed = losses[0] < losses[1] && losses[1] < losses[2];
        suite.add_test("calculate_coaxial_loss_batch", passed);
    } catch (...) {
        suite.add_test("calculate_coaxial_loss_batch", false, "exception thrown");
    }
    
    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}