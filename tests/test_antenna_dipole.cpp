#include "test_utils.hpp"
#include "src/antenna/dipole_impedance.hpp"
#include "src/antenna/dipole_radiation_pattern.hpp"
#include "src/antenna/dipole_resonant_frequency.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("Dipole Antenna Functions");

    try {
        auto z = antenna::dipole_impedance<double>::ideal_halfwave_impedance();
        bool passed = approx_equal(z.real(), 73.0, 1.0) && approx_equal(z.imag(), 42.5, 1.0);
        suite.add_test("ideal_halfwave_impedance", passed);
    } catch (...) {
        suite.add_test("ideal_halfwave_impedance", false, "exception thrown");
    }
    
    try {
        auto z = antenna::dipole_impedance<double>::calculate_impedance(5.0, 30e6, 0.001);
        bool passed = std::abs(z) > 0.0 && std::abs(z) < 1000.0;
        suite.add_test("calculate_impedance", passed);
    } catch (...) {
        suite.add_test("calculate_impedance", false, "exception thrown");
    }
    
    try {
        double z = antenna::dipole_impedance<double>::radiation_resistance(5.0, 30e6);
        bool passed = z > 0.0 && z < 1000.0;
        suite.add_test("radiation_resistance", passed);
    } catch (...) {
        suite.add_test("radiation_resistance", false, "exception thrown");
    }
    
    try {
        double x = antenna::dipole_impedance<double>::calculate_reactance(5.0, 30e6, 0.001);
        bool passed = std::abs(x) < 1000.0;
        suite.add_test("calculate_reactance", passed);
    } catch (...) {
        suite.add_test("calculate_reactance", false, "exception thrown");
    }

    try {
        double g = antenna::dipole_radiation_pattern<double>::calculate_gain(M_PI/2, 0, 5.0, 30e6);
        bool passed = g >= 0.0 && g <= 10.0;
        suite.add_test("calculate_gain", passed);
    } catch (...) {
        suite.add_test("calculate_gain", false, "exception thrown");
    }
    
    try {
        double d = antenna::dipole_radiation_pattern<double>::calculate_directivity(5.0, 30e6);
        bool passed = d > 1.0 && d < 10.0;
        suite.add_test("calculate_directivity", passed);
    } catch (...) {
        suite.add_test("calculate_directivity", false, "exception thrown");
    }
    
    try {
        double ef = antenna::dipole_radiation_pattern<double>::e_field_magnitude(1000.0, M_PI/2, 0, 1.0, 5.0, 30e6);
        bool passed = ef > 0.0 && ef < 100.0;
        suite.add_test("e_field_magnitude", passed);
    } catch (...) {
        suite.add_test("e_field_magnitude", false, "exception thrown");
    }
    
    try {
        double ef = antenna::dipole_radiation_pattern<double>::near_field_intensity(1.0, M_PI/2, 0, 1.0, 5.0, 30e6);
        bool passed = ef > 0.0;
        suite.add_test("near_field_intensity", passed);
    } catch (...) {
        suite.add_test("near_field_intensity", false, "exception thrown");
    }

    try {
        double f = antenna::dipole_resonant_frequency<double>::halfwave_resonant_frequency(5.0);
        bool passed = approx_equal(f, 30e6, 1e6);
        suite.add_test("halfwave_resonant_frequency", passed);
    } catch (...) {
        suite.add_test("halfwave_resonant_frequency", false, "exception thrown");
    }
    
    try {
        double l = antenna::dipole_resonant_frequency<double>::halfwave_physical_length(30e6);
        bool passed = approx_equal(l, 5.0, 0.1);
        suite.add_test("halfwave_physical_length", passed);
    } catch (...) {
        suite.add_test("halfwave_physical_length", false, "exception thrown");
    }
    
    try {
        double l = antenna::dipole_resonant_frequency<double>::adjusted_physical_length(30e6, 0.001, 1.0);
        bool passed = l > 4.0 && l < 6.0;
        suite.add_test("adjusted_physical_length", passed);
    } catch (...) {
        suite.add_test("adjusted_physical_length", false, "exception thrown");
    }
    
    try {
        bool resonant = antenna::dipole_resonant_frequency<double>::is_resonant(5.0, 30e6, 0.1);
        suite.add_test("is_resonant", resonant);
    } catch (...) {
        suite.add_test("is_resonant", false, "exception thrown");
    }
    
    try {
        double bw = antenna::dipole_resonant_frequency<double>::calculate_bandwidth(30e6, 5.0, 0.001);
        bool passed = bw > 0.0 && bw < 10e6;
        suite.add_test("calculate_bandwidth", passed);
    } catch (...) {
        suite.add_test("calculate_bandwidth", false, "exception thrown");
    }
    
    try {
        double q = antenna::dipole_resonant_frequency<double>::calculate_q_factor(30e6, 5.0, 0.001);
        bool passed = q > 1.0 && q < 100.0;
        suite.add_test("calculate_q_factor", passed);
    } catch (...) {
        suite.add_test("calculate_q_factor", false, "exception thrown");
    }
    
    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}