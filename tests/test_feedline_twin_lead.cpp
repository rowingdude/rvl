#include "test_utils.hpp"
#include "src/feedline/twin_lead.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("Twin Lead Functions");

    try {
        double z0 = feedline::twin_lead<double>::characteristic_impedance(0.00762, 0.00163, 2.25);
        bool passed = approx_equal(z0, 300.0, 20.0); 
        suite.add_test("characteristic_impedance", passed);
    } catch (...) {
        suite.add_test("characteristic_impedance", false, "exception thrown");
    }

    try {
        double z0 = feedline::twin_lead<double>::characteristic_impedance_precise(0.00762, 0.00163, 2.25);
        bool passed = z0 > 200.0 && z0 < 400.0;
        suite.add_test("characteristic_impedance_precise", passed);
    } catch (...) {
        suite.add_test("characteristic_impedance_precise", false, "exception thrown");
    }

    try {
        double vf = feedline::twin_lead<double>::velocity_factor(2.25);
        bool passed = approx_equal(vf, 0.6667, 0.01);
        suite.add_test("velocity_factor", passed);
    } catch (...) {
        suite.add_test("velocity_factor", false, "exception thrown");
    }

    try {
        double delay = feedline::twin_lead<double>::propagation_delay_ns_per_m(2.25);
        bool passed = delay > 4.0 && delay < 6.0;
        suite.add_test("propagation_delay_ns_per_m", passed);
    } catch (...) {
        suite.add_test("propagation_delay_ns_per_m", false, "exception thrown");
    }

    try {
        double c = feedline::twin_lead<double>::capacitance_pf_per_m(0.00762, 0.00163, 2.25);
        bool passed = c > 10.0 && c < 50.0;
        suite.add_test("capacitance_pf_per_m", passed);
    } catch (...) {
        suite.add_test("capacitance_pf_per_m", false, "exception thrown");
    }

    try {
        double l = feedline::twin_lead<double>::inductance_nh_per_m(0.00762, 0.00163);
        bool passed = l > 100.0 && l < 1000.0;
        suite.add_test("inductance_nh_per_m", passed);
    } catch (...) {
        suite.add_test("inductance_nh_per_m", false, "exception thrown");
    }

    try {
        double loss = feedline::twin_lead<double>::loss_db_per_m(100e6, 0.00762, 0.00163, 5.8e7, 2.25, 0.0001);
        bool passed = loss > 0.0 && loss < 1.0;
        suite.add_test("loss_db_per_m", passed);
    } catch (...) {
        suite.add_test("loss_db_per_m", false, "exception thrown");
    }

    try {
        double loss = feedline::twin_lead<double>::radiation_loss_db_per_m(100e6, 0.00762);
        bool passed = loss >= 0.0 && loss < 0.1;
        suite.add_test("radiation_loss_db_per_m", passed);
    } catch (...) {
        suite.add_test("radiation_loss_db_per_m", false, "exception thrown");
    }

    try {
        core::memory::simd_vector<double> spacings(3);
        core::memory::simd_vector<double> diameters(3);
        core::memory::simd_vector<double> permittivities(3);
        core::memory::simd_vector<double> impedances(3);
        
        spacings[0] = 0.00762; diameters[0] = 0.00163; permittivities[0] = 2.25;
        spacings[1] = 0.0127; diameters[1] = 0.00163; permittivities[1] = 2.25;
        spacings[2] = 0.0254; diameters[2] = 0.00163; permittivities[2] = 1.0;
        
        feedline::twin_lead<double>::characteristic_impedance_batch(spacings, diameters, permittivities, impedances);
        bool passed = impedances[0] > 200.0 && impedances[1] > 300.0 && impedances[2] > 500.0;
        suite.add_test("characteristic_impedance_batch", passed);
    } catch (...) {
        suite.add_test("characteristic_impedance_batch", false, "exception thrown");
    }

    try {
        core::memory::simd_vector<double> frequencies(2);
        core::memory::simd_vector<double> spacings(2);
        core::memory::simd_vector<double> diameters(2);
        core::memory::simd_vector<double> losses(2);
        
        frequencies[0] = 100e6; spacings[0] = 0.00762; diameters[0] = 0.00163;
        frequencies[1] = 1e9; spacings[1] = 0.00762; diameters[1] = 0.00163;
        
        feedline::twin_lead<double>::loss_batch(frequencies, spacings, diameters, losses);
        bool passed = losses[0] > 0.0 && losses[1] > losses[0];
        suite.add_test("loss_batch", passed);
    } catch (...) {
        suite.add_test("loss_batch", false, "exception thrown");
    }

    try {
        double spacing = feedline::twin_lead<double>::standard_300_ohm_spacing_mm();
        bool passed = approx_equal(spacing, 7.62, 0.01);
        suite.add_test("standard_300_ohm_spacing_mm", passed);
    } catch (...) {
        suite.add_test("standard_300_ohm_spacing_mm", false, "exception thrown");
    }
    
    try {
        double diam = feedline::twin_lead<double>::standard_300_ohm_diameter_mm();
        bool passed = approx_equal(diam, 1.63, 0.01);
        suite.add_test("standard_300_ohm_diameter_mm", passed);
    } catch (...) {
        suite.add_test("standard_300_ohm_diameter_mm", false, "exception thrown");
    }

    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}