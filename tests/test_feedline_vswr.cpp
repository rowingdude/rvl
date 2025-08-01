#include "test_utils.hpp"
#include "src/feedline/vswr.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("VSWR Functions");

    try {
        double vswr = feedline::vswr<double>::from_reflection_coefficient(0.2);
        bool passed = approx_equal(vswr, 1.5, 0.01);
        suite.add_test("vswr_from_reflection_coefficient", passed);
    } catch (...) {
        suite.add_test("vswr_from_reflection_coefficient", false, "exception thrown");
    }

    try {
        auto gamma = feedline::vswr<double>::reflection_coefficient_from_vswr(2.0);
        bool passed = approx_equal(std::abs(gamma), 0.3333, 0.001);
        suite.add_test("reflection_coefficient_from_vswr", passed);
    } catch (...) {
        suite.add_test("reflection_coefficient_from_vswr", false, "exception thrown");
    }

    try {
        std::complex<double> z_load(75.0, 0.0);
        std::complex<double> z0(50.0, 0.0);
        std::complex<double> gamma = (z_load - z0) / (z_load + z0);
        bool passed = approx_equal(std::abs(gamma), 0.2, 0.01);
        suite.add_test("complex_reflection_coefficient", passed);
    } catch (...) {
        suite.add_test("complex_reflection_coefficient", false, "exception thrown");
    }

    try {
        double vswr = feedline::vswr<double>::from_impedances(
            std::complex<double>(75.0, 0.0), 50.0);
        bool passed = approx_equal(vswr, 1.5, 0.01);
        suite.add_test("vswr_from_impedances", passed);
    } catch (...) {
        suite.add_test("vswr_from_impedances", false, "exception thrown");
    }

    try {
        double vswr_val = feedline::vswr<double>::from_reflection_coefficient(0.2);
        double rl = feedline::vswr<double>::return_loss_db(vswr_val);
        bool passed = approx_equal(rl, 13.98, 0.1);
        suite.add_test("return_loss_db", passed);
    } catch (...) {
        suite.add_test("return_loss_db", false, "exception thrown");
    }

    try {
        double rl = feedline::vswr<double>::return_loss_db(2.0);
        bool passed = rl > 9.0 && rl < 10.0;
        suite.add_test("return_loss_from_vswr_db", passed);
    } catch (...) {
        suite.add_test("return_loss_from_vswr_db", false, "exception thrown");
    }

    try {
        double ml = feedline::vswr<double>::mismatch_loss_db(2.0);
        bool passed = approx_equal(ml, 0.512, 0.01);
        suite.add_test("mismatch_loss_db", passed);
    } catch (...) {
        suite.add_test("mismatch_loss_db", false, "exception thrown");
    }

    try {
        double vswr_val = feedline::vswr<double>::from_reflection_coefficient(0.2);
        double tpf = feedline::vswr<double>::power_delivered_fraction(vswr_val);
        bool passed = approx_equal(tpf, 0.96, 0.01);
        suite.add_test("transmitted_power_fraction", passed);
    } catch (...) {
        suite.add_test("transmitted_power_fraction", false, "exception thrown");
    }

    try {
        double distance = feedline::vswr<double>::voltage_nodes_distance(2.0, 1.0);
        bool passed = distance > 0.0 && distance <= 0.5;
        suite.add_test("voltage_nodes_distance", passed);
    } catch (...) {
        suite.add_test("voltage_nodes_distance", false, "exception thrown");
    }

    try {
        core::memory::simd_vector<std::complex<double>> reflections(3);
        core::memory::simd_vector<double> vswrs(3);
        reflections[0] = std::complex<double>(0.0, 0.0); 
        reflections[1] = std::complex<double>(0.2, 0.0); 
        reflections[2] = std::complex<double>(0.5, 0.0);
        
        feedline::vswr<double>::from_reflection_coefficient_batch(reflections, vswrs);
        bool passed = approx_equal(vswrs[0], 1.0, 0.01) && 
                     approx_equal(vswrs[1], 1.5, 0.01) && 
                     approx_equal(vswrs[2], 3.0, 0.01);
        suite.add_test("vswr_batch", passed);
    } catch (...) {
        suite.add_test("vswr_batch", false, "exception thrown");
    }

    try {
        core::memory::simd_vector<double> vswrs(2);
        core::memory::simd_vector<double> return_losses(2);
        vswrs[0] = 1.5; vswrs[1] = 2.0;
        
        feedline::vswr<double>::return_loss_batch(vswrs, return_losses);
        bool passed = return_losses[0] > 13.0 && return_losses[1] > 9.0;
        suite.add_test("return_loss_batch", passed);
    } catch (...) {
        suite.add_test("return_loss_batch", false, "exception thrown");
    }
    
    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}