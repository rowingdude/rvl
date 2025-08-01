#include "test_utils.hpp"
#include "src/propagation/free_space/fspl.hpp"
#include "src/propagation/free_space/friis.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("Free Space Propagation Functions");

    try {
        double loss = propagation::free_space::fspl<double>::calculate_loss_db(1000.0, 100e6);
        bool passed = approx_equal(loss, 72.45, 0.1);
        suite.add_test("fspl::calculate_loss_db", passed);
    } catch (...) {
        suite.add_test("fspl::calculate_loss_db", false, "exception thrown");
    }
    
    try {
        double loss = propagation::free_space::fspl<double>::calculate_loss_with_gains_db(1000.0, 100e6, 0.0, 0.0);
        bool passed = approx_equal(loss, 72.45, 0.1);
        suite.add_test("fspl::calculate_loss_with_gains_db", passed);
    } catch (...) {
        suite.add_test("fspl::calculate_loss_with_gains_db", false, "exception thrown");
    }
    
    try {
        double dist = propagation::free_space::fspl<double>::distance_from_loss(80.0, 100e6);
        bool passed = dist > 0.0 && dist < 10000.0;
        suite.add_test("fspl::distance_from_loss", passed);
    } catch (...) {
        suite.add_test("fspl::distance_from_loss", false, "exception thrown");
    }
    
    try {
        double fm = propagation::free_space::fspl<double>::fade_margin_db(80.0, 72.45);
        bool passed = approx_equal(fm, 7.55, 0.01);
        suite.add_test("fspl::fade_margin_db", passed);
    } catch (...) {
        suite.add_test("fspl::fade_margin_db", false, "exception thrown");
    }

    try {
        double pr = propagation::free_space::friis<double>::calculate_received_power_dbm(20.0, 0.0, 0.0, 1000.0, 100e6);
        bool passed = pr > -80.0 && pr < -40.0;
        suite.add_test("friis::calculate_received_power_dbm", passed);
    } catch (...) {
        suite.add_test("friis::calculate_received_power_dbm", false, "exception thrown");
    }
    
    try {
        double pr = propagation::free_space::friis<double>::calculate_received_power_watts(1.0, 1.0, 1.0, 1000.0, 100e6);
        bool passed = pr > 0.0 && pr < 1.0;
        suite.add_test("friis::calculate_received_power_watts", passed);
    } catch (...) {
        suite.add_test("friis::calculate_received_power_watts", false, "exception thrown");
    }
    
    try {
        double dist = propagation::free_space::friis<double>::calculate_max_range(20.0, -80.0, 0.0, 0.0, 100e6);
        bool passed = dist > 0.0 && dist < 100000.0;
        suite.add_test("friis::calculate_max_range", passed);
    } catch (...) {
        suite.add_test("friis::calculate_max_range", false, "exception thrown");
    }
    
    try {
        double pt = propagation::free_space::friis<double>::calculate_required_transmit_power_dbm(-80.0, 0.0, 0.0, 1000.0, 100e6);
        bool passed = pt > 0.0 && pt < 50.0;
        suite.add_test("friis::calculate_required_transmit_power_dbm", passed);
    } catch (...) {
        suite.add_test("friis::calculate_required_transmit_power_dbm", false, "exception thrown");
    }
    
    try {
        double eirp = propagation::free_space::friis<double>::calculate_eirp_dbm(20.0, 3.0);
        bool passed = approx_equal(eirp, 23.0, 0.01);
        suite.add_test("friis::calculate_eirp_dbm", passed);
    } catch (...) {
        suite.add_test("friis::calculate_eirp_dbm", false, "exception thrown");
    }
    
    try {
        double ae = propagation::free_space::friis<double>::calculate_effective_aperture(3.0, 100e6);
        bool passed = ae > 0.0 && ae < 10.0;
        suite.add_test("friis::calculate_effective_aperture", passed);
    } catch (...) {
        suite.add_test("friis::calculate_effective_aperture", false, "exception thrown");
    }

    try {
        core::memory::simd_vector<double> distances(3);
        core::memory::simd_vector<double> losses(3);
        distances[0] = 100.0; distances[1] = 1000.0; distances[2] = 10000.0;
        
        propagation::free_space::fspl<double>::calculate_loss_batch(distances, 100e6, losses);
        bool passed = losses[0] < losses[1] && losses[1] < losses[2];
        suite.add_test("fspl::calculate_loss_batch", passed);
    } catch (...) {
        suite.add_test("fspl::calculate_loss_batch", false, "exception thrown");
    }
    
    try {
        core::memory::simd_vector<double> distances(2);
        core::memory::simd_vector<double> powers(2);
        distances[0] = 100.0; distances[1] = 1000.0;
        
        propagation::free_space::friis<double>::calculate_received_power_batch(
            distances, 20.0, 0.0, 0.0, 100e6, powers);
        bool passed = powers[0] > powers[1];
        suite.add_test("friis::calculate_received_power_batch", passed);
    } catch (...) {
        suite.add_test("friis::calculate_received_power_batch", false, "exception thrown");
    }
    
    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}