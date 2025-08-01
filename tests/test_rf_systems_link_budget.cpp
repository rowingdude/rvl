#include "test_utils.hpp"
#include "src/rf_systems/link_budget.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("Link Budget Functions");

    try {
        double nf = rf_systems::link_budget<double>::calculate_noise_floor_dbm(1e6, 3.0, 290.0);
        bool passed = nf > -120.0 && nf < -100.0;
        suite.add_test("calculate_noise_floor_dbm", passed);
    } catch (...) {
        suite.add_test("calculate_noise_floor_dbm", false, "exception thrown");
    }
    
    try {
        double nf = rf_systems::link_budget<double>::calculate_noise_floor_watts(1e6, 3.0, 290.0);
        bool passed = nf > 0.0 && nf < 1e-10;
        suite.add_test("calculate_noise_floor_watts", passed);
    } catch (...) {
        suite.add_test("calculate_noise_floor_watts", false, "exception thrown");
    }

    try {
        double snr = rf_systems::link_budget<double>::calculate_snr_db(-50.0, -110.0);
        bool passed = approx_equal(snr, 60.0, 0.1);
        suite.add_test("calculate_snr_db", passed);
    } catch (...) {
        suite.add_test("calculate_snr_db", false, "exception thrown");
    }
    
    try {
        double snr = rf_systems::link_budget<double>::calculate_snr_from_cnr(30.0, 10e6, 1.0);
        bool passed = snr > 20.0 && snr < 40.0;
        suite.add_test("calculate_snr_from_cnr", passed);
    } catch (...) {
        suite.add_test("calculate_snr_from_cnr", false, "exception thrown");
    }

    try {
        double margin = rf_systems::link_budget<double>::calculate_link_margin_db(-70.0, -100.0, 10.0);
        bool passed = approx_equal(margin, 20.0, 0.1);
        suite.add_test("calculate_link_margin_db", passed);
    } catch (...) {
        suite.add_test("calculate_link_margin_db", false, "exception thrown");
    }

    try {
        double sens = rf_systems::link_budget<double>::calculate_sensitivity_dbm(10.0, 1e6, 3.0, 290.0);
        bool passed = sens > -120.0 && sens < -80.0;
        suite.add_test("calculate_sensitivity_dbm", passed);
    } catch (...) {
        suite.add_test("calculate_sensitivity_dbm", false, "exception thrown");
    }

    try {
        std::vector<double> nf_db = {3.0, 6.0, 10.0};
        std::vector<double> gain_db = {20.0, 15.0, 10.0};
        double total_nf = rf_systems::link_budget<double>::calculate_cascade_noise_figure_db(nf_db, gain_db);
        bool passed = total_nf > 3.0 && total_nf < 5.0;
        suite.add_test("calculate_cascade_noise_figure_db", passed);
    } catch (...) {
        suite.add_test("calculate_cascade_noise_figure_db", false, "exception thrown");
    }

    try {
        double temp = rf_systems::link_budget<double>::calculate_total_system_noise_temp(50.0, 3.0);
        bool passed = temp > 300.0 && temp < 1000.0;
        suite.add_test("calculate_total_system_noise_temp", passed);
    } catch (...) {
        suite.add_test("calculate_total_system_noise_temp", false, "exception thrown");
    }

    try {
        double evm = rf_systems::link_budget<double>::calculate_evm_from_snr_db(30.0);
        bool passed = evm > 0.0 && evm < 10.0;
        suite.add_test("calculate_evm_from_snr_db", passed);
    } catch (...) {
        suite.add_test("calculate_evm_from_snr_db", false, "exception thrown");
    }

    try {
        double sfdr = rf_systems::link_budget<double>::calculate_sfdr_db(-110.0, -174.0, 10e6);
        bool passed = sfdr > 50.0 && sfdr < 100.0;
        suite.add_test("calculate_sfdr_db", passed);
    } catch (...) {
        suite.add_test("calculate_sfdr_db", false, "exception thrown");
    }

    try {
        rf_systems::link_budget<double>::link_parameters params;
        params.tx_power_dbm = 20.0;
        params.tx_antenna_gain_dbi = 3.0;
        params.tx_cable_loss_db = 2.0;
        params.path_loss_db = 100.0;
        params.rx_antenna_gain_dbi = 6.0;
        params.rx_cable_loss_db = 1.0;
        params.rx_noise_figure_db = 3.0;
        params.bandwidth_hz = 1e6;
        params.required_snr_db = 10.0;
        params.implementation_loss_db = 2.0;
        params.antenna_temperature_k = 290.0;
        
        auto result = rf_systems::link_budget<double>::calculate_complete_link_budget(params);
        bool passed = result.received_power_dbm < 0.0 && result.link_margin_db > -50.0;
        suite.add_test("calculate_complete_link_budget", passed);
    } catch (...) {
        suite.add_test("calculate_complete_link_budget", false, "exception thrown");
    }

    try {
        double ber = rf_systems::link_budget<double>::calculate_ber_bpsk(10.0);
        bool passed = ber > 0.0 && ber < 1.0;
        suite.add_test("calculate_ber_bpsk", passed);
    } catch (...) {
        suite.add_test("calculate_ber_bpsk", false, "exception thrown");
    }
    
    try {
        double ber = rf_systems::link_budget<double>::calculate_ber_qpsk(10.0);
        bool passed = ber > 0.0 && ber < 1.0;
        suite.add_test("calculate_ber_qpsk", passed);
    } catch (...) {
        suite.add_test("calculate_ber_qpsk", false, "exception thrown");
    }

    try {
        core::memory::simd_vector<double> powers(3);
        core::memory::simd_vector<double> snrs(3);
        powers[0] = -50.0; powers[1] = -60.0; powers[2] = -70.0;
        
        rf_systems::link_budget<double>::calculate_snr_batch(powers, -110.0, snrs);
        bool passed = snrs[0] > snrs[1] && snrs[1] > snrs[2];
        suite.add_test("calculate_snr_batch", passed);
    } catch (...) {
        suite.add_test("calculate_snr_batch", false, "exception thrown");
    }
    
    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}