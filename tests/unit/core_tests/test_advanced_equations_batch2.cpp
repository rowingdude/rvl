#include <iostream>
#include <iomanip>

#include "../src/propagation/hata_model.hpp"
#include "../src/propagation/tropospheric/ducting.hpp"

#include "../src/feedline/twin_lead.hpp"
#include "../src/feedline/waveguide.hpp"

#include "../src/rf_systems/intermodulation.hpp"
#include "../src/rf_systems/phase_noise.hpp"

using namespace rvl;

int main() {
    std::cout << "RadioVectorLib Advanced Equations Batch 2 Test\n";
    std::cout << "===============================================\n\n";

    std::cout << std::fixed << std::setprecision(3);

    try {
        std::cout << "Advanced Propagation Equations:\n";

        std::cout << "  Hata Model Tests:\n";
        double hata_urban = propagation::hata_model_d::path_loss_urban_db(900, 50, 1.5, 5); 
        std::cout << "    Urban path loss (900MHz, 5km): " << hata_urban << " dB\n";

        double hata_suburban = propagation::hata_model_d::path_loss_suburban_db(900, 50, 1.5, 5);
        std::cout << "    Suburban path loss (900MHz, 5km): " << hata_suburban << " dB\n";

        double hata_rural = propagation::hata_model_d::path_loss_rural_open_db(900, 50, 1.5, 5);
        std::cout << "    Rural path loss (900MHz, 5km): " << hata_rural << " dB\n";

        double coverage_radius = propagation::hata_model_d::coverage_radius_km(40, 15, 0, -100, 900, 50, 1.5,
                                                                              propagation::hata_model_d::environment_type::URBAN);
        std::cout << "    Coverage radius (40dBm EIRP, -100dBm sens): " << coverage_radius << " km\n";

        std::cout << "\n  Tropospheric Ducting Tests:\n";
        double duct_height = propagation::tropospheric::ducting_d::duct_height(0.1, -157e-9, 1.0003); 
        std::cout << "    Duct height (3GHz, strong gradient): " << duct_height << " m\n";

        double evap_duct = propagation::tropospheric::ducting_d::evaporation_duct_height_m(20, 15, 70, 5); 
        std::cout << "    Evaporation duct height: " << evap_duct << " m\n";

        double max_range = propagation::tropospheric::ducting_d::maximum_range_km(30, 3e9, 10); 
        std::cout << "    Maximum ducting range: " << max_range << " km\n";

        double ducting_prob = propagation::tropospheric::ducting_d::ducting_probability_percent(3, 100, 2.0); 
        std::cout << "    Ducting probability (3GHz, 100km, maritime): " << ducting_prob << "%\n";

        double signal_enhancement = propagation::tropospheric::ducting_d::signal_enhancement_db(3, 30, 5);
        std::cout << "    Signal enhancement: " << signal_enhancement << " dB\n";

        std::cout << "\nAdvanced Feedline Equations:\n";

        std::cout << "  Twin Lead Tests:\n";
        double twin_300_ohm = feedline::twin_lead_d::characteristic_impedance(0.007, 0.0016, 2.25); 
        std::cout << "    300Ω twin lead impedance: " << twin_300_ohm << " Ω\n";

        double twin_450_ohm = feedline::twin_lead_d::characteristic_impedance(0.013, 0.0016, 1.0); 
        std::cout << "    450Ω twin lead impedance: " << twin_450_ohm << " Ω\n";

        double velocity_factor = feedline::twin_lead_d::velocity_factor(2.25);
        std::cout << "    Velocity factor (polyethylene): " << velocity_factor << "\n";

        double prop_delay = feedline::twin_lead_d::propagation_delay_ns_per_m(2.25);
        std::cout << "    Propagation delay: " << prop_delay << " ns/m\n";

        double twin_loss = feedline::twin_lead_d::loss_db_per_m(100e6, 0.007, 0.0016, 5.8e7, 2.25, 0.001);
        std::cout << "    Loss at 100MHz: " << twin_loss * 1000 << " dB/km\n";

        std::cout << "\n  Waveguide Tests:\n";
        double wr90_fc = feedline::waveguide_d::cutoff_frequency_te10(0.02286); 
        std::cout << "    WR-90 cutoff frequency: " << wr90_fc / 1e9 << " GHz\n";

        double wr90_impedance = feedline::waveguide_d::characteristic_impedance_te10(10e9, 0.02286); 
        std::cout << "    WR-90 impedance at 10GHz: " << wr90_impedance << " Ω\n";

        double guide_wavelength = feedline::waveguide_d::guide_wavelength_te10(10e9, 0.02286);
        std::cout << "    Guide wavelength at 10GHz: " << guide_wavelength * 1000 << " mm\n";

        double surface_resistance = feedline::waveguide_d::surface_resistance(10e9, 5.8e7);
        double wg_attenuation = feedline::waveguide_d::attenuation_te10_db_per_m(10e9, 0.02286, 0.01016, surface_resistance);
        std::cout << "    WR-90 attenuation at 10GHz: " << wg_attenuation << " dB/m\n";

        double power_handling = feedline::waveguide_d::power_handling_kw(0.02286, 0.01016, 10e9, 30);
        std::cout << "    Power handling capability: " << power_handling << " kW\n";

        bool single_mode = feedline::waveguide_d::is_single_mode(10e9, 0.02286, 0.01016);
        std::cout << "    Single mode at 10GHz: " << (single_mode ? "Yes" : "No") << "\n";

        std::cout << "\nSystem Parameter Equations:\n";

        std::cout << "  Intermodulation (IP3) Tests:\n";
        double input_ip3 = rf_systems::intermodulation_d::input_ip3_from_imd(-10, 60); 
        std::cout << "    Input IP3 (-10dBm input, 60dB IMD): " << input_ip3 << " dBm\n";

        double output_ip3 = rf_systems::intermodulation_d::output_ip3_from_input_ip3(input_ip3, 20); 
        std::cout << "    Output IP3 (20dB gain): " << output_ip3 << " dBm\n";

        double imd_level = rf_systems::intermodulation_d::imd_level(-20, input_ip3);
        std::cout << "    IMD level at -20dBm input: " << imd_level << " dB\n";

        double sfdr = rf_systems::intermodulation_d::spurious_free_dynamic_range(input_ip3, -130, 20);
        std::cout << "    SFDR (IP3=" << input_ip3 << "dBm, NF=-130dBm): " << sfdr << " dB\n";

        double cascade_ip3 = rf_systems::intermodulation_d::cascade_ip3_two_stages(20, 30, 15); 
        std::cout << "    Cascaded IP3 (20dBm + 30dBm stages): " << cascade_ip3 << " dBm\n";

        std::cout << "\n  Phase Noise Tests:\n";
        double phase_noise_dbc = rf_systems::phase_noise_d::power_spectral_density_dbc_per_hz(1e-15, 1e-3); 
        std::cout << "    Phase noise (1fW/1mW): " << phase_noise_dbc << " dBc/Hz\n";

        double phase_jitter_rad = rf_systems::phase_noise_d::rms_phase_jitter_rad(-120, 1e6); 
        std::cout << "    RMS phase jitter: " << phase_jitter_rad * 1000 << " mrad\n";

        double timing_jitter_fs = rf_systems::phase_noise_d::timing_jitter_femtoseconds(-120, 1e9, 1e6); 
        std::cout << "    Timing jitter: " << timing_jitter_fs << " fs\n";

        double pn_from_q = rf_systems::phase_noise_d::phase_noise_from_q_factor(1e6, 100e6, 1000); 
        std::cout << "    Phase noise from Q (Q=1M, 100MHz, 1kHz offset): " << pn_from_q << " dBc/Hz\n";

        double additive_pn = rf_systems::phase_noise_d::additive_phase_noise_dbm_per_hz(-120, 10); 
        std::cout << "    Additive phase noise: " << additive_pn << " dBm/Hz\n";

        std::cout << "\nBatch Processing Tests:\n";

        core::memory::simd_vector<double> frequencies = {150, 450, 900, 1200};
        core::memory::simd_vector<double> tx_heights = {50, 50, 50, 50};
        core::memory::simd_vector<double> rx_heights = {1.5, 1.5, 1.5, 1.5};
        core::memory::simd_vector<double> distances = {2, 5, 10, 15};
        core::memory::simd_vector<double> path_losses(4);

        propagation::hata_model_d::path_loss_batch(frequencies, tx_heights, rx_heights, distances, path_losses);
        std::cout << "  Hata path losses (dB): ";
        for (size_t i = 0; i < path_losses.size(); ++i) {
            std::cout << path_losses[i] << " ";
        }
        std::cout << "\n";

        core::memory::simd_vector<double> spacings = {0.005, 0.007, 0.010, 0.013};
        core::memory::simd_vector<double> diameters = {0.0016, 0.0016, 0.0016, 0.0016};
        core::memory::simd_vector<double> permittivities = {2.25, 2.25, 1.0, 1.0};
        core::memory::simd_vector<double> impedances(4);

        feedline::twin_lead_d::characteristic_impedance_batch(spacings, diameters, permittivities, impedances);
        std::cout << "  Twin lead impedances (Ω): ";
        for (size_t i = 0; i < impedances.size(); ++i) {
            std::cout << impedances[i] << " ";
        }
        std::cout << "\n";

        core::memory::simd_vector<double> wg_dimensions = {0.02286, 0.01580, 0.01067, 0.00711}; 
        core::memory::simd_vector<double> cutoff_freqs(4);

        feedline::waveguide_d::cutoff_frequency_batch(wg_dimensions, cutoff_freqs);
        std::cout << "  Waveguide cutoff frequencies (GHz): ";
        for (size_t i = 0; i < cutoff_freqs.size(); ++i) {
            std::cout << cutoff_freqs[i] / 1e9 << " ";
        }
        std::cout << "\n";

        core::memory::simd_vector<double> offset_freqs = {100, 1000, 10000, 100000};
        core::memory::simd_vector<double> pn_profile(4);

        rf_systems::phase_noise_d::phase_noise_profile_batch(offset_freqs, pn_profile, 100e6, 1e6);
        std::cout << "  Phase noise profile (dBc/Hz): ";
        for (size_t i = 0; i < pn_profile.size(); ++i) {
            std::cout << pn_profile[i] << " ";
        }
        std::cout << "\n";

        std::cout << "\nAll advanced equations batch 2 tested successfully!\n";

    } catch (const std::exception& e) {
        std::cout << " Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}