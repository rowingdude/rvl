#include <iostream>
#include <iomanip>

#include "../src/antenna/helical_antenna.hpp"
#include "../src/antenna/microstrip.hpp"
#include "../src/antenna/near_field.hpp"

#include "../src/propagation/free_space/friis_polarization.hpp"

using namespace rvl;

int main() {
    std::cout << "RadioVectorLib Advanced Equations Batch 1 Test\n";
    std::cout << "===============================================\n\n";

    std::cout << std::fixed << std::setprecision(3);

    try {
        std::cout << "Advanced Antenna Equations:\n";

        std::cout << "  Helical Antenna Tests:\n";
        double helix_gain = antenna::helical_antenna_d::gain_axial_mode(10, 0.075, 0.3, 0.3); 
        double helix_gain_db = antenna::helical_antenna_d::gain_axial_mode_db(10, 0.075, 0.3, 0.3);
        std::cout << "    Axial mode gain (10 turns): " << helix_gain << " (" << helix_gain_db << " dB)\n";

        double opt_circumference = antenna::helical_antenna_d::optimal_circumference(0.3);
        double opt_spacing = antenna::helical_antenna_d::optimal_turn_spacing(0.3);
        std::cout << "    Optimal circumference (λ=0.3m): " << opt_circumference << " m\n";
        std::cout << "    Optimal turn spacing (λ=0.3m): " << opt_spacing << " m\n";

        double pitch_angle = antenna::helical_antenna_d::pitch_angle_deg(0.075, 0.3);
        std::cout << "    Pitch angle: " << pitch_angle << "°\n";

        double beamwidth = antenna::helical_antenna_d::beamwidth_deg(10, 0.075, 0.3);
        std::cout << "    Beamwidth: " << beamwidth << "°\n";

        std::cout << "\n  Microstrip Patch Antenna Tests:\n";
        double eff_permittivity = antenna::microstrip_d::effective_dielectric_constant(4.4, 0.01, 0.0016); 
        std::cout << "    Effective permittivity (FR4): " << eff_permittivity << "\n";

        double patch_freq = antenna::microstrip_d::resonant_frequency(0.02, eff_permittivity); 
        std::cout << "    Resonant frequency (20mm patch): " << patch_freq / 1e9 << " GHz\n";

        double patch_length = antenna::microstrip_d::patch_length_for_frequency(2.4e9, eff_permittivity);
        std::cout << "    Patch length for 2.4GHz: " << patch_length * 1000 << " mm\n";

        double patch_width = antenna::microstrip_d::patch_width(2.4e9, 4.4, 0.0016);
        std::cout << "    Patch width for 2.4GHz: " << patch_width * 1000 << " mm\n";

        double inset_pos = antenna::microstrip_d::inset_feed_position(patch_length, 50.0);
        std::cout << "    Inset feed position (50Ω): " << inset_pos * 1000 << " mm\n";

        std::cout << "\n  Near Field Distance Tests:\n";
        double antenna_size = 1.0; 
        double wavelength = 0.1;   

        double fraunhofer_dist = antenna::near_field_d::fraunhofer_distance(antenna_size, wavelength);
        std::cout << "    Fraunhofer distance (1m antenna, 3GHz): " << fraunhofer_dist << " m\n";

        double fresnel_dist = antenna::near_field_d::fresnel_distance(antenna_size, wavelength);
        std::cout << "    Fresnel distance: " << fresnel_dist << " m\n";

        double reactive_boundary = antenna::near_field_d::reactive_near_field_boundary(antenna_size);
        std::cout << "    Reactive near field boundary: " << reactive_boundary << " m\n";

        auto field_region_5m = antenna::near_field_d::classify_field_region(5.0, antenna_size, wavelength);
        std::cout << "    Field region at 5m: " << (field_region_5m == antenna::near_field_d::field_region::FAR_FIELD ? "Far field" :
                                                     field_region_5m == antenna::near_field_d::field_region::RADIATING_NEAR_FIELD ? "Radiating near field" :
                                                     "Reactive near field") << "\n";

        double min_measurement_dist = antenna::near_field_d::minimum_measurement_distance(antenna_size, wavelength, 1.0);
        std::cout << "    Minimum measurement distance (1dB accuracy): " << min_measurement_dist << " m\n";

        std::cout << "\nAdvanced Propagation Equations:\n";

        std::cout << "  Friis with Polarization Mismatch Tests:\n";
        double rx_power_perfect = propagation::free_space::friis_polarization_d::received_power_dbm(
            30, 3, 3, 100e6, 1000, 0); 
        std::cout << "    Received power (perfect match): " << rx_power_perfect << " dBm\n";

        double rx_power_45deg = propagation::free_space::friis_polarization_d::received_power_dbm(
            30, 3, 3, 100e6, 1000, M_PI/4); 
        std::cout << "    Received power (45° mismatch): " << rx_power_45deg << " dBm\n";

        double rx_power_90deg = propagation::free_space::friis_polarization_d::received_power_dbm(
            30, 3, 3, 100e6, 1000, M_PI/2); 
        std::cout << "    Received power (90° mismatch): " << rx_power_90deg << " dBm\n";

        double pol_loss_45deg = propagation::free_space::friis_polarization_d::polarization_mismatch_loss_db(M_PI/4);
        std::cout << "    Polarization loss (45°): " << pol_loss_45deg << " dB\n";

        double pol_loss_90deg = propagation::free_space::friis_polarization_d::polarization_mismatch_loss_db(M_PI/2);
        std::cout << "    Polarization loss (90°): " << pol_loss_90deg << " dB\n";

        using pol_type = propagation::free_space::friis_polarization_d::polarization_type;
        double linear_to_circular_angle = propagation::free_space::friis_polarization_d::polarization_mismatch_angle(
            pol_type::LINEAR_VERTICAL, pol_type::CIRCULAR_RIGHT);
        double linear_to_circular_loss = propagation::free_space::friis_polarization_d::polarization_mismatch_loss_db(linear_to_circular_angle);
        std::cout << "    Linear to circular polarization loss: " << linear_to_circular_loss << " dB\n";

        std::cout << "\nBatch Processing Tests:\n";

        core::memory::simd_vector<double> turns = {8, 10, 12, 15};
        core::memory::simd_vector<double> spacings = {0.06, 0.075, 0.09, 0.11};
        core::memory::simd_vector<double> circumferences = {0.24, 0.3, 0.36, 0.45};
        core::memory::simd_vector<double> wavelengths = {0.3, 0.3, 0.3, 0.3};
        core::memory::simd_vector<double> helix_gains(4);

        antenna::helical_antenna_d::gain_axial_db_batch(turns, spacings, circumferences, wavelengths, helix_gains);
        std::cout << "  Helical gains (dB): ";
        for (size_t i = 0; i < helix_gains.size(); ++i) {
            std::cout << helix_gains[i] << " ";
        }
        std::cout << "\n";

        core::memory::simd_vector<double> patch_lengths = {0.015, 0.02, 0.025, 0.03};
        core::memory::simd_vector<double> rel_perms = {4.4, 4.4, 4.4, 4.4};
        core::memory::simd_vector<double> widths = {0.01, 0.01, 0.01, 0.01};
        core::memory::simd_vector<double> heights = {0.0016, 0.0016, 0.0016, 0.0016};
        core::memory::simd_vector<double> patch_freqs(4);

        antenna::microstrip_d::resonant_frequency_batch(patch_lengths, rel_perms, widths, heights, patch_freqs);
        std::cout << "  Microstrip frequencies (GHz): ";
        for (size_t i = 0; i < patch_freqs.size(); ++i) {
            std::cout << patch_freqs[i] / 1e9 << " ";
        }
        std::cout << "\n";

        core::memory::simd_vector<double> antenna_sizes = {0.5, 1.0, 1.5, 2.0};
        core::memory::simd_vector<double> test_wavelengths = {0.1, 0.1, 0.1, 0.1};
        core::memory::simd_vector<double> fraunhofer_distances(4);

        antenna::near_field_d::fraunhofer_distance_batch(antenna_sizes, test_wavelengths, fraunhofer_distances);
        std::cout << "  Fraunhofer distances (m): ";
        for (size_t i = 0; i < fraunhofer_distances.size(); ++i) {
            std::cout << fraunhofer_distances[i] << " ";
        }
        std::cout << "\n";

        core::memory::simd_vector<double> mismatch_angles = {0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};
        core::memory::simd_vector<double> pol_losses(5);

        propagation::free_space::friis_polarization_d::polarization_loss_batch(mismatch_angles, pol_losses);
        std::cout << "  Polarization losses (dB): ";
        for (size_t i = 0; i < pol_losses.size(); ++i) {
            std::cout << pol_losses[i] << " ";
        }
        std::cout << "\n";

        std::cout << "\nAll advanced equations batch 1 tested successfully!\n";

    } catch (const std::exception& e) {
        std::cout << " Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}