#include <iostream>
#include <iomanip>

#include "../src/antenna/dipole_resonant_frequency.hpp"
#include "../src/antenna/dipole_impedance.hpp"
#include "../src/antenna/antenna_gain.hpp"
#include "../src/antenna/dipole_radiation_pattern.hpp"
#include "../src/antenna/mutual_impedance.hpp"

#include "../src/rf_systems/scattering_parameters.hpp"

#include "../src/propagation/free_space/fspl.hpp"
#include "../src/propagation/free_space/friis.hpp"
#include "../src/propagation/ground_wave/ground_wave_loss.hpp"

#include "../src/feedline/transmission_line_impedance.hpp"
#include "../src/feedline/vswr.hpp"
#include "../src/feedline/transmission_line_loss.hpp"

using namespace rvl;

int main() {
    std::cout << "RadioVectorLib Basic Equations Test\n";
    std::cout << "===================================\n\n";

    std::cout << std::fixed << std::setprecision(3);

    try {

        std::cout << "Antenna Equations:\n";

        double freq = antenna::dipole_resonant_frequency_d::calculate_frequency(0.15); 
        std::cout << "  Half-wave dipole (15cm): " << freq / 1e6 << " MHz\n";

        auto impedance = antenna::dipole_impedance_d::ideal_halfwave_impedance();
        std::cout << "  Half-wave dipole impedance: " << impedance.real() << " + j" << impedance.imag() << " Ω\n";

        double gain_dbi = antenna::antenna_gain_d::numeric_to_dbi(1.64);
        std::cout << "  Half-wave dipole gain: " << gain_dbi << " dBi\n";

        double pattern_90deg = antenna::dipole_radiation_pattern_d::electric_field_pattern(M_PI/2);
        std::cout << "  Dipole pattern at 90°: " << pattern_90deg << "\n";

        std::cout << "\nPropagation Equations:\n";

        double fspl = propagation::free_space::fspl_d::calculate_db(1000, 100e6); 
        std::cout << "  FSPL (1km, 100MHz): " << fspl << " dB\n";

        double rx_power = propagation::free_space::friis_d::received_power_dbm(30, 3, 3, 100e6, 1000);
        std::cout << "  Friis received power: " << rx_power << " dBm\n";

        double gw_loss = propagation::ground_wave::ground_wave_loss_d::calculate_db(10000, 3e6);
        std::cout << "  Ground wave loss (10km, 3MHz): " << gw_loss << " dB\n";

        std::cout << "\nFeedline Equations:\n";

        double z0_coax = feedline::transmission_line_impedance_d::coaxial_cable_impedance(0.001, 0.0046, 2.3);
        std::cout << "  RG-58 coax impedance: " << z0_coax << " Ω\n";

        std::complex<double> z_load(75, 25);
        double vswr_val = feedline::vswr_d::from_impedances(z_load, 50);
        std::cout << "  VSWR (75+j25Ω load on 50Ω line): " << vswr_val << "\n";

        double tl_loss = feedline::transmission_line_loss_d::calculate_db(0.1, 30); 
        std::cout << "  Transmission line loss (30m): " << tl_loss << " dB\n";

        std::cout << "\nRF System Equations:\n";

        std::complex<double> s11 = rf_systems::scattering_parameters_d::s11_from_impedance(z_load, 50);
        double return_loss = rf_systems::scattering_parameters_d::return_loss_db(s11);
        std::cout << "  S11 magnitude: " << std::abs(s11) << "\n";
        std::cout << "  Return loss: " << return_loss << " dB\n";

        std::cout << "\nAll basic equations tested successfully!\n";

    } catch (const std::exception& e) {
        std::cout << " Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}