#include "src/propagation/ionospheric_ray_tracing.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

using namespace rvl::propagation;

void test_ionospheric_profile_creation() {
    std::cout << "Testing Ionospheric Profile Creation..." << std::endl;

    auto profile = ionospheric_ray_tracing_d::create_chapman_profile(10e6, 300.0, 50.0);

    std::cout << "  Profile heights: " << profile.height_km.size() << " points" << std::endl;
    std::cout << "  Height range: " << profile.height_km.front() << " - " 
              << profile.height_km.back() << " km" << std::endl;

    auto max_ne_it = std::max_element(profile.electron_density.begin(), 
                                     profile.electron_density.end());
    size_t peak_idx = std::distance(profile.electron_density.begin(), max_ne_it);

    std::cout << "  Peak electron density: " << *max_ne_it << " m⁻³" << std::endl;
    std::cout << "  Peak height: " << profile.height_km[peak_idx] << " km" << std::endl;

    std::cout << "  Magnetic field: " << profile.earth_magnetic_field_t * 1e6 << " μT" << std::endl;
    std::cout << "  Inclination: " << profile.magnetic_inclination_rad * 180.0 / M_PI << "°" << std::endl;

    assert(profile.height_km.size() == profile.electron_density.size());
    assert(profile.height_km.size() == profile.collision_frequency.size());
    assert(profile.height_km.front() < profile.height_km.back());
    assert(*max_ne_it > 0.0);

    std::cout << "  Ionospheric profile tests passed" << std::endl;
}

void test_ray_integration() {
    std::cout << "Testing Ray Integration..." << std::endl;

    auto profile = ionospheric_ray_tracing_d::create_chapman_profile(12e6, 320.0, 45.0);

    ionospheric_ray_tracing_d::ray_state initial_state;
    initial_state.position = {6371.0, 0.0, 0.0}; 
    initial_state.direction = {0.707, 0.0, 0.707}; 
    initial_state.group_path = 0.0;
    initial_state.phase_path = 0.0;
    initial_state.time_delay = 0.0;
    initial_state.is_valid = true;

    double test_frequency = 14.2e6;

    auto next_state = ionospheric_ray_tracing_d::integrate_ray_step(
        initial_state, profile, test_frequency, 1.0);

    std::cout << "  Initial position: [" << initial_state.position[0] 
              << ", " << initial_state.position[1] 
              << ", " << initial_state.position[2] << "] km" << std::endl;
    std::cout << "  Next position: [" << next_state.position[0] 
              << ", " << next_state.position[1] 
              << ", " << next_state.position[2] << "] km" << std::endl;

    double initial_height = std::sqrt(initial_state.position[0]*initial_state.position[0] + 
                                    initial_state.position[1]*initial_state.position[1] + 
                                    initial_state.position[2]*initial_state.position[2]) - 6371.0;
    double next_height = std::sqrt(next_state.position[0]*next_state.position[0] + 
                                 next_state.position[1]*next_state.position[1] + 
                                 next_state.position[2]*next_state.position[2]) - 6371.0;

    std::cout << "  Height change: " << initial_height << " → " << next_height << " km" << std::endl;
    std::cout << "  Group path: " << next_state.group_path << " km" << std::endl;
    std::cout << "  Phase path: " << next_state.phase_path << " km" << std::endl;

    assert(next_state.is_valid);
    assert(next_height > initial_height); 
    assert(next_state.group_path > 0.0);
    assert(next_state.phase_path > 0.0);

    std::cout << "  Ray integration tests passed" << std::endl;
}

void test_ray_tracing() {
    std::cout << "Testing Complete Ray Tracing..." << std::endl;

    auto profile = ionospheric_ray_tracing_d::create_chapman_profile(10e6, 300.0, 50.0);

    ionospheric_ray_tracing_d::vector3_type tx_pos = {6371.0, 0.0, 0.0}; 
    ionospheric_ray_tracing_d::vector3_type rx_pos = {6371.0 + 1500.0, 0.0, 0.0}; 

    double frequency = 7.1e6;

    std::vector<double> test_elevations = {5.0, 15.0, 30.0, 45.0, 60.0}; 

    for (double elev_deg : test_elevations) {
        double elevation_rad = elev_deg * M_PI / 180.0;

        auto ray_path = ionospheric_ray_tracing_d::trace_ray_path(
            tx_pos, rx_pos, elevation_rad, 0.0, frequency, profile, 5.0, 2000);

        if (!ray_path.empty()) {

            double max_height = 0.0;
            double final_distance = 0.0;

            for (const auto& state : ray_path) {
                double height = std::sqrt(state.position[0]*state.position[0] + 
                                        state.position[1]*state.position[1] + 
                                        state.position[2]*state.position[2]) - 6371.0;
                max_height = std::max(max_height, height);
            }

            if (!ray_path.empty()) {
                auto final_pos = ray_path.back().position;
                final_distance = std::sqrt(std::pow(final_pos[0] - rx_pos[0], 2) +
                                         std::pow(final_pos[1] - rx_pos[1], 2) +
                                         std::pow(final_pos[2] - rx_pos[2], 2));
            }

            std::cout << "  Elevation " << std::setw(2) << (int)elev_deg << "°: "
                      << std::setw(4) << ray_path.size() << " steps, "
                      << "max height " << std::setw(6) << std::setprecision(1) << std::fixed << max_height << " km, "
                      << "final distance " << std::setw(6) << std::setprecision(1) << final_distance << " km";

            if (ray_path.back().is_valid) {
                std::cout << " [VALID]";
            } else {
                std::cout << " [INVALID]";
            }
            std::cout << std::endl;
        } else {
            std::cout << "  Elevation " << std::setw(2) << (int)elev_deg << "°: No ray path generated" << std::endl;
        }
    }

    std::cout << "  Ray tracing tests completed" << std::endl;
}

void test_muf_calculation() {
    std::cout << "Testing MUF Calculation..." << std::endl;

    auto profile = ionospheric_ray_tracing_d::create_chapman_profile(12e6, 300.0, 50.0);

    ionospheric_ray_tracing_d::vector3_type tx_pos = {6371.0, 0.0, 0.0};
    ionospheric_ray_tracing_d::vector3_type rx_pos = {6371.0 + 1800.0, 0.0, 0.0}; 

    auto muf_analysis = ionospheric_ray_tracing_d::calculate_muf(tx_pos, rx_pos, profile, 90);

    std::cout << "  Path Analysis (2000 km):" << std::endl;
    std::cout << "    Critical frequency (foF2): " 
              << muf_analysis.critical_frequency_hz / 1e6 << " MHz" << std::endl;
    std::cout << "    Maximum Usable Frequency: " 
              << muf_analysis.maximum_usable_frequency_hz / 1e6 << " MHz" << std::endl;
    std::cout << "    Optimum Working Frequency: " 
              << muf_analysis.optimum_working_frequency_hz / 1e6 << " MHz" << std::endl;
    std::cout << "    Lowest Usable Frequency: " 
              << muf_analysis.lowest_usable_frequency_hz / 1e6 << " MHz" << std::endl;
    std::cout << "    Penetration Frequency: " 
              << muf_analysis.penetration_frequency_hz / 1e6 << " MHz" << std::endl;

    double muf_to_fof2_ratio = muf_analysis.maximum_usable_frequency_hz / 
                              muf_analysis.critical_frequency_hz;
    std::cout << "    MUF/foF2 ratio: " << muf_to_fof2_ratio << std::endl;

    assert(muf_analysis.critical_frequency_hz > 0.0);
    assert(muf_analysis.maximum_usable_frequency_hz > muf_analysis.critical_frequency_hz);
    assert(muf_analysis.optimum_working_frequency_hz < muf_analysis.maximum_usable_frequency_hz);
    assert(muf_to_fof2_ratio > 1.0 && muf_to_fof2_ratio < 10.0); 

    std::cout << "  MUF calculation tests passed" << std::endl;
}

void test_multi_hop_calculation() {
    std::cout << "Testing Multi-Hop Calculation..." << std::endl;

    auto profile = ionospheric_ray_tracing_d::create_chapman_profile(11e6, 320.0, 45.0);

    ionospheric_ray_tracing_d::vector3_type tx_pos = {6371.0, 0.0, 0.0};
    ionospheric_ray_tracing_d::vector3_type rx_pos = {6371.0 + 7500.0, 0.0, 0.0};

    double frequency = 14.2e6;

    auto hops = ionospheric_ray_tracing_d::calculate_multi_hop_path(
        tx_pos, rx_pos, frequency, profile, 5);

    std::cout << "  Long Distance Path Analysis (8000 km @ 14.2 MHz):" << std::endl;
    std::cout << "    Number of hops: " << hops.size() << std::endl;

    double total_distance = 0.0;
    double total_delay = 0.0;
    double total_loss = 0.0;

    for (size_t i = 0; i < hops.size(); ++i) {
        const auto& hop = hops[i];

        std::cout << "    Hop " << (i + 1) << ": ";

        if (hop.is_valid_hop) {
            std::cout << std::setprecision(0) << hop.hop_distance_km << " km, "
                      << "elevation " << std::setprecision(1) << hop.launch_angle_rad * 180.0 / M_PI << "°, "
                      << "max height " << std::setprecision(0) << hop.maximum_height_km << " km, "
                      << "delay " << std::setprecision(2) << hop.group_delay_ms << " ms";

            total_distance += hop.hop_distance_km;
            total_delay += hop.group_delay_ms;
            total_loss += hop.absorption_loss_db;
        } else {
            std::cout << "FAILED";
        }
        std::cout << std::endl;
    }

    if (!hops.empty() && hops.back().is_valid_hop) {
        std::cout << "    Total distance: " << std::setprecision(0) << total_distance << " km" << std::endl;
        std::cout << "    Total delay: " << std::setprecision(2) << total_delay << " ms" << std::endl;
        std::cout << "    Total path loss: " << std::setprecision(1) << total_loss << " dB" << std::endl;
    }

    if (!hops.empty()) {
        assert(hops[0].hop_distance_km > 0.0);

        if (total_distance > 5000.0) {
            assert(hops.size() >= 2);
        }
    }

    std::cout << "  Multi-hop calculation tests completed" << std::endl;
}

void test_ham_radio_scenarios() {
    std::cout << "Testing Ham Radio Scenarios..." << std::endl;

    std::vector<std::pair<std::string, ionospheric_ray_tracing_d::ionospheric_profile>> conditions = {
        {"Daytime High Solar", ionospheric_ray_tracing_d::create_chapman_profile(15e6, 350.0, 40.0)},
        {"Daytime Low Solar", ionospheric_ray_tracing_d::create_chapman_profile(8e6, 280.0, 60.0)},
        {"Nighttime", ionospheric_ray_tracing_d::create_chapman_profile(5e6, 250.0, 70.0)}
    };

    std::vector<std::pair<std::string, double>> ham_bands = {
        {"80m", 3.7e6},
        {"40m", 7.1e6},
        {"20m", 14.2e6},
        {"15m", 21.2e6},
        {"10m", 28.3e6}
    };

    ionospheric_ray_tracing_d::vector3_type tx_pos = {6371.0, 0.0, 0.0};
    ionospheric_ray_tracing_d::vector3_type rx_pos = {6371.0 + 2800.0, 0.0, 0.0};

    std::cout << "  Propagation Analysis (3000 km path):" << std::endl;
    std::cout << "  " << std::setw(12) << "Condition" 
              << std::setw(8) << "Band" 
              << std::setw(12) << "MUF (MHz)" 
              << std::setw(12) << "OWF (MHz)" 
              << std::setw(10) << "Usable?" << std::endl;

    for (const auto& condition : conditions) {
        for (const auto& band : ham_bands) {
            auto muf_analysis = ionospheric_ray_tracing_d::calculate_muf(
                tx_pos, rx_pos, condition.second, 45);

            bool band_usable = (band.second <= muf_analysis.maximum_usable_frequency_hz) &&
                              (band.second >= muf_analysis.lowest_usable_frequency_hz);

            std::cout << "  " << std::setw(12) << condition.first 
                      << std::setw(8) << band.first
                      << std::setw(12) << std::setprecision(1) << std::fixed 
                      << muf_analysis.maximum_usable_frequency_hz / 1e6
                      << std::setw(12) << std::setprecision(1) 
                      << muf_analysis.optimum_working_frequency_hz / 1e6
                      << std::setw(10) << (band_usable ? "YES" : "NO") << std::endl;
        }
    }

    std::cout << "  Ham radio scenario tests completed" << std::endl;
}

int main() {
    std::cout << "=== Ionospheric Ray Tracing Test Suite ===" << std::endl;
    std::cout << std::endl;

    try {
        test_ionospheric_profile_creation();
        std::cout << std::endl;

        test_ray_integration();
        std::cout << std::endl;

        test_ray_tracing();
        std::cout << std::endl;

        test_muf_calculation();
        std::cout << std::endl;

        test_multi_hop_calculation();
        std::cout << std::endl;

        test_ham_radio_scenarios();
        std::cout << std::endl;

        std::cout << "=== ALL IONOSPHERIC RAY TRACING TESTS PASSED! ===" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}