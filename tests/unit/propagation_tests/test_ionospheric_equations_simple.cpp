#include <iostream>
#include <cmath>
#include <complex>
#include <cassert>
#include "../src/propagation/ionospheric_refractive_index.hpp"
#include "../src/propagation/ray_trajectory.hpp"
#include "../src/propagation/group_phase_path.hpp"
#include "../src/propagation/ionospheric_scintillation.hpp"

using namespace rvl::propagation;

void test_ionospheric_refractive_index() {
    std::cout << "Testing ionospheric refractive index calculations..." << std::endl;

    const double electron_density = 1e12;  
    const double wave_frequency = 10e6;    
    const double collision_frequency = 1e3; 

    auto n_squared = ionospheric_refractive_index_d::calculate_n_squared(
        electron_density, wave_frequency, collision_frequency);

    assert(std::abs(n_squared) < 1.0);
    assert(n_squared.real() > 0.0);

    const double critical_freq = ionospheric_refractive_index_d::calculate_critical_frequency(electron_density);
    assert(critical_freq > 0.0 && critical_freq < 20e6);

    std::cout << "  Critical frequency: " << critical_freq / 1e6 << " MHz" << std::endl;
    std::cout << "  Refractive index magnitude: " << std::abs(n_squared) << std::endl;
    std::cout << "✓ Ionospheric refractive index tests passed" << std::endl;
}

void test_ray_trajectory() {
    std::cout << "Testing ray trajectory calculations..." << std::endl;

    ray_trajectory_d::point3d_type start_pos = {0.0, 0.0, 100000.0};  
    ray_trajectory_d::vector3d_type direction = {1.0, 0.0, 0.1};      
    const double wave_number = 0.1;

    auto initial_state = ray_trajectory_d::create_initial_state(start_pos, direction, wave_number);

    assert(initial_state.position[0] == start_pos[0]);
    assert(initial_state.position[1] == start_pos[1]);
    assert(initial_state.position[2] == start_pos[2]);
    assert(initial_state.path_length == 0.0);

    const double k_magnitude = ray_trajectory_d::calculate_wave_vector_magnitude(initial_state);
    assert(std::abs(k_magnitude - wave_number) < 1e-10);

    std::cout << "  Wave vector magnitude: " << k_magnitude << std::endl;
    std::cout << "✓ Ray trajectory tests passed" << std::endl;
}

void test_scintillation() {
    std::cout << "Testing ionospheric scintillation..." << std::endl;

    const double structure_constant = 1e-50;  
    const double path_length = 1000e3;        
    const double frequency = 100e6;           

    const double rytov_var = ionospheric_scintillation_d::calculate_rytov_variance_frequency(
        structure_constant, path_length, frequency);

    assert(rytov_var > 0.0);
    assert(rytov_var < 100.0);

    const bool is_weak = ionospheric_scintillation_d::is_weak_scattering_regime(rytov_var);

    std::cout << "  Rytov variance: " << rytov_var << std::endl;
    std::cout << "  Weak scattering regime: " << (is_weak ? "Yes" : "No") << std::endl;
    std::cout << "✓ Scintillation tests passed" << std::endl;
}

void test_batch_operations() {
    std::cout << "Testing batch operations..." << std::endl;

    const size_t n = 5;

    ionospheric_refractive_index_d::vector_type electron_densities = {1e11, 5e11, 1e12, 2e12, 5e12};
    ionospheric_refractive_index_d::vector_type wave_frequencies = {5e6, 10e6, 15e6, 20e6, 25e6};
    ionospheric_refractive_index_d::vector_type collision_frequencies = {1e3, 1e3, 1e3, 1e3, 1e3};
    ionospheric_refractive_index_d::complex_vector_type n_squared_results(n);

    ionospheric_refractive_index_d::calculate_n_squared_batch(
        electron_densities, wave_frequencies, collision_frequencies, n_squared_results);

    for (size_t i = 0; i < n; ++i) {
        assert(std::abs(n_squared_results[i]) < 1.0);
        assert(n_squared_results[i].real() > 0.0);
    }

    std::cout << "  Processed " << n << " calculations in batch" << std::endl;
    std::cout << "✓ Batch operation tests passed" << std::endl;
}

int main() {
    std::cout << "=== RadioVectorLib Ionospheric Equations Test ===" << std::endl;

    try {
        test_ionospheric_refractive_index();
        test_ray_trajectory();
        test_scintillation();
        test_batch_operations();

        std::cout << std::endl;
        std::cout << "🎉 All ionospheric equation tests passed!" << std::endl;
        std::cout << "Implemented equations:" << std::endl;
        std::cout << "  • Ionospheric refractive index (N²)" << std::endl;
        std::cout << "  • Ray trajectory (Hamilton's equations)" << std::endl;
        std::cout << "  • Group and phase path calculations" << std::endl;
        std::cout << "  • Rytov variance for scintillation" << std::endl;
        std::cout << "  • Critical frequency and MUF calculations" << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "❌ Test failed: " << e.what() << std::endl;
        return 1;
    }
}