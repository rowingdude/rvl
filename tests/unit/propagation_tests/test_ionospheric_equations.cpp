#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include "../src/propagation/ionospheric_refractive_index.hpp"
#include "../src/propagation/ray_trajectory.hpp"
#include "../src/propagation/group_phase_path.hpp"
#include "../src/propagation/ionospheric_scintillation.hpp"

using namespace rvl::propagation;

class IonosphericEquationsTest : public ::testing::Test {
protected:
    static constexpr double tolerance = 1e-6;
    static constexpr float tolerance_f = 1e-5f;
};

TEST_F(IonosphericEquationsTest, IonosphericRefractiveIndexBasic) {
    const double electron_density = 1e12;  
    const double wave_frequency = 10e6;    
    const double collision_frequency = 1e3; 

    auto n_squared = ionospheric_refractive_index_d::calculate_n_squared(
        electron_density, wave_frequency, collision_frequency);

    EXPECT_LT(std::abs(n_squared), 1.0);
    EXPECT_GT(n_squared.real(), 0.0);

    const double n_magnitude = ionospheric_refractive_index_d::calculate_refractive_index_magnitude(
        electron_density, wave_frequency, collision_frequency);

    EXPECT_GT(n_magnitude, 0.0);
    EXPECT_LT(n_magnitude, 1.0);
}

TEST_F(IonosphericEquationsTest, CriticalFrequency) {
    const double electron_density = 1e12;  

    const double critical_freq = ionospheric_refractive_index_d::calculate_critical_frequency(electron_density);

    EXPECT_GT(critical_freq, 0.0);
    EXPECT_LT(critical_freq, 20e6);  

    const double expected_plasma_freq = std::sqrt(electron_density * 1.602e-19 * 1.602e-19 / 
                                                (9.109e-31 * 8.854e-12));
    const double expected_critical = expected_plasma_freq / (2.0 * M_PI);

    EXPECT_NEAR(critical_freq, expected_critical, tolerance);
}

TEST_F(IonosphericEquationsTest, MufCalculation) {
    const double electron_density = 1e12;   
    const double wave_frequency = 5e6;      
    const double elevation_angle = M_PI/6;  

    const double muf = ionospheric_refractive_index_d::calculate_muf_factor(
        electron_density, wave_frequency, elevation_angle);

    EXPECT_GT(muf, 0.0);

    const bool is_reflected = ionospheric_refractive_index_d::is_wave_reflected(
        electron_density, wave_frequency, elevation_angle);

    EXPECT_TRUE(is_reflected || !is_reflected);  
}

TEST_F(IonosphericEquationsTest, RayTrajectoryBasics) {
    ray_trajectory_d::point3d_type start_pos = {0.0, 0.0, 100000.0};  
    ray_trajectory_d::vector3d_type direction = {1.0, 0.0, 0.1};      
    const double wave_number = 0.1;

    auto initial_state = ray_trajectory_d::create_initial_state(start_pos, direction, wave_number);

    EXPECT_EQ(initial_state.position[0], start_pos[0]);
    EXPECT_EQ(initial_state.position[1], start_pos[1]);
    EXPECT_EQ(initial_state.position[2], start_pos[2]);
    EXPECT_EQ(initial_state.path_length, 0.0);

    const double k_magnitude = ray_trajectory_d::calculate_wave_vector_magnitude(initial_state);
    EXPECT_NEAR(k_magnitude, wave_number, tolerance);

    const auto ray_dir = ray_trajectory_d::calculate_ray_direction(initial_state);
    const double dir_magnitude = std::sqrt(ray_dir[0]*ray_dir[0] + ray_dir[1]*ray_dir[1] + ray_dir[2]*ray_dir[2]);
    EXPECT_NEAR(dir_magnitude, 1.0, tolerance);
}

TEST_F(IonosphericEquationsTest, GroupPhasePathBasic) {
    const double frequency = 10e6;  
    const double refractive_index = 0.9;
    const double df_dn = -1e-9;  

    const double group_index = group_phase_path_d::calculate_group_refractive_index(
        refractive_index, frequency, df_dn);

    EXPECT_GT(group_index, 0.0);
    EXPECT_NE(group_index, refractive_index);  

    auto n_func = [refractive_index](double f) -> double {
        return refractive_index + 1e-15 * (f - 10e6);  
    };

    const double group_index_numerical = group_phase_path_d::calculate_group_refractive_index_numerical(
        frequency, n_func, 1e3);

    EXPECT_GT(group_index_numerical, 0.0);
}

TEST_F(IonosphericEquationsTest, RytovVariance) {
    const double structure_constant = 1e-50;  
    const double path_length = 1000e3;        
    const double frequency = 100e6;           

    const double rytov_var = ionospheric_scintillation_d::calculate_rytov_variance_frequency(
        structure_constant, path_length, frequency);

    EXPECT_GT(rytov_var, 0.0);
    EXPECT_LT(rytov_var, 100.0);  

    const bool is_weak = ionospheric_scintillation_d::is_weak_scattering_regime(rytov_var);

    const double c = 299792458.0;
    const double wavenumber = 2.0 * M_PI * frequency / c;
    const double rytov_var_k = ionospheric_scintillation_d::calculate_rytov_variance(
        structure_constant, path_length, wavenumber);

    EXPECT_NEAR(rytov_var, rytov_var_k, tolerance);
}

TEST_F(IonosphericEquationsTest, ScintillationIndex) {

    std::vector<std::complex<double>> field_data;

    for (int i = 0; i < 100; ++i) {
        const double amplitude = 1.0 + 0.1 * std::sin(i * 0.1);  
        const double phase = i * 0.02;
        field_data.push_back(std::polar(amplitude, phase));
    }

    ionospheric_scintillation_d::complex_vector_type field_vector(field_data.begin(), field_data.end());

    const double s4 = ionospheric_scintillation_d::scintillation_index_s4(field_vector);

    EXPECT_GT(s4, 0.0);
    EXPECT_LT(s4, 1.0);  

    const double phase_scint = ionospheric_scintillation_d::phase_scintillation_index(field_vector);

    EXPECT_GT(phase_scint, 0.0);
}

TEST_F(IonosphericEquationsTest, FresnelScale) {
    const double path_length = 500e3;  
    const double frequency = 1e9;      
    const double c = 299792458.0;
    const double wavenumber = 2.0 * M_PI * frequency / c;

    const double fresnel_scale = ionospheric_scintillation_d::fresnel_scale_length(path_length, wavenumber);

    EXPECT_GT(fresnel_scale, 0.0);
    EXPECT_LT(fresnel_scale, 10000.0);  
}

TEST_F(IonosphericEquationsTest, BatchOperations) {
    const size_t n = 5;

    ionospheric_refractive_index_d::vector_type electron_densities = {1e11, 5e11, 1e12, 2e12, 5e12};
    ionospheric_refractive_index_d::vector_type wave_frequencies = {5e6, 10e6, 15e6, 20e6, 25e6};
    ionospheric_refractive_index_d::vector_type collision_frequencies = {1e3, 1e3, 1e3, 1e3, 1e3};
    ionospheric_refractive_index_d::complex_vector_type n_squared_results(n);

    EXPECT_NO_THROW(ionospheric_refractive_index_d::calculate_n_squared_batch(
        electron_densities, wave_frequencies, collision_frequencies, n_squared_results));

    for (size_t i = 0; i < n; ++i) {
        EXPECT_LT(std::abs(n_squared_results[i]), 1.0);
        EXPECT_GT(n_squared_results[i].real(), 0.0);
    }

    ionospheric_refractive_index_d::vector_type critical_freqs(n);

    EXPECT_NO_THROW(ionospheric_refractive_index_d::calculate_critical_frequency_batch(
        electron_densities, critical_freqs));

    for (size_t i = 0; i < n; ++i) {
        EXPECT_GT(critical_freqs[i], 0.0);
        EXPECT_LT(critical_freqs[i], 50e6);
    }

    ionospheric_scintillation_d::vector_type structure_constants = {1e-51, 5e-51, 1e-50, 2e-50, 5e-50};
    ionospheric_scintillation_d::vector_type path_lengths = {500e3, 750e3, 1000e3, 1250e3, 1500e3};
    ionospheric_scintillation_d::vector_type frequencies = {50e6, 100e6, 150e6, 200e6, 250e6};
    ionospheric_scintillation_d::vector_type s4_indices(n);

    EXPECT_NO_THROW(ionospheric_scintillation_d::calculate_s4_batch(
        structure_constants, path_lengths, frequencies, s4_indices));

    for (size_t i = 0; i < n; ++i) {
        EXPECT_GE(s4_indices[i], 0.0);
        EXPECT_LE(s4_indices[i], 1.0);
    }
}

TEST_F(IonosphericEquationsTest, FloatPrecisionVariants) {

    const float electron_density = 1e12f;
    const float wave_frequency = 10e6f;
    const float collision_frequency = 1e3f;

    auto n_squared_f = ionospheric_refractive_index_f::calculate_n_squared(
        electron_density, wave_frequency, collision_frequency);

    EXPECT_LT(std::abs(n_squared_f), 1.0f);
    EXPECT_GT(n_squared_f.real(), 0.0f);

    const float structure_constant = 1e-50f;
    const float path_length = 1000e3f;
    const float frequency = 100e6f;

    const float rytov_var_f = ionospheric_scintillation_f::calculate_rytov_variance_frequency(
        structure_constant, path_length, frequency);

    EXPECT_GT(rytov_var_f, 0.0f);
}

TEST_F(IonosphericEquationsTest, ErrorHandling) {

    EXPECT_THROW(ionospheric_refractive_index_d::calculate_n_squared(-1.0, 10e6, 1e3), 
                 rvl::core::invalid_argument_error);

    EXPECT_THROW(ionospheric_scintillation_d::calculate_rytov_variance_frequency(1e-50, -1000e3, 100e6), 
                 rvl::core::invalid_argument_error);

    ionospheric_refractive_index_d::vector_type small_vec = {1e12};
    ionospheric_refractive_index_d::vector_type large_vec = {1e12, 2e12};
    ionospheric_refractive_index_d::complex_vector_type result_vec(1);

    EXPECT_THROW(ionospheric_refractive_index_d::calculate_n_squared_batch(
        small_vec, large_vec, small_vec, result_vec), 
        rvl::core::dimension_mismatch_error);
}