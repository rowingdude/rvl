#include "test_utils.hpp"
#include "src/propagation/ionospheric_ray_tracing.hpp"
#include "src/propagation/ionospheric_refractive_index.hpp"

using namespace rvl;
using namespace rvl_test;

int main() {
    TestSuite suite("Ionospheric Propagation Functions");

    try {
        propagation::ionospheric_ray_tracing<double>::ionospheric_profile profile;
        profile.foF2 = 10e6;
        profile.hmF2 = 300.0;
        profile.foE = 3e6;
        profile.hmE = 110.0;
        profile.critical_frequency_factor = 1.0;
        
        double ne = propagation::ionospheric_ray_tracing<double>::calculate_electron_density(
            250.0, profile);
        bool passed = ne > 1e10 && ne < 1e13;
        suite.add_test("calculate_electron_density", passed);
    } catch (...) {
        suite.add_test("calculate_electron_density", false, "exception thrown");
    }
    
    try {
        propagation::ionospheric_ray_tracing<double>::ray_state state;
        state.position = {0.0, 0.0, 200.0};
        state.direction = {0.8, 0.0, 0.6};
        state.time = 0.0;
        
        propagation::ionospheric_ray_tracing<double>::ionospheric_profile profile;
        profile.foF2 = 10e6;
        profile.hmF2 = 300.0;
        
        auto new_state = propagation::ionospheric_ray_tracing<double>::integrate_ray_step(
            state, profile, 14e6, 1.0);
        bool passed = new_state.time > state.time;
        suite.add_test("integrate_ray_step", passed);
    } catch (...) {
        suite.add_test("integrate_ray_step", false, "exception thrown");
    }
    
    try {
        propagation::ionospheric_ray_tracing<double>::ray_parameters params;
        params.frequency_hz = 14e6;
        params.elevation_angle_deg = 30.0;
        params.azimuth_angle_deg = 0.0;
        params.max_height_km = 500.0;
        params.step_size_km = 1.0;
        
        propagation::ionospheric_ray_tracing<double>::ionospheric_profile profile;
        profile.foF2 = 10e6;
        profile.hmF2 = 300.0;
        
        auto trajectory = propagation::ionospheric_ray_tracing<double>::trace_ray(params, profile);
        bool passed = trajectory.points.size() > 0;
        suite.add_test("trace_ray", passed);
    } catch (...) {
        suite.add_test("trace_ray", false, "exception thrown");
    }
    
    try {
        double muf = propagation::ionospheric_ray_tracing<double>::calculate_muf(
            1000.0, 10e6, 300.0);
        bool passed = muf > 10e6 && muf < 50e6;
        suite.add_test("calculate_muf", passed);
    } catch (...) {
        suite.add_test("calculate_muf", false, "exception thrown");
    }
    
    try {
        double angle = propagation::ionospheric_ray_tracing<double>::calculate_critical_angle(
            14e6, 10e6);
        bool passed = angle >= 0.0 && angle <= 90.0;
        suite.add_test("calculate_critical_angle", passed);
    } catch (...) {
        suite.add_test("calculate_critical_angle", false, "exception thrown");
    }

    try {
        double n = propagation::ionospheric_refractive_index<double>::calculate_refractive_index(
            14e6, 1e12, 0.0, 0.0);
        bool passed = n > 0.0 && n <= 1.0;
        suite.add_test("calculate_refractive_index", passed);
    } catch (...) {
        suite.add_test("calculate_refractive_index", false, "exception thrown");
    }
    
    try {
        double fp = propagation::ionospheric_refractive_index<double>::calculate_plasma_frequency(1e12);
        bool passed = fp > 1e6 && fp < 20e6;
        suite.add_test("calculate_plasma_frequency", passed);
    } catch (...) {
        suite.add_test("calculate_plasma_frequency", false, "exception thrown");
    }
    
    try {
        double ne = propagation::ionospheric_refractive_index<double>::electron_density_from_plasma_frequency(10e6);
        bool passed = ne > 1e11 && ne < 1e13;
        suite.add_test("electron_density_from_plasma_frequency", passed);
    } catch (...) {
        suite.add_test("electron_density_from_plasma_frequency", false, "exception thrown");
    }
    
    try {
        auto indices = propagation::ionospheric_refractive_index<double>::calculate_complex_refractive_index(
            14e6, 1e12, 0.0, 1e-9);
        bool passed = indices.ordinary.real() > 0.0 && indices.extraordinary.real() > 0.0;
        suite.add_test("calculate_complex_refractive_index", passed);
    } catch (...) {
        suite.add_test("calculate_complex_refractive_index", false, "exception thrown");
    }
    
    try {
        double delay = propagation::ionospheric_refractive_index<double>::calculate_group_delay(
            14e6, 1e12, 1000.0);
        bool passed = delay > 0.0 && delay < 0.01;
        suite.add_test("calculate_group_delay", passed);
    } catch (...) {
        suite.add_test("calculate_group_delay", false, "exception thrown");
    }
    
    try {
        double vg = propagation::ionospheric_refractive_index<double>::calculate_group_velocity(
            14e6, 1e12);
        bool passed = vg > 0.0 && vg <= 3e8;
        suite.add_test("calculate_group_velocity", passed);
    } catch (...) {
        suite.add_test("calculate_group_velocity", false, "exception thrown");
    }
    
    try {
        bool crit = propagation::ionospheric_refractive_index<double>::is_critical_frequency(
            10e6, 1.24e12);
        suite.add_test("is_critical_frequency", crit);
    } catch (...) {
        suite.add_test("is_critical_frequency", false, "exception thrown");
    }
    
    suite.print_summary();
    return suite.all_passed() ? 0 : 1;
}