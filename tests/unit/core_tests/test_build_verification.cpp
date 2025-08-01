#include <iostream>
#include <iomanip>
#include "radio_vector_lib.hpp"
#include "../src/core/memory_allocator.hpp"
#include "../src/core/aligned_buffer.hpp"
#include "../src/core/simd_detection.hpp"
#include "../src/core/constants.hpp"
#include "../src/core/units.hpp"
#include "../src/core/error.hpp"
#include "../src/core/logging.hpp"

using namespace rvl;

int main() {
    std::cout << "RadioVectorLib Build Verification Test\n";
    std::cout << "=====================================\n\n";

    std::cout << "Library Version: " << version_info::string << "\n";
    std::cout << "  Major: " << version_info::major << "\n";
    std::cout << "  Minor: " << version_info::minor << "\n";
    std::cout << "  Patch: " << version_info::patch << "\n\n";

    std::cout << "CPU SIMD Features:\n";
    std::cout << "  Best available: " << core::simd::cpu_features::instance().best_available_name() << "\n";
    std::cout << "  SSE:     " << (core::simd::has_sse() ? "Yes" : "No") << "\n";
    std::cout << "  SSE2:    " << (core::simd::has_sse2() ? "Yes" : "No") << "\n";
    std::cout << "  SSE3:    " << (core::simd::has_sse3() ? "Yes" : "No") << "\n";
    std::cout << "  SSSE3:   " << (core::simd::has_ssse3() ? "Yes" : "No") << "\n";
    std::cout << "  SSE4.1:  " << (core::simd::has_sse41() ? "Yes" : "No") << "\n";
    std::cout << "  SSE4.2:  " << (core::simd::has_sse42() ? "Yes" : "No") << "\n";
    std::cout << "  AVX:     " << (core::simd::has_avx() ? "Yes" : "No") << "\n";
    std::cout << "  AVX2:    " << (core::simd::has_avx2() ? "Yes" : "No") << "\n";
    std::cout << "  FMA3:    " << (core::simd::has_fma3() ? "Yes" : "No") << "\n";
    std::cout << "  AVX512F: " << (core::simd::has_avx512f() ? "Yes" : "No") << "\n";
    std::cout << "  NEON:    " << (core::simd::has_neon() ? "Yes" : "No") << "\n\n";

    std::cout << "Testing Memory Allocation:\n";
    try {
        core::memory::aligned_buffer<float> buffer(1024);
        std::cout << "  ✓ Aligned buffer allocation (1024 floats)\n";
        std::cout << "    Size: " << buffer.size() << "\n";
        std::cout << "    Capacity: " << buffer.capacity() << "\n";
        std::cout << "    Alignment: " << (reinterpret_cast<uintptr_t>(buffer.data()) % 64 == 0 ? "64-byte aligned" : "NOT aligned!") << "\n";

        core::memory::simd_vector<double> vec;
        vec.resize(512);
        std::cout << "  ✓ SIMD vector allocation (512 doubles)\n";
    } catch (const std::exception& e) {
        std::cout << "   Memory allocation failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "\nTesting Constants and Units:\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  Speed of light: " << constants::physical_d::c << " m/s\n";
    std::cout << "  10 dBm to Watts: " << units::dbm_to_watts(10.0) << " W\n";
    std::cout << "  1 GHz wavelength: " << units::frequency_to_wavelength(1e9) << " m\n";
    std::cout << "  VSWR 2.0 return loss: " << units::vswr_to_return_loss_db(2.0) << " dB\n";

    std::cout << "\nTesting Logging:\n";
    core::get_logger().set_level(core::log_level::TRACE);
    RVL_LOG_TRACE("This is a trace message");
    RVL_LOG_DEBUG("This is a debug message");
    RVL_LOG_INFO("This is an info message");
    RVL_LOG_WARNING("This is a warning message");

    std::cout << "\nTesting Error Handling:\n";
    try {
        throw core::invalid_argument_error("Test error");
    } catch (const core::rvl_error& e) {
        std::cout << "  ✓ Caught error: " << e.what() << "\n";
        std::cout << "    Error code: " << static_cast<int>(e.code()) << "\n";
    }

    std::cout << "\nBuild verification PASSED!\n";
    return 0;
}