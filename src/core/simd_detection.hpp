/*
================================================================================
||                                                                            ||
||                              RadioVectorLibrary                            ||
||                                   (RVL)                                    ||
||                                                                            ||
||                                                                            ||
||                         Developer: Benjamin Cance, KC8BWS                  ||
||                         Email: kc8bws@kc8bws.com                           ||
||                         MIT License                                        ||
||                                                                            ||
||                                                                            ||
================================================================================
*/

#ifndef RVL_CORE_SIMD_DETECTION_HPP
#define RVL_CORE_SIMD_DETECTION_HPP

#include <cstdint>

#ifdef _MSC_VER
    #include <intrin.h>
#else
    #include <cpuid.h>
#endif

namespace rvl {
namespace core {
namespace simd {

enum class instruction_set : uint32_t {
    NONE      = 0,
    SSE       = 1 << 0,
    SSE2      = 1 << 1,
    SSE3      = 1 << 2,
    SSSE3     = 1 << 3,
    SSE41     = 1 << 4,
    SSE42     = 1 << 5,
    AVX       = 1 << 6,
    AVX2      = 1 << 7,
    FMA3      = 1 << 8,
    AVX512F   = 1 << 9,
    AVX512DQ  = 1 << 10,
    AVX512BW  = 1 << 11,
    AVX512VL  = 1 << 12,
    NEON      = 1 << 13,
    SVE       = 1 << 14,
    SVE2      = 1 << 15
};

inline instruction_set operator|(instruction_set a, instruction_set b) {
    return static_cast<instruction_set>(static_cast<uint32_t>(a) | static_cast<uint32_t>(b));
}

inline instruction_set operator&(instruction_set a, instruction_set b) {
    return static_cast<instruction_set>(static_cast<uint32_t>(a) & static_cast<uint32_t>(b));
}

inline bool has_instruction_set(instruction_set available, instruction_set query) {
    return (available & query) == query;
}

class cpu_features {
private:
    instruction_set available_sets_ = instruction_set::NONE;

    cpu_features() {
        detect_features();
    }

    void detect_features() {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
        int info[4];

#ifdef _MSC_VER
        __cpuid(info, 0);
        int max_id = info[0];

        if (max_id >= 1) {
            __cpuid(info, 1);
            if (info[3] & (1 << 25)) available_sets_ = available_sets_ | instruction_set::SSE;
            if (info[3] & (1 << 26)) available_sets_ = available_sets_ | instruction_set::SSE2;
            if (info[2] & (1 << 0))  available_sets_ = available_sets_ | instruction_set::SSE3;
            if (info[2] & (1 << 9))  available_sets_ = available_sets_ | instruction_set::SSSE3;
            if (info[2] & (1 << 19)) available_sets_ = available_sets_ | instruction_set::SSE41;
            if (info[2] & (1 << 20)) available_sets_ = available_sets_ | instruction_set::SSE42;
            if (info[2] & (1 << 28)) available_sets_ = available_sets_ | instruction_set::AVX;
            if (info[2] & (1 << 12)) available_sets_ = available_sets_ | instruction_set::FMA3;
        }

        if (max_id >= 7) {
            __cpuidex(info, 7, 0);
            if (info[1] & (1 << 5))  available_sets_ = available_sets_ | instruction_set::AVX2;
            if (info[1] & (1 << 16)) available_sets_ = available_sets_ | instruction_set::AVX512F;
            if (info[1] & (1 << 17)) available_sets_ = available_sets_ | instruction_set::AVX512DQ;
            if (info[1] & (1 << 30)) available_sets_ = available_sets_ | instruction_set::AVX512BW;
            if (info[1] & (1 << 31)) available_sets_ = available_sets_ | instruction_set::AVX512VL;
        }
#else
        __cpuid(0, info[0], info[1], info[2], info[3]);
        int max_id = info[0];

        if (max_id >= 1) {
            __cpuid(1, info[0], info[1], info[2], info[3]);
            if (info[3] & (1 << 25)) available_sets_ = available_sets_ | instruction_set::SSE;
            if (info[3] & (1 << 26)) available_sets_ = available_sets_ | instruction_set::SSE2;
            if (info[2] & (1 << 0))  available_sets_ = available_sets_ | instruction_set::SSE3;
            if (info[2] & (1 << 9))  available_sets_ = available_sets_ | instruction_set::SSSE3;
            if (info[2] & (1 << 19)) available_sets_ = available_sets_ | instruction_set::SSE41;
            if (info[2] & (1 << 20)) available_sets_ = available_sets_ | instruction_set::SSE42;
            if (info[2] & (1 << 28)) available_sets_ = available_sets_ | instruction_set::AVX;
            if (info[2] & (1 << 12)) available_sets_ = available_sets_ | instruction_set::FMA3;
        }

        if (max_id >= 7) {
            __cpuid_count(7, 0, info[0], info[1], info[2], info[3]);
            if (info[1] & (1 << 5))  available_sets_ = available_sets_ | instruction_set::AVX2;
            if (info[1] & (1 << 16)) available_sets_ = available_sets_ | instruction_set::AVX512F;
            if (info[1] & (1 << 17)) available_sets_ = available_sets_ | instruction_set::AVX512DQ;
            if (info[1] & (1 << 30)) available_sets_ = available_sets_ | instruction_set::AVX512BW;
            if (info[1] & (1 << 31)) available_sets_ = available_sets_ | instruction_set::AVX512VL;
        }
#endif

#elif defined(__aarch64__) || defined(_M_ARM64)
        available_sets_ = available_sets_ | instruction_set::NEON;

#ifdef __ARM_FEATURE_SVE
        available_sets_ = available_sets_ | instruction_set::SVE;
#endif
#ifdef __ARM_FEATURE_SVE2
        available_sets_ = available_sets_ | instruction_set::SVE2;
#endif
#endif
    }

public:
    static const cpu_features& instance() {
        static cpu_features features;
        return features;
    }

    bool has(instruction_set set) const {
        return has_instruction_set(available_sets_, set);
    }

    instruction_set available() const {
        return available_sets_;
    }

    const char* best_available_name() const {
        if (has(instruction_set::AVX512F)) return "AVX512";
        if (has(instruction_set::AVX2)) return "AVX2";
        if (has(instruction_set::AVX)) return "AVX";
        if (has(instruction_set::SSE42)) return "SSE4.2";
        if (has(instruction_set::SSE41)) return "SSE4.1";
        if (has(instruction_set::SSSE3)) return "SSSE3";
        if (has(instruction_set::SSE3)) return "SSE3";
        if (has(instruction_set::SSE2)) return "SSE2";
        if (has(instruction_set::SSE)) return "SSE";
        if (has(instruction_set::SVE2)) return "SVE2";
        if (has(instruction_set::SVE)) return "SVE";
        if (has(instruction_set::NEON)) return "NEON";
        return "Scalar";
    }
};

inline bool has_sse() { return cpu_features::instance().has(instruction_set::SSE); }
inline bool has_sse2() { return cpu_features::instance().has(instruction_set::SSE2); }
inline bool has_sse3() { return cpu_features::instance().has(instruction_set::SSE3); }
inline bool has_ssse3() { return cpu_features::instance().has(instruction_set::SSSE3); }
inline bool has_sse41() { return cpu_features::instance().has(instruction_set::SSE41); }
inline bool has_sse42() { return cpu_features::instance().has(instruction_set::SSE42); }
inline bool has_avx() { return cpu_features::instance().has(instruction_set::AVX); }
inline bool has_avx2() { return cpu_features::instance().has(instruction_set::AVX2); }
inline bool has_fma3() { return cpu_features::instance().has(instruction_set::FMA3); }
inline bool has_avx512f() { return cpu_features::instance().has(instruction_set::AVX512F); }
inline bool has_neon() { return cpu_features::instance().has(instruction_set::NEON); }

} 
} 
} 

#endif 