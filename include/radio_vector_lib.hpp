#ifndef RADIO_VECTOR_LIB_HPP
#define RADIO_VECTOR_LIB_HPP

#include "rvl_core.hpp"

#define RVL_VERSION_MAJOR 1
#define RVL_VERSION_MINOR 0
#define RVL_VERSION_PATCH 0

#define RVL_VERSION_STRING "1.0.0"

namespace rvl {

struct version_info {
    static constexpr int major = RVL_VERSION_MAJOR;
    static constexpr int minor = RVL_VERSION_MINOR;
    static constexpr int patch = RVL_VERSION_PATCH;
    static constexpr const char* string = RVL_VERSION_STRING;
};

} 

#endif 