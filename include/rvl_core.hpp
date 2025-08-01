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

#ifndef RVL_CORE_HPP
#define RVL_CORE_HPP

#include <cstddef>
#include <cstdint>

namespace rvl {

    namespace core {
        namespace simd {
            namespace x86 {}
            namespace arm {}
        }
        namespace memory {}
    }

    namespace antenna {
        namespace patterns {}
        namespace arrays {}
    }

    namespace propagation {
        namespace ionospheric {}
        namespace tropospheric {}
        namespace ground_wave {}
    }

    namespace rf_systems {}

    namespace utils {}

    namespace detail {}

} 

#endif 