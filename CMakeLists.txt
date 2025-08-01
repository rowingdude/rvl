cmake_minimum_required(VERSION 3.20)
project(RadioVectorLib VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(RVL_BUILD_TESTS "Build unit tests" ON)
option(RVL_BUILD_EXAMPLES "Build examples" ON)
option(RVL_BUILD_BENCHMARKS "Build benchmarks" ON)
option(RVL_ENABLE_OPENMP "Enable OpenMP support" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_library(radio_vector_lib INTERFACE)
add_library(rvl::rvl ALIAS radio_vector_lib)

target_include_directories(radio_vector_lib
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(radio_vector_lib INTERFACE
        -march=native
        -Wall
        -Wextra
        -Wpedantic
        -Wno-unused-parameter
        $<$<CONFIG:Release>:-O3 -DNDEBUG>
        $<$<CONFIG:Debug>:-O0 -g -DRVL_DEBUG>
        $<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
    )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(radio_vector_lib INTERFACE
        /arch:AVX2
        /W4
        /wd4101  # Suppress unused variable warnings  
        $<$<CONFIG:Release>:/O2 /DNDEBUG>
        $<$<CONFIG:Debug>:/Od /Zi /DRVL_DEBUG>
        $<$<CONFIG:RelWithDebInfo>:/O2 /Zi /DNDEBUG>
    )
endif()

if(RVL_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(radio_vector_lib INTERFACE OpenMP::OpenMP_CXX)
        target_compile_definitions(radio_vector_lib INTERFACE RVL_USE_OPENMP)
    endif()
endif()

find_package(Threads REQUIRED)
target_link_libraries(radio_vector_lib INTERFACE Threads::Threads)

if(RVL_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
    
    # Add convenience target to build and run all tests
    add_custom_target(test_all
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Building and running all RadioVectorLib tests"
        VERBATIM
    )
endif()

if(RVL_BUILD_EXAMPLES AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples")
    add_subdirectory(examples)
    message(STATUS "Building examples")
elseif(RVL_BUILD_EXAMPLES)
    message(STATUS "Examples directory not found, skipping examples")
endif()

if(RVL_BUILD_BENCHMARKS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks")
    add_subdirectory(benchmarks)
    message(STATUS "Building benchmarks")
elseif(RVL_BUILD_BENCHMARKS)
    message(STATUS "Benchmarks directory not found, skipping benchmarks")
endif()

include(GNUInstallDirs)
install(TARGETS radio_vector_lib
    EXPORT RadioVectorLibTargets
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY src/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rvl
    FILES_MATCHING PATTERN "*.hpp"
)

install(EXPORT RadioVectorLibTargets
    FILE RadioVectorLibTargets.cmake
    NAMESPACE rvl::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RadioVectorLib
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RadioVectorLibConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/RadioVectorLibConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RadioVectorLib
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/RadioVectorLibConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/RadioVectorLibConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/RadioVectorLibConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RadioVectorLib
)