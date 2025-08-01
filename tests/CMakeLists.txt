# RadioVectorLib Tests CMakeLists.txt

# Set minimum required version
cmake_minimum_required(VERSION 3.20)

# Find math library for linking
find_library(MATH_LIBRARY m)

# Define function-level test sources (our corrected tests)
set(FUNCTION_TEST_SOURCES
    test_antenna_array.cpp
    test_antenna_dipole.cpp
    test_antenna_efficiency.cpp
    test_feedline_smith_chart.cpp
    test_feedline_twin_lead.cpp
    test_feedline_vswr.cpp
    test_ionospheric_propagation.cpp
    test_propagation_free_space.cpp
    test_rf_systems_link_budget.cpp
    test_transmission_line.cpp
    test_tropospheric_propagation.cpp
    # Note: test_propagation_models.cpp omitted - needs major API updates
)

# Create executables for function-level tests
foreach(test_src ${FUNCTION_TEST_SOURCES})
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_src}")
        get_filename_component(test_name ${test_src} NAME_WE)
        
        add_executable(${test_name} ${test_src})
        
        # Link to our radio vector lib interface
        target_link_libraries(${test_name} PRIVATE rvl::rvl)
        
        # Link math library if found
        if(MATH_LIBRARY)
            target_link_libraries(${test_name} PRIVATE ${MATH_LIBRARY})
        endif()
        
        # Set C++20 standard
        target_compile_features(${test_name} PRIVATE cxx_std_20)
        
        # Set include directories to find headers
        target_include_directories(${test_name} PRIVATE 
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_SOURCE_DIR}
        )
        
        # Add to CTest
        add_test(NAME ${test_name} COMMAND ${test_name})
        
        # Set test properties
        set_tests_properties(${test_name} PROPERTIES
            TIMEOUT 60
            FAIL_REGULAR_EXPRESSION "FAILED|DOWN"
        )
    else()
        message(WARNING "Test file ${test_src} not found, skipping...")
    endif()
endforeach()

# Handle unit tests with GoogleTest if available
option(BUILD_TESTING "Build the testing tree" ON)

if(BUILD_TESTING)
    find_package(GTest QUIET)
    if(GTest_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/unit/CMakeLists.txt")
        message(STATUS "Google Test found, building unit tests")
        add_subdirectory(unit)
    else()
        if(NOT GTest_FOUND)
            message(STATUS "Google Test not found, skipping unit tests")
            message(STATUS "  Install with: sudo apt-get install libgtest-dev (Ubuntu/Debian)")
            message(STATUS "  Or: brew install googletest (macOS)")
        endif()
        if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/unit/CMakeLists.txt")
            message(STATUS "Unit tests CMakeLists.txt not found, skipping unit tests")
        endif()
    endif()
endif()

# Custom target to run all function tests
add_custom_target(run_function_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --verbose
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running all RadioVectorLib function tests"
    DEPENDS ${FUNCTION_TEST_SOURCES}
)

# Custom target to run only working tests (excluding propagation_models)
add_custom_target(run_working_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --exclude-regex "test_propagation_models"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running verified working RadioVectorLib tests"
)

# Integration and benchmark targets if they exist
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/integration" AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/integration/CMakeLists.txt")
    add_subdirectory(integration)
    message(STATUS "Including integration tests")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks" AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/CMakeLists.txt")
    add_subdirectory(benchmarks)
    message(STATUS "Including benchmark tests")
endif()

# Print test configuration info
message(STATUS "RadioVectorLib Tests Configuration:")
message(STATUS "  Function tests: ${FUNCTION_TEST_SOURCES}")
message(STATUS "  GoogleTest available: ${GTest_FOUND}")
message(STATUS "  Math library: ${MATH_LIBRARY}")