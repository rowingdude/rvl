#!/bin/bash

# RadioVectorLib Complete Test Suite Runner

echo "========================================"
echo "RadioVectorLib Complete Test Suite"
echo "========================================"

# Compile flags
CXX=g++
CXXFLAGS="-std=c++20 -O3 -march=native -I. -Iinclude -lm"

# Arrays to track results
declare -a test_names
declare -a test_results

# Function to compile and run a test
run_test() {
    local test_file=$1
    local test_name=$(basename $test_file .cpp)
    local exe_name="build_test_${test_name}"
    
    echo -n "Running ${test_name}... "
    
    # Compile
    if $CXX $CXXFLAGS $test_file -o $exe_name 2>/dev/null; then
        # Run
        if ./$exe_name > /tmp/${test_name}.log 2>&1; then
            echo "PASSED"
            test_names+=("$test_name")
            test_results+=("PASSED")
        else
            echo "FAILED"
            test_names+=("$test_name")
            test_results+=("FAILED")
            echo "  Error output:"
            tail -20 /tmp/${test_name}.log | sed 's/^/    /'
        fi
        rm -f $exe_name
    else
        echo "COMPILATION FAILED"
        test_names+=("$test_name")
        test_results+=("COMPILE_FAIL")
    fi
}

# Find and run all test files
for test_file in tests/test_*.cpp; do
    if [ -f "$test_file" ]; then
        run_test "$test_file"
    fi
done

# Print summary
echo ""
echo "========================================"
echo "Test Summary"
echo "========================================"

passed=0
failed=0
compile_failed=0

for i in "${!test_names[@]}"; do
    printf "%-40s %s\n" "${test_names[$i]}" "${test_results[$i]}"
    case "${test_results[$i]}" in
        "PASSED") ((passed++)) ;;
        "FAILED") ((failed++)) ;;
        "COMPILE_FAIL") ((compile_failed++)) ;;
    esac
done

echo "========================================"
echo "Total tests: ${#test_names[@]}"
echo "Passed: $passed"
echo "Failed: $failed"
echo "Compilation failed: $compile_failed"
echo "========================================"

# Exit with error if any tests failed
if [ $failed -gt 0 ] || [ $compile_failed -gt 0 ]; then
    exit 1
else
    exit 0
fi