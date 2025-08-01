#ifndef RVL_TEST_UTILS_HPP
#define RVL_TEST_UTILS_HPP

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <complex>

namespace rvl_test {

template<typename T>
bool approx_equal(T a, T b, T tolerance = T(1e-6)) {
    return std::abs(a - b) < tolerance;
}

template<typename T>
bool approx_equal(std::complex<T> a, std::complex<T> b, T tolerance = T(1e-6)) {
    return std::abs(a - b) < tolerance;
}

class TestResult {
public:
    std::string function_name;
    bool passed;
    std::string error_message;
    
    TestResult(const std::string& name, bool pass, const std::string& error = "")
        : function_name(name), passed(pass), error_message(error) {}
};

class TestSuite {
private:
    std::vector<TestResult> results;
    std::string suite_name;
    
public:
    TestSuite(const std::string& name) : suite_name(name) {}
    
    void add_test(const std::string& function_name, bool passed, const std::string& error = "") {
        results.emplace_back(function_name, passed, error);
    }
    
    void print_summary() const {
        std::cout << "\n=== " << suite_name << " ===" << std::endl;
        for (const auto& result : results) {
            std::cout << result.function_name << ": " 
                     << (result.passed ? "UP" : "DOWN");
            if (!result.passed && !result.error_message.empty()) {
                std::cout << " (" << result.error_message << ")";
            }
            std::cout << std::endl;
        }
        
        int passed = 0;
        for (const auto& result : results) {
            if (result.passed) passed++;
        }
        
        std::cout << "\nTotal: " << passed << "/" << results.size() << " passed" << std::endl;
    }
    
    bool all_passed() const {
        for (const auto& result : results) {
            if (!result.passed) return false;
        }
        return true;
    }
};

} 

#endif 