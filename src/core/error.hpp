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

#ifndef RVL_CORE_ERROR_HPP
#define RVL_CORE_ERROR_HPP

#include <stdexcept>
#include <string>
#include <sstream>
#include <source_location>

namespace rvl {
namespace core {

enum class error_code {
    SUCCESS = 0,
    INVALID_ARGUMENT,
    OUT_OF_RANGE,
    ALLOCATION_FAILED,
    DIMENSION_MISMATCH,
    NUMERICAL_ERROR,
    NOT_IMPLEMENTED,
    IO_ERROR,
    INVALID_STATE,
    CONVERGENCE_FAILED
};

class rvl_error : public std::runtime_error {
private:
    error_code code_;
    std::source_location location_;

public:
    rvl_error(error_code code, const std::string& message,
              const std::source_location& loc = std::source_location::current())
        : std::runtime_error(format_error(code, message, loc)),
          code_(code),
          location_(loc) {}

    error_code code() const noexcept { return code_; }
    const std::source_location& location() const noexcept { return location_; }

private:
    static std::string format_error(error_code code, const std::string& message,
                                   const std::source_location& loc) {
        std::ostringstream oss;
        oss << "[RVL Error " << static_cast<int>(code) << "] "
            << message
            << " (at " << loc.file_name() << ":" << loc.line() << ")";
        return oss.str();
    }
};

class invalid_argument_error : public rvl_error {
public:
    explicit invalid_argument_error(const std::string& message,
                                  const std::source_location& loc = std::source_location::current())
        : rvl_error(error_code::INVALID_ARGUMENT, message, loc) {}
};

class out_of_range_error : public rvl_error {
public:
    explicit out_of_range_error(const std::string& message,
                               const std::source_location& loc = std::source_location::current())
        : rvl_error(error_code::OUT_OF_RANGE, message, loc) {}
};

class allocation_error : public rvl_error {
public:
    explicit allocation_error(const std::string& message,
                            const std::source_location& loc = std::source_location::current())
        : rvl_error(error_code::ALLOCATION_FAILED, message, loc) {}
};

class dimension_mismatch_error : public rvl_error {
public:
    explicit dimension_mismatch_error(const std::string& message,
                                    const std::source_location& loc = std::source_location::current())
        : rvl_error(error_code::DIMENSION_MISMATCH, message, loc) {}
};

class numerical_error : public rvl_error {
public:
    explicit numerical_error(const std::string& message,
                           const std::source_location& loc = std::source_location::current())
        : rvl_error(error_code::NUMERICAL_ERROR, message, loc) {}
};

class not_implemented_error : public rvl_error {
public:
    explicit not_implemented_error(const std::string& message = "Feature not yet implemented",
                                 const std::source_location& loc = std::source_location::current())
        : rvl_error(error_code::NOT_IMPLEMENTED, message, loc) {}
};

#ifdef RVL_ENABLE_CHECKS
    #define RVL_ASSERT(condition, message) \
        do { \
            if (!(condition)) { \
                throw rvl::core::rvl_error(rvl::core::error_code::INVALID_STATE, \
                                          std::string("Assertion failed: ") + message); \
            } \
        } while(0)

    #define RVL_REQUIRE(condition, error_type, message) \
        do { \
            if (!(condition)) { \
                throw error_type(message); \
            } \
        } while(0)
#else
    #define RVL_ASSERT(condition, message) ((void)0)
    #define RVL_REQUIRE(condition, error_type, message) ((void)0)
#endif

#define RVL_NOT_IMPLEMENTED() \
    throw rvl::core::not_implemented_error()

#define RVL_NOT_IMPLEMENTED_MSG(msg) \
    throw rvl::core::not_implemented_error(msg)

template<typename T>
inline void check_finite(T value, const std::string& name) {
    if (!std::isfinite(value)) {
        throw numerical_error(name + " must be finite");
    }
}

template<typename T>
inline void check_positive(T value, const std::string& name) {
    if (value <= T(0)) {
        throw invalid_argument_error(name + " must be positive");
    }
}

template<typename T>
inline void check_non_negative(T value, const std::string& name) {
    if (value < T(0)) {
        throw invalid_argument_error(name + " must be non-negative");
    }
}

template<typename T>
inline void check_range(T value, T min, T max, const std::string& name) {
    if (value < min || value > max) {
        std::ostringstream oss;
        oss << name << " must be in range [" << min << ", " << max << "]";
        throw out_of_range_error(oss.str());
    }
}

} 
} 

#endif 