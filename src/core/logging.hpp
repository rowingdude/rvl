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

#ifndef RVL_CORE_LOGGING_HPP
#define RVL_CORE_LOGGING_HPP

#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <mutex>

namespace rvl {
namespace core {

enum class log_level {
    TRACE = 0,
    DEBUG = 1,
    INFO = 2,
    WARNING = 3,
    ERROR = 4,
    CRITICAL = 5,
    OFF = 6
};

class logger {
private:
    log_level min_level_ = log_level::INFO;
    mutable std::mutex mutex_;
    std::ostream* output_ = &std::cerr;

    logger() = default;

    const char* level_string(log_level level) const {
        switch (level) {
            case log_level::TRACE:    return "TRACE";
            case log_level::DEBUG:    return "DEBUG";
            case log_level::INFO:     return "INFO ";
            case log_level::WARNING:  return "WARN ";
            case log_level::ERROR:    return "ERROR";
            case log_level::CRITICAL: return "CRIT ";
            default:                  return "?????";
        }
    }

    void log_impl(log_level level, const std::string& message) {
        if (level < min_level_) return;

        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;

        std::lock_guard<std::mutex> lock(mutex_);
        *output_ << "[" << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S")
                 << "." << std::setfill('0') << std::setw(3) << ms.count()
                 << "] [" << level_string(level) << "] "
                 << message << std::endl;
    }

public:
    static logger& instance() {
        static logger instance;
        return instance;
    }

    void set_level(log_level level) {
        std::lock_guard<std::mutex> lock(mutex_);
        min_level_ = level;
    }

    void set_output(std::ostream& output) {
        std::lock_guard<std::mutex> lock(mutex_);
        output_ = &output;
    }

    log_level get_level() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return min_level_;
    }

    template<typename... Args>
    void log(log_level level, Args&&... args) {
        if (level < min_level_) return;

        std::ostringstream oss;
        ((oss << std::forward<Args>(args)), ...);
        log_impl(level, oss.str());
    }

    template<typename... Args>
    void trace(Args&&... args) { log(log_level::TRACE, std::forward<Args>(args)...); }

    template<typename... Args>
    void debug(Args&&... args) { log(log_level::DEBUG, std::forward<Args>(args)...); }

    template<typename... Args>
    void info(Args&&... args) { log(log_level::INFO, std::forward<Args>(args)...); }

    template<typename... Args>
    void warning(Args&&... args) { log(log_level::WARNING, std::forward<Args>(args)...); }

    template<typename... Args>
    void error(Args&&... args) { log(log_level::ERROR, std::forward<Args>(args)...); }

    template<typename... Args>
    void critical(Args&&... args) { log(log_level::CRITICAL, std::forward<Args>(args)...); }
};

inline logger& get_logger() { return logger::instance(); }

#define RVL_LOG_TRACE(...) ::rvl::core::get_logger().trace(__VA_ARGS__)
#define RVL_LOG_DEBUG(...) ::rvl::core::get_logger().debug(__VA_ARGS__)
#define RVL_LOG_INFO(...) ::rvl::core::get_logger().info(__VA_ARGS__)
#define RVL_LOG_WARNING(...) ::rvl::core::get_logger().warning(__VA_ARGS__)
#define RVL_LOG_ERROR(...) ::rvl::core::get_logger().error(__VA_ARGS__)
#define RVL_LOG_CRITICAL(...) ::rvl::core::get_logger().critical(__VA_ARGS__)

#ifdef RVL_ENABLE_TRACE
    #define RVL_TRACE(...) RVL_LOG_TRACE(__VA_ARGS__)
#else
    #define RVL_TRACE(...) ((void)0)
#endif

} 
} 

#endif 