#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <iomanip>
#include <memory>
#include <mutex>
#include <thread>

// Define a macro to get the function name in a cross-platform way
#if defined(__GNUC__) || (defined(__MWERKS__) && (__MWERKS__ >= 0x3000)) || (defined(__ICC) && (__ICC >= 600)) || defined(__ghs__)
    #define LOG_FUNC_NAME __PRETTY_FUNCTION__
#elif defined(__DMC__) && (__DMC__ >= 0x810)
    #define LOG_FUNC_NAME __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
    #define LOG_FUNC_NAME __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600)) || (defined(__IBMCPP__) && (__IBMCPP__ >= 500))
    #define LOG_FUNC_NAME __FUNCTION__
#elif defined(__BORLANDC__) && (__BORLANDC__ >= 0x550)
    #define LOG_FUNC_NAME __FUNC__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
    #define LOG_FUNC_NAME __func__
#elif defined(__cplusplus) && (__cplusplus >= 201103)
    #define LOG_FUNC_NAME __func__
#else
    #define LOG_FUNC_NAME "(unknown function)"
#endif


/**
 * @enum LogLevel
 * @brief Defines the severity levels for log messages.
 */
enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR,
    CRITICAL
};

/**
 * @class Logger
 * @brief A thread-safe singleton logger for C++ applications.
 *
 * Provides a flexible logging framework with different levels, configurable outputs
 * (console and/or file), and detailed message formatting including timestamp,
 * log level, and source location.
 */
class Logger {
public:
    // Delete copy constructor and assignment operator for singleton
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    /**
     * @brief Get the singleton instance of the Logger.
     * @return Reference to the single Logger instance.
     */
    static Logger& get_instance () {
        static Logger instance;
        return instance;
    }

    /**
     * @brief Sets the minimum log level to be processed.
     * Messages with a lower severity will be ignored.
     * @param level The minimum log level.
     */
    void set_level (LogLevel level) {
        m_level = level;
    }

    /**
     * @brief Enables or disables logging to the console.
     * @param enable True to enable, false to disable.
     */
    void enable_console_logging (bool enable) {
        m_log_to_console = enable;
    }

    /**
     * @brief Sets the file path for log output.
     * If the path is set, logging to file is automatically enabled.
     * @param filepath Path to the log file.
     * @return True if the file was opened successfully, false otherwise.
     */
    bool set_logfile (const std::string& filepath) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_logfile.is_open()) {
            m_logfile.close();
        }
        m_logfile.open(filepath, std::ios::out | std::ios::app);
        m_log_to_file = m_logfile.is_open();
        if (!m_log_to_file) {
            std::cerr << "Error: Could not open log file: " << filepath << std::endl;
        }
        return m_log_to_file;
    }

    /**
     * @brief The core logging method.
     * @param level The severity level of the message.
     * @param message The log message content.
     * @param file The source file where the log was triggered.
     * @param line The line number in the source file.
     * @param function The function name where the log was triggered.
     */
    void log (LogLevel level, const std::string& message, const char* file, int line, const char* function) {
        if (level < m_level) {
            return;
        }

        std::lock_guard<std::mutex> lock(m_mutex);

        std::ostringstream formatted_message;
        formatted_message << get_timestamp() << " "
                          << level_to_string(level)
                          << " [" << get_filename(file) << ":" << line << " (" << function << ")] "
                          << message;

        if (m_log_to_console) {
            std::cout << formatted_message.str() << std::endl;
        }

        if (m_log_to_file && m_logfile.is_open()) {
            m_logfile << formatted_message.str() << std::endl;
        }
    }

private:
    /**
     * @brief Private constructor for singleton pattern.
     */
    Logger () : m_level(LogLevel::INFO), m_log_to_console(true), m_log_to_file(false) {}

    /**
     * @brief Private destructor. Closes the log file if open.
     */
    ~Logger () {
        if (m_logfile.is_open()) {
            m_logfile.close();
        }
    }

    /**
     * @brief Converts a LogLevel enum to its string representation.
     * @param level The log level to convert.
     * @return A string representing the log level, padded for alignment.
     */
    std::string level_to_string (LogLevel level) const {
        switch (level) {
            case LogLevel::DEBUG:    return "[DEBUG]  ";
            case LogLevel::INFO:     return "[INFO]   ";
            case LogLevel::WARNING:  return "[WARNING]";
            case LogLevel::ERROR:    return "[ERROR]  ";
            case LogLevel::CRITICAL: return "[CRITICAL]";
            default:                 return "[UNKNOWN]";
        }
    }

    /**
     * @brief Generates a current timestamp string.
     * @return A string with the current time in ISO 8601 format with milliseconds.
     */
    std::string get_timestamp () const {
        auto now = std::chrono::system_clock::now();
        auto now_time_t = std::chrono::system_clock::to_time_t(now);
        auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

        std::stringstream ss;
        // Use std::put_time for safe, formatted time output
        ss << std::put_time(std::localtime(&now_time_t), "%Y-%m-%d %H:%M:%S");
        ss << '.' << std::setfill('0') << std::setw(3) << now_ms.count();
        return ss.str();
    }

    /**
     * @brief Extracts the filename from a full path.
     * @param path The full path string (from __FILE__).
     * @return The basename of the file.
     */
    std::string get_filename (const std::string& path) const {
        size_t last_slash = path.find_last_of("/\\");
        return (last_slash == std::string::npos) ? path : path.substr(last_slash + 1);
    }

    // Member variables
    std::mutex m_mutex;
    LogLevel m_level;
    bool m_log_to_console;
    bool m_log_to_file;
    std::ofstream m_logfile;
};


/**
 * @class LogStream
 * @brief A helper class to enable stream-like logging syntax (e.g., LOG_INFO << "Value: " << 42;).
 *
 * This class captures the streamed data into a string buffer and, upon destruction
 * (when the statement ends), sends the complete message to the Logger.
 */
class LogStream {
public:
    LogStream (LogLevel level, const char* file, int line, const char* function)
        : m_level(level), m_file(file), m_line(line), m_function(function) {}

    // Destructor that logs the accumulated message
    ~LogStream () {
        Logger::get_instance().log(m_level, m_stream.str(), m_file, m_line, m_function);
    }

    // Overload the << operator to capture data of any type
    template <typename T>
    LogStream& operator<< (const T& msg) {
        m_stream << msg;
        return *this;
    }

private:
    std::ostringstream m_stream;
    LogLevel m_level;
    const char* m_file;
    int m_line;
    const char* m_function;
};

// --- User-facing Logging Macros ---
// These macros create a temporary LogStream object which captures the message
// and then logs it upon destruction at the end of the full expression.

#define LOG_DEBUG    LogStream(LogLevel::DEBUG,    __FILE__, __LINE__, LOG_FUNC_NAME)
#define LOG_INFO     LogStream(LogLevel::INFO,     __FILE__, __LINE__, LOG_FUNC_NAME)
#define LOG_WARNING  LogStream(LogLevel::WARNING,  __FILE__, __LINE__, LOG_FUNC_NAME)
#define LOG_ERROR    LogStream(LogLevel::ERROR,    __FILE__, __LINE__, LOG_FUNC_NAME)
#define LOG_CRITICAL LogStream(LogLevel::CRITICAL, __FILE__, __LINE__, LOG_FUNC_NAME)


/*
// =================================================================================
//                                 EXAMPLE USAGE
// =================================================================================
//
// To use the logger, simply include this header file ("logging.h").
// No .cpp file is needed.
//
// int main () {
//     // --- Configuration (Optional) ---
//
//     // Get the logger instance
//     auto& logger = Logger::get_instance();
//
//     // Set a log file. Logging to this file will be enabled automatically.
//     if (!logger.set_logfile("simulation.log")) {
//         return 1; // Exit if file can't be opened
//     }
//
//     // Set the minimum level of messages to log.
//     // For example, to see everything including debug messages:
//     logger.set_level(LogLevel::DEBUG);
//
//     // To disable console output and only log to the file:
//     // logger.enable_console_logging(false);
//
//
//     // --- Logging Examples ---
//
//     LOG_INFO << "Starting soft particle simulation.";
//
//     double timestep = 0.005;
//     int total_steps = 10000;
//     LOG_INFO << "Simulation parameters: dt = " << timestep << ", steps = " << total_steps;
//
//     LOG_DEBUG << "This is a detailed debug message for developers.";
//
//     for (int i = 0; i <= total_steps; ++i) {
//         if (i % 1000 == 0) {
//             LOG_INFO << "Progress: Step " << i << " of " << total_steps;
//         }
//
//         // Simulate a potential issue
//         if (i == 5000) {
//             LOG_WARNING << "Particle 123 has unusually high velocity.";
//         }
//     }
//
//     // Simulate an error condition
//     double energy = -1.0; // Invalid energy
//     if (energy < 0) {
//         LOG_ERROR << "System energy is negative (" << energy << "), which is non-physical.";
//     }
//
//     // Simulate a critical failure
//     LOG_CRITICAL << "Simulation diverged. Halting execution.";
//
//
//     // The logger will automatically clean up when the program exits.
//     return 0;
// }
//
*/

#endif // LOGGER_H
