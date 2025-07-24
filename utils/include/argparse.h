#ifndef ARGPARSE_H
#define ARGPARSE_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <algorithm>

/**
 * @class ArgParser
 * @brief A simple, header-only C++ command line argument parser inspired by Python's argparse.
 *
 * Provides a fluent interface to define expected arguments, parse them from
 * argc/argv, and retrieve their values in a type-safe manner.
 */
class ArgParser {
private:
    // Abstract base class for type erasure of argument details
    struct ArgumentBase {
        virtual ~ArgumentBase() = default;
        virtual void parse(const std::string& value) = 0;
        virtual std::string get_type_name() const = 0;
        std::string name;
        std::string short_name;
        std::string help_text;
        bool is_required = false;
        bool is_flag = false;
        bool has_been_set = false;
        std::string default_value_str;
    };

    // Templated derived class to hold argument value of a specific type
    template <typename T>
    struct Argument : public ArgumentBase {
        T value;

        void parse (const std::string& value_str) override {
            std::istringstream iss(value_str);
            if (!(iss >> value)) {
                throw std::invalid_argument("Invalid value for argument '" + name + "': " + value_str);
            }
        }

        std::string get_type_name () const override {
            if constexpr (std::is_same_v<T, int>) return "integer";
            if constexpr (std::is_same_v<T, double>) return "float";
            if constexpr (std::is_same_v<T, float>) return "float";
            if constexpr (std::is_same_v<T, std::string>) return "string";
            if constexpr (std::is_same_v<T, bool>) return "flag";
            return "unknown";
        }
    };

    // Specialization for bool (flags)
    template <>
    struct Argument<bool> : public ArgumentBase {
        bool value = false;

        Argument () {
            is_flag = true;
        }

        void parse (const std::string& value_str) override {
            // Flags are set to true if present, no value needed.
            value = true;
        }

        std::string get_type_name() const override {
            return "flag";
        }
    };


public:
    /**
     * @brief Constructs an argument parser.
     * @param app_name The name of the application, used in the help message.
     */
    ArgParser (const std::string& app_name = "") : m_app_name(app_name) {
        add_argument<bool>("--help", "-h").help("Show this help message and exit");
    }

    /**
     * @brief Adds a new argument to be parsed.
     * @tparam T The type of the argument's value (e.g., int, double, string).
     * @param name The long name of the argument (e.g., "--file").
     * @param short_name The optional short name (e.g., "-f").
     * @return A reference to the argument being configured, allowing for chained calls.
     */
    template <typename T>
    ArgParser& add_argument (const std::string& name, const std::string& short_name = "") {
        auto arg = std::make_unique<Argument<T>>();
        arg->name = name;
        arg->short_name = short_name;

        m_arg_map[name] = arg.get();
        if (!short_name.empty()) {
            m_arg_map[short_name] = arg.get();
        }

        m_args_in_order.push_back(std::move(arg));
        return *this;
    }

    /**
     * @brief Sets the help text for the most recently added argument.
     */
    ArgParser& help (const std::string& help_text) {
        if (m_args_in_order.empty()) throw std::logic_error("Cannot set help for a non-existent argument.");
        m_args_in_order.back()->help_text = help_text;
        return *this;
    }

    /**
     * @brief Marks the most recently added argument as required.
     */
    ArgParser& required () {
        if (m_args_in_order.empty()) throw std::logic_error("Cannot set required for a non-existent argument.");
        m_args_in_order.back()->is_required = true;
        return *this;
    }

    /**
     * @brief Sets the default value for the most recently added argument.
     */
    template<typename T>
    ArgParser& default_value (const T& val) {
        if (m_args_in_order.empty()) throw std::logic_error("Cannot set default value for a non-existent argument.");
        std::ostringstream oss;
        oss << val;
        m_args_in_order.back()->default_value_str = oss.str();
        return *this;
    }


    /**
     * @brief Parses the command line arguments.
     * @param argc The argument count from main().
     * @param argv The argument vector from main().
     */
    void parse (int argc, char** argv) {
        m_app_name = (m_app_name.empty() && argc > 0) ? argv[0] : m_app_name;
        std::vector<std::string> tokens(argv + 1, argv + argc);

        for (size_t i = 0; i < tokens.size(); ++i) {
            const std::string& token = tokens[i];

            if (m_arg_map.find(token) == m_arg_map.end()) {
                throw std::runtime_error("Unknown argument: " + token);
            }

            ArgumentBase* arg = m_arg_map[token];
            arg->has_been_set = true;

            if(arg->name == "--help" || arg->short_name == "--help") {
                print_help();
                exit(0);
            }

            if (arg->is_flag) {
                arg->parse("true");
            }
            else {
                if (i + 1 >= tokens.size()) {
                    throw std::runtime_error("Argument requires a value: " + token);
                }
                arg->parse(tokens[++i]);
            }
        }

        // Apply defaults and check for required arguments
        for (const auto& arg_ptr : m_args_in_order) {
            if (!arg_ptr->has_been_set) {
                if (arg_ptr->is_required) {
                    throw std::runtime_error("Required argument missing: " + arg_ptr->name);
                }
                if (!arg_ptr->default_value_str.empty()) {
                    arg_ptr->parse(arg_ptr->default_value_str);
                }
            }
        }
    }

    /**
     * @brief Retrieves the value of a parsed argument.
     * @tparam T The expected type of the argument.
     * @param name The name of the argument (e.g., "--file").
     * @return The value of the argument.
     */
    template <typename T>
    T get (const std::string& name) const {
        if (m_arg_map.find(name) == m_arg_map.end()) {
            throw std::runtime_error("Request for unknown argument: " + name);
        }

        auto* base_ptr = m_arg_map.at(name);
        auto* derived_ptr = dynamic_cast<Argument<T>*>(base_ptr);

        if (!derived_ptr) {
            throw std::runtime_error("Type mismatch for argument: " + name);
        }
        return derived_ptr->value;
    }

private:
    void print_help () const {
        std::cout << "Usage: " << m_app_name << " [options]\n\n";
        std::cout << "Options:\n";

        for (const auto& arg_ptr : m_args_in_order) {
            std::string names = "  " + arg_ptr->name;
            if (!arg_ptr->short_name.empty()) {
                names += ", " + arg_ptr->short_name;
            }

            std::cout << "    " << std::left << std::setw(24) << names;

            std::string help = arg_ptr->help_text;
            if(!arg_ptr->default_value_str.empty()){
                help += " (default: " + arg_ptr->default_value_str + ")";
            }
            std::cout << help << "\n";
        }
    }

    std::string m_app_name;
    std::vector<std::unique_ptr<ArgumentBase>> m_args_in_order;
    std::map<std::string, ArgumentBase*> m_arg_map;
};


/*
// =================================================================================
//                                 EXAMPLE USAGE
// =================================================================================
//
// #include "argparse.h"
// #include "logging.h" // Assuming you have the logger from previous example
//
// int main(int argc, char** argv) {
//     // 1. Define the parser and the arguments it should expect.
//     ArgParser parser("particle_simulator");
//
//     parser.add_argument<std::string>("--input", "-i")
//           .help("Input file for particle starting positions.")
//           .required();
//
//     parser.add_argument<std::string>("--output", "-o")
//           .help("Output file for final particle positions.")
//           .default_value("output.dat");
//
//     parser.add_argument<int>("--steps", "-s")
//           .help("Number of simulation steps to run.")
//           .default_value(10000);
//
//     parser.add_argument<double>("--dt")
//           .help("Timestep for the simulation.")
//           .default_value(0.001);
//
//     parser.add_argument<bool>("--verbose", "-v")
//           .help("Enable verbose logging output.");
//
//     // 2. Parse the arguments.
//     // This will automatically handle --help and exit if it's present.
//     // It will throw an exception on error (e.g., missing required arg).
//     try {
//         parser.parse(argc, argv);
//     } catch (const std::exception& e) {
//         std::cerr << "Error: " << e.what() << std::endl;
//         std::cerr << "Use --help for more information." << std::endl;
//         return 1;
//     }
//
//     // 3. Retrieve the parsed values.
//     std::string input_file = parser.get<std::string>("--input");
//     std::string output_file = parser.get<std::string>("--output");
//     int steps = parser.get<int>("--steps");
//     double dt = parser.get<double>("--dt");
//     bool verbose = parser.get<bool>("--verbose");
//
//     // 4. Use the values in your application.
//     // Example of integrating with the logger.
//     auto& logger = Logger::get_instance();
//     if (verbose) {
//         logger.set_level(LogLevel::DEBUG);
//     }
//     else {
//         logger.set_level(LogLevel::INFO);
//     }
//
//     LOG_INFO << "Simulation configured to run.";
//     LOG_INFO << "  Input file: " << input_file;
//     LOG_INFO << "  Output file: " << output_file;
//     LOG_INFO << "  Steps: " << steps;
//     LOG_INFO << "  Timestep: " << dt;
//     LOG_DEBUG << "Verbose logging is enabled.";
//
//     // ... Your simulation code would go here ...
//
//     LOG_INFO << "Simulation finished.";
//
//     return 0;
// }
//
// --- Example Command Line Invocations ---
//
// ./particle_simulator --input initial.dat
// ./particle_simulator --input initial.dat --output final.dat -s 50000
// ./particle_simulator --dt 0.0005 -i particles.in -v
// ./particle_simulator --help
//
*/

#endif // ARGPARSE_H
