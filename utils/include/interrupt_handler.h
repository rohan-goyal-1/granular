#ifndef INTERRUPT_HANDLER_H
#define INTERRUPT_HANDLER_H

#include <iostream>
#include <csignal>
#include <atomic>
#include <functional>

/**
 * @class InterruptHandler
 * @brief A singleton to handle process signals (like SIGINT from Ctrl+C) for graceful shutdown.
 *
 * This class sets up a signal handler that safely sets an atomic flag upon
 * receiving an interrupt signal. Long-running loops in the application can
 * poll this flag to know when to terminate gracefully.
 */
class InterruptHandler {
public:
    /**
     * @brief Get the singleton instance of the InterruptHandler.
     * The first call to this function will set up the signal handlers.
     * @return Reference to the single InterruptHandler instance.
     */
    static InterruptHandler& get_instance () {
        static InterruptHandler instance;
        return instance;
    }

    /**
     * @brief Checks if an interrupt signal has been received.
     * This function is thread-safe and designed to be called frequently
     * from within long-running loops.
     * @return True if an interrupt has been caught, false otherwise.
     */
    static bool is_interrupted () {
        return s_interrupt_flag.load();
    }

    // Delete copy/move constructors and assignment operators for singleton
    InterruptHandler(const InterruptHandler&) = delete;
    InterruptHandler& operator=(const InterruptHandler&) = delete;
    InterruptHandler(InterruptHandler&&) = delete;
    InterruptHandler& operator=(InterruptHandler&&) = delete;

private:
    /**
     * @brief Private constructor to set up the signal handling.
     * This is called only once when the singleton instance is created.
     */
    InterruptHandler () {
        // Register our static signal_handler function to be called for SIGINT and SIGTERM
        std::signal(SIGINT, InterruptHandler::signal_handler);
        std::signal(SIGTERM, InterruptHandler::signal_handler);
    }

    /**
     * @brief The actual signal handler function.
     * IMPORTANT: This function must be static. It has C-style linkage to be
     * compatible with std::signal. It should do the minimum work possible,
     * such as setting a flag. Avoid any complex operations, memory allocation,
     * or I/O within the handler.
     * @param signal The signal number received.
     */
    static void signal_handler (int signal) {
        // We received an interrupt signal. Set the atomic flag to true.
        // The main application loop will pick this up.
        s_interrupt_flag.store(true);
    }

    // The static atomic flag that is shared across all parts of the application.
    // 'static' ensures there's only one instance of this flag.
    // 'atomic' ensures that reads/writes are safe from data races in a multi-threaded context.
    static inline std::atomic<bool> s_interrupt_flag{false};
};


/*
// =================================================================================
//                                 EXAMPLE USAGE
// =================================================================================
//
// This shows how to integrate InterruptHandler with the Logger and a typical
// multi-threaded simulation structure.
//
// #include "logging.h"
// #include "interrupt_handler.h"
// #include <thread>
// #include <vector>
//
// // A function that represents the work a thread does in a simulation
// void simulation_worker(int thread_id, int steps) {
//     LOG_INFO << "Worker thread " << thread_id << " starting.";
//
//     for (int i = 0; i < steps; ++i) {
//         // --- CHECK FOR INTERRUPT ---
//         // The most important part: check if a shutdown has been requested.
//         if (InterruptHandler::is_interrupted()) {
//             LOG_WARNING << "Worker " << thread_id << ": Interrupt detected at step " << i << ". Stopping early.";
//             break; // Exit the loop gracefully
//         }
//
//         // Simulate some work
//         std::this_thread::sleep_for(std::chrono::milliseconds(2));
//
//         if (i % 1000 == 0) {
//             LOG_DEBUG << "Worker " << thread_id << " progress: step " << i;
//         }
//     }
//
//     LOG_INFO << "Worker thread " << thread_id << " finished its work.";
// }
//
// void perform_final_save() {
//     LOG_INFO << "Performing final data save before exiting...";
//     // Simulate saving data to a file
//     std::this_thread::sleep_for(std::chrono::milliseconds(500));
//     LOG_INFO << "Final data save complete.";
// }
//
// int main() {
//     // --- INITIALIZATION ---
//
//     // 1. Set up the logger
//     auto& logger = Logger::get_instance();
//     logger.set_logfile("interrupt_sim.log");
//     logger.set_level(LogLevel::DEBUG);
//
//     // 2. Initialize the interrupt handler (just by calling get_instance)
//     InterruptHandler::get_instance();
//
//     LOG_INFO << "Main: Simulation starting. Press Ctrl+C to interrupt gracefully.";
//
//     // --- SIMULATION SETUP ---
//     const int num_threads = 4;
//     const int total_steps = 5000; // ~10 seconds of work
//     std::vector<std::thread> threads;
//
//     for (int i = 0; i < num_threads; ++i) {
//         threads.emplace_back(simulation_worker, i, total_steps);
//     }
//
//     // --- MAIN THREAD WAITS ---
//     // The main thread can do other work or just wait for threads to complete.
//     for (auto& t : threads) {
//         if (t.joinable()) {
//             t.join();
//         }
//     }
//
//     // --- GRACEFUL SHUTDOWN ---
//
//     // This part of the code is reached either when the simulation finishes normally
//     // or when it was interrupted and the worker threads exited their loops.
//
//     if (InterruptHandler::is_interrupted()) {
//         LOG_CRITICAL << "Main: Shutdown initiated by user (Ctrl+C).";
//     }
//     else {
//         LOG_INFO << "Main: Simulation completed normally.";
//     }
//
//     // Perform final cleanup, like saving aggregate data.
//     // This is guaranteed to happen before the program exits.
//     perform_final_save();
//
//     LOG_INFO << "Main: All tasks are complete. Exiting.";
//
//     return 0; // The logger's destructor will be called here, ensuring the log file is closed.
// }
//
*/

#endif // INTERRUPT_HANDLER_H
