#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <omp.h>

#include "include/system.h"
#include "protocol/include/diffusion.h"
#include "utils/include/file_handler.h"
#include "utils/include/logger.h"

int main (int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: <program_name> <output_file> <run_time>\n"
                  << "The program expects an output file to store the result and a run time for which to check diffusion." << std::endl;
        return EXIT_FAILURE;
    }

    auto& out_logger = Logger::get_instance();
    out_logger.set_level(LogLevel::INFO);

    const std::string out_path = argv[1];
    const size_t run_time = std::stoi(argv[2]);

    std::vector<double> phis = {0.10, 0.20, 0.30, 0.40, 0.50, 0.60};
    PlotLogger logger(out_path);
    logger.add_plot("diffusion", "Diffusion of 3 particles", "log(Δt)", "log(r²)");

    System system(SYSTEM_TYPE::SOFT_MONO, 3, 3, 1.0, 0.0, 1e-2, 1e-3);
    LOG_INFO << "Starting diffusion runs...";
    auto data = find_diffusion(&system, 0.1, 1e3, run_time, 30, phis);
    for (size_t i = 0; i < phis.size(); i++) {
        logger.add_dataset("diffusion", std::to_string(phis[i]), data[i], "φ=" + std::to_string(phis[i]));
    }

    return EXIT_SUCCESS;
}
