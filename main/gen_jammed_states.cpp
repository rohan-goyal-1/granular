// TODO: come up with shutdown protocol and run time error protocol

#include <iostream>
#include <iomanip>
#include <string>

#include <omp.h>

#include "include/system.h"
#include "utils/include/file_handler.h"

int main (int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: <program_name> <input_file> <output_file> <run_times>\n"
                  << "The program expects an input file to parse, an output file to store the result, and an iteration count." << std::endl;
        return EXIT_FAILURE;
    }

    const std::string in_path = argv[1];
    const std::string out_path = argv[2];
    const size_t ITERS = std::stoi(argv[3]);

    Metadata md = load_from_file(in_path);
    SystemLogger logger(out_path);

    #pragma omp parallel for
    for (size_t i = 0; i < ITERS; i++) {
        System system(md);
        system.send_to_jamming();

        #pragma omp critical
        {
            logger.log_system(system);
        }
    }

    return EXIT_SUCCESS;
}
