#include "include/system.h"
#include "utils/include/hdf5_logger.h"
#include "utils/include/logger.h"
#include <iomanip>   // For std::setprecision, std::fixed
#include <limits>    // For std::numeric_limits

std::string format_string (size_t i, size_t sz) {
    std::ostringstream oss;
    oss << std::setw(sz) << std::setfill('0') << i;
    return oss.str();
}

int main (int argc, char** argv) {
    HDF5Logger file_logger("run.h5");
    file_logger.create_group("/frames");
    System sys(PARTICLE_TYPE::BUMPY, 30, 8, 0.0, 1e-3, 1/sqrt(3.0));
    LOG_INFO << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << sys.particles[0]->area;

    LOG_INFO << "starting";
    // sys.send_to_jamming();
    sys.grow_to_phi(0.50);
    sys.set_energy(1e-3);

    LOG_INFO << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << sys.L;
    // LOG_INFO << sys.active_p;

    for (size_t i = 0; i < 1e6; i++) {
        if (i % 10000 == 0) {
            std::ostringstream oss;
            oss << std::setw(7) << std::setfill('0') << i / 10000;
            std::string name = "/frames/" + oss.str();
            auto particles = sys.get_verts();
            file_logger.write_dataset(name, {sys.active_p, 8, 2}, particles);
            file_logger.write_attribute(name, "sigma", sys.get_sigma());
        }
        sys.update();
    }

    LOG_INFO << "ended";

    return 0;
}
