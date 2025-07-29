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
    System sys(PARTICLE_TYPE::BUMPY, 1, 5, 0.0, 1e-3, 0.5773502692);
    // auto& logger = Logger::get_instance();
    // logger.set_level(LogLevel::INFO);

    LOG_INFO << "area " << sys.get_particle_area();
    LOG_INFO << "sigma " << sys.get_sigma();
    // LOG_INFO << "starting jamming";
    // sys.send_to_jamming();
    // LOG_INFO << "ended";
    // sys.purge_rattlers();
    // // for (int i = 0; i < 3; i++) {
    // //     auto* bumpy = dynamic_cast<BumpyParticle*>(sys.particles[i].get());
    // //     for (int j = 0; j < 3; j++)
    // //         std::cout << bumpy->verts[j] << '\n';
    // //     std::cout << '\n';
    // // }
    // // std::cout << sys.L << '\n';
    // // for (size_t i = 0; i < 5; i++) {
    // //     auto* p = dynamic_cast<BumpyParticle*>(sys.particles[i].get());
    // //     std::cout << p->moi << ' ' << p->com.x() << ' ' << p->com.y() << ' ' << p->sigma << ' ' << p->theta << '\n';
    // // }
    //
    // // sys.grow_to_phi(0.3);
    // // sys.set_ke(1e-5);
    // // LOG_INFO << sys.get_sigma();
    // // // LOG_INFO << sys.phi;
    // LOG_INFO << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << sys.L;
    // //
    // // LOG_INFO << "pe " << sys.get_pe();
    // for (size_t i = 0; i < 1e0; i++) {
    //     if (i % 10000 == 0) {
    //         std::ostringstream oss;
    //         oss << std::setw(7) << std::setfill('0') << i / 1000;
    //         std::string name = "/frames/" + oss.str();
    //         auto particles = sys.get_verts();
    //         file_logger.write_dataset(name, {sys.active_p, 8, 2}, particles);
    //         file_logger.write_attribute(name, "sigma", sys.get_sigma());
    //     }
    //     // sys.update();
    // }
    //
    // LOG_INFO << "ended";

    return 0;
}
