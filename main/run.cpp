#include "include/system.h"
#include "utils/include/hdf5_logger.h"
#include "utils/include/logger.h"

std::string format_string (size_t i, size_t sz) {
    std::ostringstream oss;
    oss << std::setw(sz) << std::setfill('0') << i;
    return oss.str();
}

int main (int argc, char** argv) {
    HDF5Logger file_logger("run.h5");
    file_logger.create_group("/frames");
    System sys(PARTICLE_TYPE::BUMPY, 5, 8, 0.0, 1e-3, 1/sqrt(15.0));
    // auto& logger = Logger::get_instance();
    // logger.set_level(LogLevel::INFO);

    sys.send_to_jamming();
    std::cout << sys.N_c << '\n';
    // for (size_t i = 0; i < 5; i++) {
    //     auto* p = dynamic_cast<BumpyParticle*>(sys.particles[i].get());
    //     std::cout << p->moi << ' ' << p->com.x() << ' ' << p->com.y() << ' ' << p->sigma << ' ' << p->theta << '\n';
    // }

    // LOG_INFO << system.phi;
    auto particles = sys.get_verts();
    file_logger.write_dataset("/frames/0", {5, 8, 2}, particles);
    file_logger.write_attribute("/frames/0", "sigma", sys.get_sigma());

    // for (size_t i = 0; i < 1e6; i++) {
    //     if (i % 1000 == 0) {
    //         std::ostringstream oss;
    //         oss << std::setw(7) << std::setfill('0') << i / 1000;
    //         std::string name = "/frames/" + oss.str();
    //     }
    //     system.update();
    // }

    return 0;
}
