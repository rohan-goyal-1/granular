#include "include/system.h"
#include "utils/include/logger.h"
#include "utils/include/hdf5_logger.h"

#include "protocol/include/cage_sampling.h"

int main (int argc, char** argv) {
    HDF5Logger logger("sample3.h5");

    System sys(PARTICLE_TYPE::BUMPY, 3, 5, 1.0, 0.0, 0, 1e-3);
    sys.send_to_jamming();
    LOG_INFO << "Generated jammed configuration";
    sys.remove_rattlers();

    // for (double d_phi = 0.001; sys.get_pe() != 0.05; ) {
    //     sys.rescale(sys.phi - d_phi);
    // }

    std::unordered_map<size_t, size_t> inds;
    std::vector<Eigen::VectorXd> particles;
    for (size_t k = 0; k < 3; k++) if (sys.active[k]) {
        inds[particles.size()] = k;
        particles.push_back(sys.particles[k]->com);
    }

    // Log jammed configuration
    logger.write_dataset("/particles", {particles.size(), 5, 2}, sys.get_verts());
    logger.write_attribute("/particles", "sigma", sys.get_sigma());

    sys.rescale(sys.phi - 0.2);

    // Log after decompression
    logger.write_dataset("/n_particles", {particles.size(), 5, 2}, sys.get_verts());
    logger.write_attribute("/n_particles", "sigma", sys.get_sigma());

    auto [cage, neighbors] = generate_cage(particles, 0);
    if (cage.empty()) return 1;

    // Log cage
    std::vector<double> flat_cage;
    for (const auto& v : cage) {
        flat_cage.push_back(v.x());
        flat_cage.push_back(v.y());
    }
    logger.write_dataset("/cage", {cage.size(), 2}, flat_cage);

    size_t num_samples = 100'000;
    auto sample_pts = generate_sample_points(cage, num_samples);

    std::vector<double> data;
    data.reserve(num_samples * 3);
    size_t cnt = 0;
    double area = 0.0;
    for (const auto& pt : sample_pts) {
        sys.particles[inds[0]]->move(pt - sys.particles[inds[0]]->com);
        // bool works = true;
        // for (const size_t neigh : neighbors) if (neigh != -1) {
        //     works &= (sys.particles[inds[neigh]]->get_energy_interaction(sys.particles[inds[0]].get()) == 0.0);
        // }
        size_t angles_worked = 0;
        size_t angles_sampled = 1'000;
        double d_theta = 2 * M_PI / 5 / angles_sampled;
        for (size_t j = 0; j < angles_sampled; j++) {
            sys.particles[inds[0]]->rotate(d_theta);
            angles_worked += (sys.get_pe() == 0.0);
        }
        double a = (double) angles_worked / angles_sampled;
        data.push_back(pt.x());
        data.push_back(pt.y());
        data.push_back(a);
        LOG_INFO << "Finished point: " << cnt++;
        area += a;
    }
    area /= num_samples;
    LOG_INFO << "A=" << area;

    logger.write_dataset("/samples", {sample_pts.size(), 3}, data);

    return 0;
}
