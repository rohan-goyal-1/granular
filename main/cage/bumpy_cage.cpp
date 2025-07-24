#include "include/system.h"
#include "utils/include/logger.h"
#include "utils/include/hdf5_logger.h"

#include "protocol/include/cage_sampling.h"

#include <omp.h>

std::string format_string (size_t i, size_t sz) {
    std::ostringstream oss;
    oss << std::setw(sz) << std::setfill('0') << i;
    return oss.str();
}

int main (int argc, char** argv) {
    HDF5Logger logger("bumpy_out.h5");

    std::vector<size_t> N_s = {5, 10, 15};

    for (size_t x = 0; x < 3; x++) {
        size_t N = N_s[x];

        std::string n_name  = "/" + format_string(x, 7);
        #pragma omp critical
        {
            logger.create_group(n_name);
            logger.write_attribute(n_name, "N", N);
        }

        #pragma omp parallel for
        for (int i = 0; i < 1e2; i++) {
            std::string name = n_name + "/" + format_string(i, 7);

            System sys(PARTICLE_TYPE::BUMPY, 5, N, 1.0, 0.0, 0, 1e-3);
            sys.send_to_jamming();

            #pragma omp critical
            LOG_INFO << "Finished jamming: N=" << N << ";\ttrial: " << i;

            double phi_j = sys.phi;

            std::vector<Eigen::VectorXd> particles;
            for (auto& p : sys.particles)
                particles.push_back(p->com);

            auto [cage, neighbors] = generate_cage(particles, 0);
            if (cage.empty()) continue;
            #pragma omp critical
            LOG_INFO << "Finished generating Voronoi: N=" << N << ";\ttrial: " << i;

            size_t num_samples = 1'000'000;
            auto sample_pts = generate_sample_points(cage, num_samples);

            // Write trial group
            #pragma omp critical
            logger.create_group(name);

            int step = 0;
            for (double d_phi = 0; d_phi <= 0.50; d_phi += 1e-3, step++) {
                std::string trial_name = name + "/" + format_string(step, 7);

                sys.rescale(phi_j - d_phi);
                double f = evaluate(sys, sample_pts);

                #pragma omp critical
                {
                    logger.create_group(trial_name);
                    logger.write_attribute(trial_name, "fraction", f);
                    logger.write_attribute(trial_name, "d_phi", d_phi);
                }
            }
            #pragma omp critical
            LOG_INFO << "Finished trial: N=" << N << ";\ttrial: " << i;
        }
        #pragma omp critical
        LOG_INFO << "Finished N: N=" << N;
    }

    return 0;
}
