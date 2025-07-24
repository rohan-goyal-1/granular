#include "include/system.h"
#include "utils/include/hdf5_logger.h"
#include "utils/include/logger.h"
#include "utils/include/argparse.h"
#include "protocol/include/cage_sampling.h"

#include <unordered_map>
#include <algorithm>

std::string format_string (size_t i, size_t sz) {
    std::ostringstream oss;
    oss << std::setw(sz) << std::setfill('0') << i;
    return oss.str();
}

int main (int argc, char** argv) {
    ArgParser parser("iter_decomp");
    parser.add_argument<bool>("--bumpy", "-b").help("Use bumpy particles");
    parser.add_argument<size_t>("--nv", "-v").help("Number of vertices for bumpy particles.").default_value(1);
    parser.add_argument<std::string>("--output", "-o").help("The path to the file to write results to.").default_value("out.h5");
    parser.add_argument<size_t>("--np", "-n").help("The number of particles to use.");
    parser.add_argument<size_t>("--points", "-p").help("Number of points to test.").default_value(1e5);
    parser.add_argument<double>("--dt", "-d").help("The timestep for the simulation.").default_value(1e-3);
    parser.add_argument<double>("--dphi").help("Amount of final decompression from jamming.");
    parser.add_argument<size_t>("--steps").help("Number of log-scale decompression steps.");
    parser.add_argument<size_t>("--angles", "-a").help("Number of angle samples to try for bumpy particles.").default_value(1e3);
    parser.add_argument<double>("--mu", "-m").help("Effective friction coefficient for the particles.");

    try {
        parser.parse(argc, argv);
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Use --help for more information." << std::endl;
        return 1;
    }

    PARTICLE_TYPE p_type = parser.get<bool>("--bumpy") ? PARTICLE_TYPE::BUMPY : PARTICLE_TYPE::SMOOTH;
    std::string output_file = parser.get<std::string>("--output");
    size_t Nv = parser.get<size_t>("--nv");
    size_t N = parser.get<size_t>("--np");
    size_t points = parser.get<size_t>("--points");
    double dt = parser.get<double>("--dt");
    double d_phi = parser.get<double>("--dphi");
    size_t num_steps = parser.get<size_t>("--steps");
    size_t angles_sampled = parser.get<size_t>("--angles");
    double mu_eff = parser.get<double>("--mu");

    auto& std_logger = Logger::get_instance();
    std_logger.set_level(LogLevel::INFO);

    HDF5Logger logger(output_file);
    System sys(p_type, N, Nv, 0.0, dt, mu);
    sys.send_to_jamming();
    sys.remove_rattlers();
    LOG_INFO << "Jammed system and removed rattlers.";
    double phi_j = sys.phi;

    std::vector<Eigen::Vector2d> active_particle_COMs;
    std::unordered_map<size_t, size_t> global_to_local_map;
    std::unordered_map<size_t, size_t> local_to_global_map;
    size_t num_active = 0;
    for (size_t i = 0; i < N; ++i) {
        if (sys.active[i]) {
            global_to_local_map[i] = num_active;
            local_to_global_map[num_active] = i;
            active_particle_COMs.push_back(sys.particles[i]->com);
            num_active++;
        }
    }

    logger.write_attribute("/", "Np", int(num_active));
    logger.write_attribute("/", "Nv", int(Nv));
    logger.write_attribute("/", "bumpy", (int) (p_type == PARTICLE_TYPE::BUMPY));

    std::vector<double> flat_dynam_matrix;
    flat_dynam_matrix.reserve(sys.adj_contacts.size());
    for (int i = 0; i < N * Nv; i++) {
        for (int j = 0; j < N * Nv; j++) {
            flat_dynam_matrix.push_back(sys.adj_contacts(i, j));
        }
    }
    logger.write_dataset("/adj", {N * Nv, N * Nv}, flat_dynam_matrix);

    std::vector<double> vertex_data;
    for (size_t i = 0; i < N; ++i) if (sys.active[i]) {
        if (p_type == PARTICLE_TYPE::BUMPY) {
            auto* p = dynamic_cast<BumpyParticle*>(sys.particles[i].get());
            vertex_data.push_back(p->com.x());
            vertex_data.push_back(p->com.y());
            vertex_data.push_back(p->theta);
        }
        else {
            auto* p = dynamic_cast<SmoothParticle*>(sys.particles[i].get());
            vertex_data.push_back(p->com.x());
            vertex_data.push_back(p->com.y());
            vertex_data.push_back(0);
        }
    }
    logger.write_dataset("/vertices", {num_active, 3}, vertex_data);
    logger.write_attribute("/vertices", "sigma", sys.get_sigma());

    logger.create_group("/cages");
    std::vector<std::vector<Eigen::Vector2d>> all_cages(num_active);
    std::vector<std::vector<Eigen::Vector2d>> all_sample_points(num_active);

    LOG_INFO << "Generating cages and sample points for all active particles...";
    for (size_t local_i = 0; local_i < num_active; ++local_i) {
        all_cages[local_i] = generate_minkowski_cage(
            local_i,
            active_particle_COMs,
            local_to_global_map,
            sys.particles
        );

        if (!all_cages[local_i].empty()) {
            std::vector<double> cage_data;
            for (const auto& pt : all_cages[local_i]) {
                cage_data.push_back(pt.x());
                cage_data.push_back(pt.y());
            }
            logger.write_dataset("/cages/" + format_string(local_i, 6),
                                 {all_cages[local_i].size(), 2}, cage_data);

            all_sample_points[local_i] = generate_sample_points(all_cages[local_i], points);
        }
    }
    LOG_INFO << "Done generating cages.";

    auto find_relevant_neighbors = [&] (size_t center_local_idx, const std::vector<Eigen::Vector2d>& cage) {
        std::vector<size_t> neighbors;
        size_t center_global_idx = local_to_global_map.at(center_local_idx);
        const auto& center_particle = sys.particles[center_global_idx];
        double R_center = center_particle->radius + 0.5 * center_particle->sigma;
        double L = sys.L;

        auto point_to_segment_dist = [] (const Eigen::Vector2d& p, const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
            Eigen::Vector2d ab = b - a;
            Eigen::Vector2d ap = p - a;
            double t = ab.dot(ap) / ab.squaredNorm();
            t = std::clamp(t, 0.0, 1.0);
            return (a + t * ab - p).norm();
        };

        for (size_t j = 0; j < num_active; ++j) {
            if (j == center_local_idx) continue;
            size_t j_global = local_to_global_map.at(j);
            const auto& pj = sys.particles[j_global];
            double Rj = pj->radius + 0.5 * pj->sigma;
            double effective_radius = R_center + Rj;

            // Find minimum image COM of pj relative to the cage center
            Eigen::Vector2d cage_center(0, 0);
            for (const auto& pt : cage) cage_center += pt;
            cage_center /= cage.size();

            Eigen::Vector2d shifted_com = minimum_image(pj->com - cage_center, L) + cage_center;

            for (size_t i = 0; i < cage.size(); ++i) {
                const auto& a = cage[i];
                const auto& b = cage[(i + 1) % cage.size()];
                if (point_to_segment_dist(shifted_com, a, b) < effective_radius) {
                    neighbors.push_back(j_global);
                    break;
                }
            }
        }

        return neighbors;
    };

    std::vector<double> delta_phis = {0.0};
    if (num_steps > 1) {
        double log_min = std::min(std::log10(d_phi), -6.0);
        double log_max = std::log10(d_phi);
        double log_step = (log_max - log_min + 1) / (num_steps - 1);
        for (size_t i = 0; i < num_steps - 1; ++i) {
            delta_phis.push_back(std::pow(10.0, log_min + i * log_step));
        }
    }

    for (size_t decompress_idx = 0; decompress_idx < delta_phis.size(); ++decompress_idx) {
        double current_dphi = delta_phis[decompress_idx];
        double phi_target = phi_j - current_dphi;
        sys.rescale(phi_target);

        LOG_INFO << "Processing decompression step " << decompress_idx << " (dphi = " << current_dphi << ")";

        std::string decompress_group_path = "/" + format_string(decompress_idx, 6);
        logger.create_group(decompress_group_path);
        logger.write_attribute(decompress_group_path, "L", sys.L);
        logger.write_attribute(decompress_group_path, "dphi", current_dphi);

        for (size_t local_i = 0; local_i < num_active; ++local_i) {
            size_t global_i = local_to_global_map[local_i];
            const auto& sample_pts = all_sample_points[local_i];
            if (sample_pts.empty()) continue;

            std::vector<size_t> neighbor_globals = find_relevant_neighbors(local_i, all_cages[local_i]);
            std::cout << neighbor_globals.size() << '\n';
            Eigen::Vector2d original_pos = sys.particles[global_i]->com;
            double original_theta = sys.particles[global_i]->theta;

            std::vector<double> freedom_data;
            freedom_data.reserve(sample_pts.size() * 3);

            for (const auto& pt : sample_pts) {
                sys.particles[global_i]->move(pt - sys.particles[global_i]->com);
                size_t worked = 0;

                if (p_type == PARTICLE_TYPE::BUMPY) {
                    double d_theta = 2 * M_PI / Nv / angles_sampled;
                    for (size_t ang_idx = 0; ang_idx < angles_sampled; ++ang_idx) {
                        sys.particles[global_i]->rotate(d_theta);
                        bool interaction_free = std::all_of(neighbor_globals.begin(), neighbor_globals.end(), [&] (size_t gj) {
                            return sys.particles[gj]->get_energy_interaction(sys.particles[global_i].get()) <= 0;
                        });
                        worked += interaction_free;
                    }
                    sys.particles[global_i]->rotate(original_theta - sys.particles[global_i]->theta);
                }
                else {
                    bool interaction_free = std::all_of(neighbor_globals.begin(), neighbor_globals.end(), [&] (size_t gj) {
                        return sys.particles[gj]->get_energy_interaction(sys.particles[global_i].get()) <= 0;
                    });
                    worked = interaction_free;
                }

                freedom_data.push_back(pt.x());
                freedom_data.push_back(pt.y());
                freedom_data.push_back((double) worked / (p_type == PARTICLE_TYPE::BUMPY ? angles_sampled : 1.0));
            }
            sys.particles[global_i]->move(original_pos - sys.particles[global_i]->com);

            std::string particle_dataset_path = decompress_group_path + "/" + format_string(local_i, 6);
            logger.write_dataset(particle_dataset_path, {freedom_data.size() / 3, 3}, freedom_data);
        }
    }

    LOG_INFO << "Simulation finished and data written to " << output_file;
    return 0;
}
