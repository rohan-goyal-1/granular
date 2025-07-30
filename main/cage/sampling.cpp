#include "include/system.h"
#include "utils/include/hdf5_logger.h"
#include "utils/include/logger.h"
#include "utils/include/argparse.h"
#include "protocol/include/cage_sampling.h"

#include <algorithm>
#include <iomanip>
#include <vector>

std::string format_string (size_t i, size_t sz) {
    std::ostringstream oss;
    oss << std::setw(sz) << std::setfill('0') << i;
    return oss.str();
}

int main (int argc, char** argv) {
    // --- Argument Parsing ---
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
    } catch (const std::exception& e) {
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
    if (mu_eff == 0.0) { p_type = PARTICLE_TYPE::SMOOTH; Nv = 1; }

    auto& std_logger = Logger::get_instance();
    std_logger.set_level(LogLevel::INFO);

    // --- Parameter Logging ---
    LOG_INFO << "=====================================================";
    LOG_INFO << "Starting Iterative Decompression Simulation";
    LOG_INFO << "=====================================================";
    LOG_INFO << "Parameters:";
    LOG_INFO << "  - Output File: " << output_file;
    LOG_INFO << "  - Particle Type: " << (p_type == PARTICLE_TYPE::BUMPY ? "Bumpy" : "Smooth");
    LOG_INFO << "  - Initial Particles (Np): " << N;
    if (p_type == PARTICLE_TYPE::BUMPY) {
        LOG_INFO << "  - Vertices per Particle (Nv): " << Nv;
        LOG_INFO << "  - Angle Samples: " << angles_sampled;
    }
    LOG_INFO << "  - Friction (mu): " << mu_eff;
    LOG_INFO << "  - Timestep (dt): " << dt;
    LOG_INFO << "  - Decompression Target (dphi): " << d_phi;
    LOG_INFO << "  - Decompression Steps: " << num_steps;
    LOG_INFO << "  - Sample Points per Cage: " << points;
    LOG_INFO << "-----------------------------------------------------";

    HDF5Logger logger(output_file);
    System sys(p_type, N, Nv, 0.0, dt, mu_eff);

    // --- Jamming and Hard Rattler Purge ---
    LOG_INFO << "Initializing system: sending to jamming point...";
    sys.send_to_jamming();
    double phi_j = sys.phi;
    LOG_INFO << "System jammed at packing fraction phi_j = " << std::fixed << std::setprecision(6) << phi_j;

    LOG_INFO << "Purging rattlers from system...";
    sys.purge_rattlers(); // This permanently removes rattlers
    LOG_INFO << "Rattlers purged. " << sys.num_p << " particles remain.";
    LOG_INFO << "-----------------------------------------------------";

    // --- Collect centers of mass ---
    std::vector<Eigen::Vector2d> active_particle_COMs;
    active_particle_COMs.reserve(sys.num_p);
    for (size_t i = 0; i < sys.num_p; ++i) {
        active_particle_COMs.push_back(sys.particles[i]->com);
    }

    // --- Writing Initial Jammed State Data ---
    logger.write_attribute("/", "Np", int(sys.num_p));
    logger.write_attribute("/", "Nv", int(Nv));
    logger.write_attribute("/", "mu", mu_eff);
    logger.write_attribute("/", "bumpy", (int)(p_type == PARTICLE_TYPE::BUMPY));

    std::vector<double> flat_dynam_matrix;
    flat_dynam_matrix.reserve(sys.adj_contacts.nonZeros());
    for (int i = 0; i < sys.num_p * Nv; i++) {
        for (int j = 0; j < sys.num_p * Nv; j++) {
            flat_dynam_matrix.push_back(sys.adj_contacts.coeffRef(i, j));
        }
    }
    logger.write_dataset("/adj", {sys.num_p * Nv, sys.num_p * Nv}, flat_dynam_matrix);

    std::vector<double> vertex_data;
    vertex_data.reserve(sys.num_p * 3);
    for (size_t i = 0; i < sys.num_p; ++i) {
        vertex_data.push_back(sys.particles[i]->com.x());
        vertex_data.push_back(sys.particles[i]->com.y());
        vertex_data.push_back(sys.particles[i]->theta);
    }
    logger.write_dataset("/vertices", {sys.num_p, 3}, vertex_data);
    logger.write_attribute("/vertices", "sigma", (p_type == PARTICLE_TYPE::BUMPY ? sys.get_sigma() : 2.0));

    // --- Cage and Sample Point Generation ---
    std::vector<std::vector<Eigen::Vector2d>> all_cages(sys.num_p);
    std::vector<std::vector<Eigen::Vector2d>> all_sample_points(sys.num_p);
    size_t degenerate_cages = 0;

    LOG_INFO << "Generating initial Voronoi cages for " << sys.num_p << " particles...";
    for (size_t i = 0; i < sys.num_p; ++i) {
        if ((i + 1) % std::max<size_t>(1, sys.num_p / 10) == 0 || i == sys.num_p - 1) {
            LOG_INFO << "  ... cage generation " << std::fixed << std::setprecision(1)
                     << (100.0 * (i + 1) / sys.num_p) << "% complete ("
                     << (i + 1) << "/" << sys.num_p << ").";
        }
        all_cages[i] = make_cage(i, active_particle_COMs, sys.L);
        if (!all_cages[i].empty()) {
            all_sample_points[i] = generate_sample_points(all_cages[i], points);
        } else {
            degenerate_cages++;
        }
    }
    LOG_INFO << "Initial cage generation finished. " << sys.num_p - degenerate_cages
             << " valid cages created (" << degenerate_cages << " degenerate).";
    LOG_INFO << "-----------------------------------------------------";


    // --- Neighbor Finding ---
    auto find_relevant_neighbors = [&] (size_t center_idx, const std::vector<Eigen::Vector2d>& cage) {
        std::vector<size_t> neighbors;
        const auto& center_particle = sys.particles[center_idx];
        double R_center = 1 + 0.5 * center_particle->sigma;
        double L = sys.L;

        auto point_to_segment_dist = [] (const Eigen::Vector2d& p, const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
            Eigen::Vector2d ab = b - a; Eigen::Vector2d ap = p - a;
            double t = ab.dot(ap) / ab.squaredNorm(); t = std::clamp(t, 0.0, 1.0);
            return (a + t * ab - p).norm();
        };

        for (size_t j = 0; j < sys.num_p; ++j) {
            if (j == center_idx) continue;
            const auto& pj = sys.particles[j];
            double Rj = 1 + 0.5 * pj->sigma;
            double effective_radius = R_center + Rj;

            Eigen::Vector2d cage_center(0, 0);
            for (const auto& pt : cage) cage_center += pt;
            cage_center /= cage.size();
            Eigen::Vector2d shifted_com = pj->minimum_image(pj->com - cage_center, L) + cage_center;

            for (size_t i = 0; i < cage.size(); ++i) {
                const auto& a = cage[i]; const auto& b = cage[(i + 1) % cage.size()];
                if (point_to_segment_dist(shifted_com, a, b) < effective_radius) {
                    neighbors.push_back(j); break;
                }
            }
        }
        return neighbors;
    };

    // --- Decompression Steps ---
    std::vector<double> delta_phis = {0.0};
    if (num_steps > 1) {
        double log_min = std::min(std::log10(d_phi), -5.0);
        double log_max = std::log10(d_phi);
        double log_step = (log_max - log_min + 1) / (num_steps - 1);
        for (size_t i = 1; i < num_steps; ++i) {
            delta_phis.push_back(std::pow(10.0, log_min + (i-1) * log_step));
        }
    }

    for (size_t decompress_idx = 0; decompress_idx < delta_phis.size(); ++decompress_idx) {
        double current_dphi = delta_phis[decompress_idx];
        double phi_target = phi_j - current_dphi;
        double rat = sys.rescale(phi_target);

        LOG_INFO << "[Step " << (decompress_idx + 1) << "/" << delta_phis.size() << "] "
                 << "Decompressing to phi = " << std::fixed << std::setprecision(6) << phi_target
                 << " (dphi=" << current_dphi << ")";

        std::string decompress_group_path = "/" + format_string(decompress_idx, 6);
        logger.create_group(decompress_group_path);
        logger.write_attribute(decompress_group_path, "L", sys.L);
        logger.write_attribute(decompress_group_path, "dphi", current_dphi);

        for (auto& cage : all_cages) { for (auto& p : cage) p *= rat; }
        for (auto& samples : all_sample_points) { for (auto& p : samples) p *= rat; }

        std::string cages_group_path = decompress_group_path + "/cages";
        logger.create_group(cages_group_path);
        for (size_t i = 0; i < all_cages.size(); ++i) {
            if (!all_cages[i].empty()) {
                std::vector<double> cage_data; cage_data.reserve(all_cages[i].size() * 2);
                for (const auto& pt : all_cages[i]) { cage_data.push_back(pt.x()); cage_data.push_back(pt.y()); }
                logger.write_dataset(cages_group_path + "/" + format_string(i, 6), {all_cages[i].size(), 2}, cage_data);
            }
        }

        for (size_t i = 0; i < sys.num_p; ++i) {
            if ((i + 1) % std::max<size_t>(1, sys.num_p / 10) == 0 || i == sys.num_p - 1) {
                 LOG_INFO << "  ... calculating freedom " << std::fixed << std::setprecision(1)
                          << (100.0 * (i + 1) / sys.num_p) << "% complete (" << (i + 1) << "/" << sys.num_p << ").";
            }

            const auto& sample_pts = all_sample_points[i];
            if (sample_pts.empty()) continue;

            std::vector<size_t> neighbors = find_relevant_neighbors(i, all_cages[i]);
            Eigen::Vector2d original_pos = sys.particles[i]->com;
            double original_theta = sys.particles[i]->theta;
            std::vector<double> freedom_data;
            freedom_data.reserve(sample_pts.size() * 3);

            for (const auto& pt : sample_pts) {
                sys.particles[i]->move(pt - sys.particles[i]->com);
                size_t worked = 0;
                if (p_type == PARTICLE_TYPE::BUMPY) {
                    double d_theta = 2 * M_PI / Nv / angles_sampled;
                    for (size_t ang_idx = 0; ang_idx < angles_sampled; ++ang_idx) {
                        sys.particles[i]->rotate(d_theta);
                        bool interaction_free = std::all_of(neighbors.begin(), neighbors.end(), [&](size_t j) {
                            return sys.particles[j]->get_energy_interaction(sys.particles[i].get()) <= 0;
                        });
                        worked += interaction_free;
                    }
                    sys.particles[i]->rotate(original_theta - sys.particles[i]->theta);
                } else {
                    worked = std::all_of(neighbors.begin(), neighbors.end(), [&](size_t j) {
                        return sys.particles[j]->get_energy_interaction(sys.particles[i].get()) <= 0;
                    });
                }
                freedom_data.push_back(pt.x());
                freedom_data.push_back(pt.y());
                freedom_data.push_back((double)worked / (p_type == PARTICLE_TYPE::BUMPY ? angles_sampled : 1.0));
            }
            sys.particles[i]->move(original_pos - sys.particles[i]->com);
            std::string particle_dataset_path = decompress_group_path + "/" + format_string(i, 6);
            logger.write_dataset(particle_dataset_path, {freedom_data.size() / 3, 3}, freedom_data);
        }
        LOG_INFO << "  Finished step " << (decompress_idx + 1) << ".";
        LOG_INFO << "-----------------------------------------------------";
    }

    LOG_INFO << "=====================================================";
    LOG_INFO << "Success! Simulation finished.";
    LOG_INFO << "All data written to " << output_file;
    LOG_INFO << "=====================================================";
    return 0;
}
