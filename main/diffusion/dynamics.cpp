#include "include/system.h"
#include "utils/include/logger.h"
#include "utils/include/argparse.h"
#include "utils/include/hdf5_logger.h"

std::string format_string (size_t i, size_t sz) {
    std::ostringstream oss;
    oss << std::setw(sz) << std::setfill('0') << i;
    return oss.str();
}

const double KE = 1e-5;

int main (int argc, char** argv) {
    // --- Argument Parsing ---
    ArgParser parser("compression_dynamics");
    parser.add_argument<bool>("--bumpy", "-b").help("Use bumpy particles.").default_value(false);
    parser.add_argument<size_t>("--nv", "-v").help("Number of vertices for bumpy particles.").default_value(1);
    parser.add_argument<std::string>("--output", "-o").help("The path to the file to write results to.").default_value("dyn.h5");
    parser.add_argument<size_t>("--np", "-n").help("The number of particles to use.").required();
    parser.add_argument<double>("--phi_s").help("Starting packing fraction.").required();
    parser.add_argument<double>("--phi_e").help("Ending packing fraction.").required();
    parser.add_argument<size_t>("--steps").help("Number of compression steps from phi_s to phi_e.").required();
    parser.add_argument<double>("--mu", "-m").help("Effective friction coefficient for the particles.").required();
    parser.add_argument<double>("--dt").help("The timestep for the simulation.").default_value(1e-3);
    parser.add_argument<double>("--dt_min").help("The time interval for logging system state.").required();
    parser.add_argument<double>("--dt_max").help("The total simulation time to run at each compression step.").required();

    try {
        parser.parse(argc, argv);
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Use --help for more information." << std::endl;
        return 1;
    }

    // --- Variable Initialization from Arguments ---
    PARTICLE_TYPE p_type = parser.get<bool>("--bumpy") ? PARTICLE_TYPE::BUMPY : PARTICLE_TYPE::SMOOTH;
    std::string output_file = parser.get<std::string>("--output");
    size_t Nv = parser.get<size_t>("--nv");
    size_t N = parser.get<size_t>("--np");
    double dt = parser.get<double>("--dt");
    double phi_s = parser.get<double>("--phi_s");
    double phi_e = parser.get<double>("--phi_e"); // Fixed bug: was reading --phi_s again
    size_t num_steps = parser.get<size_t>("--steps");
    double mu_eff = parser.get<double>("--mu");
    double dt_min = parser.get<double>("--dt_min");
    double dt_max = parser.get<double>("--dt_max");
    if (mu_eff == 0.0) { p_type = PARTICLE_TYPE::SMOOTH; Nv = 1; }

    auto& std_logger = Logger::get_instance();
    std_logger.set_level(LogLevel::INFO);

    // --- Parameter Logging ---
    LOG_INFO << "=====================================================";
    LOG_INFO << "Starting Iterative Compression Dynamics Simulation";
    LOG_INFO << "=====================================================";
    LOG_INFO << "Parameters:";
    LOG_INFO << "  - Output File: " << output_file;
    LOG_INFO << "  - Particle Type: " << (p_type == PARTICLE_TYPE::BUMPY ? "Bumpy" : "Smooth");
    LOG_INFO << "  - Initial Particles (Np): " << N;
    if (p_type == PARTICLE_TYPE::BUMPY) {
        LOG_INFO << "  - Vertices per Particle (Nv): " << Nv;
    }
    LOG_INFO << "  - Friction (mu): " << mu_eff;
    LOG_INFO << "  - Starting Phi (phi_s): " << phi_s;
    LOG_INFO << "  - Ending Phi (phi_e): " << phi_e;
    LOG_INFO << "  - Compression Steps: " << num_steps;
    LOG_INFO << "  - Timestep (dt): " << dt;
    LOG_INFO << "  - Logging Interval (dt_min): " << dt_min;
    LOG_INFO << "  - Simulation Duration (dt_max): " << dt_max;
    LOG_INFO << "-----------------------------------------------------";

    HDF5Logger logger(output_file);
    System sys(p_type, N, Nv, 0.0, dt, mu_eff);

    // --- Initial System Preparation ---
    LOG_INFO << "Initializing system: growing to starting density...";
    sys.grow_to_phi(phi_s);
    LOG_INFO << "System packed to phi = " << phi_s;

    // --- Writing Global Metadata to HDF5 ---
    logger.write_attribute("/", "Np", int(sys.num_p));
    logger.write_attribute("/", "Nv", int(Nv));
    logger.write_attribute("/", "mu", mu_eff);
    logger.write_attribute("/", "bumpy", (int)(p_type == PARTICLE_TYPE::BUMPY));

    // --- Compression Schedule ---
    std::vector<double> phi_schedule;
    if (num_steps > 1) {
        double phi_step_size = (phi_e - phi_s) / (num_steps - 1);
        for (size_t i = 0; i < num_steps; ++i) {
            phi_schedule.push_back(phi_s + i * phi_step_size);
        }
    } else if (num_steps == 1) {
        phi_schedule.push_back(phi_e); // If only one step, go straight to the end
    }
    // If num_steps is 0, the schedule is empty and the main loop is skipped.

    // --- Compression and Dynamics Loop ---
    for (size_t i = 0; i < phi_schedule.size(); ++i) {
        double phi_target = phi_schedule[i];

        // Compress the system to the target packing fraction
        sys.rescale(phi_target);

        LOG_INFO << "[Step " << (i + 1) << "/" << phi_schedule.size() << "] "
                 << "Compressed to phi = " << std::fixed << std::setprecision(6) << phi_target;

        // Inject kinetic energy to start the dynamics
        sys.set_energy(KE);
        LOG_INFO << "Kinetic energy added. KE = " << KE << ".";

        // Create a unique group in the HDF5 file for this phi value
        std::string group_path = "/" + format_string(i, 6);
        logger.create_group(group_path);
        logger.write_attribute(group_path, "L", sys.L);
        logger.write_attribute(group_path, "phi", phi_target);

        // --- Run Dynamics Simulation ---
        LOG_INFO << "  Running dynamics for T = " << dt_max << "...";
        size_t total_steps = static_cast<size_t>(dt_max / dt);
        size_t log_interval_steps = (dt_min > 0) ? static_cast<size_t>(dt_min / dt) : 1;
        if (log_interval_steps == 0) log_interval_steps = 1;

        std::vector<double> trajectory_data;
        size_t num_frames = total_steps / log_interval_steps;
        trajectory_data.reserve(num_frames * sys.num_p * 3); // Pre-allocate memory

        for (size_t step = 0; step < total_steps; ++step) {
            sys.update(); // Advance simulation by one timestep `dt`

            if (step % log_interval_steps == 0) {
                for (size_t p_idx = 0; p_idx < sys.num_p; ++p_idx) {
                    trajectory_data.push_back(sys.particles[p_idx]->com.x());
                    trajectory_data.push_back(sys.particles[p_idx]->com.y());
                    if (auto* p = dynamic_cast<BumpyParticle*>(sys.particles[p_idx].get())) {
                        trajectory_data.push_back(p->theta);
                    } else {
                        trajectory_data.push_back(0.0);
                    }
                }
            }
        }

        // Write the collected trajectory data for this compression step
        if (!trajectory_data.empty()) {
            size_t actual_frames = trajectory_data.size() / (sys.num_p * 3);
            logger.write_dataset(group_path + "/positions", {actual_frames, sys.num_p, 3}, trajectory_data);
        }
        LOG_INFO << "  Finished step " << (i + 1) << ". Logged " << (trajectory_data.size() / (sys.num_p * 3)) << " frames.";
        LOG_INFO << "-----------------------------------------------------";
    }

    LOG_INFO << "=====================================================";
    LOG_INFO << "Success! Simulation finished.";
    LOG_INFO << "All data written to " << output_file;
    LOG_INFO << "=====================================================";

    return 0;
}
