#include "protocol/include/diffusion.h"
#include "utils/include/logger.h"

// std::vector<std::vector<std::pair<double, double>>> find_diffusion (System* sys, double start, double end, double time, size_t count, std::vector<double>& phis) {
//     const double DT = sys->dt;
//     const size_t NP = sys->particles.size();
//     size_t total_frames = static_cast<size_t>(time / DT);
//
//     // Logarithmically spaced Δt values
//     std::vector<double> delta_ts; delta_ts.reserve(count);
//     double log_min = std::log10(start);
//     double log_max = std::log10(end);
//     double step = (log_max - log_min) / (count - 1);
//
//     for (size_t i = 0; i < count; i++) {
//         delta_ts.push_back(std::pow(10.0, log_min + i * step));
//     }
//
//     size_t phi_count = phis.size();
//     std::vector<std::vector<std::pair<double, double>>> results; results.reserve(phi_count);
//     for (size_t i = 0; i < phis.size(); i++) {
//         LOG_INFO << "Starting diffusion run for φ = " << phis[i];
//
//         sys->grow_to_phi(phis[i]);
//         sys->set_temp(sys->E_i - sys->get_pe());
//         results.push_back(std::vector<std::pair<double, double>>());
//         auto& result = results.back();
//
//         std::vector<double> average_msd(count, 0.0);
//         std::vector<std::vector<Eigen::VectorXd>> last_pos(count, std::vector<Eigen::VectorXd>(NP));
//
//         // Save initial positions
//         for (size_t i = 0; i < count; i++) {
//             for (size_t j = 0; j < NP; j++) {
//                 last_pos[i][j] = sys->particles[j].com;
//             }
//         }
//
//         for (size_t cur_frame = 0; cur_frame <= total_frames; cur_frame++) {
//             for (size_t i = 0; i < count; ++i) {
//                 size_t frame_freq = static_cast<size_t>(delta_ts[i] / DT);
//                 if (cur_frame % frame_freq == 0 && frame_freq > 0) {
//                     double msd = 0.0;
//                     for (size_t j = 0; j < NP; j++) {
//                         Eigen::VectorXd delta = sys->particles[j].com - last_pos[i][j];
//                         msd += delta.squaredNorm();
//                         last_pos[i][j] = sys->particles[j].com;
//                     }
//                     msd /= NP;
//
//                     average_msd[i] = ((cur_frame / frame_freq) * average_msd[i] + msd) / (cur_frame / frame_freq + 1);
//                 }
//             }
//
//             sys->update();
//         }
//
//         // Store log10(Δt), log10(MSD)
//         for (size_t i = 0; i < count; ++i) {
//             if (average_msd[i] == 0.0) continue;
//             result.emplace_back(std::log10(delta_ts[i]), std::log10(average_msd[i]));
//         }
//
//         LOG_INFO << "Finished diffusion run for φ = " << phis[i];
//     }
//
//     return results;
// }
