#include "include/system.h"
#include "utils/include/logger.h"

System::System (PARTICLE_TYPE particle_type, size_t num_p, size_t num_v, double E_i, double dt, double mu, double phi)
    : adj_contacts(num_v * num_p, num_v * num_p),
      phi(phi),
      mu(mu),
      L(1.0),
      dt(dt),
      gen((std::random_device{}())),
      dist_angle(0.0, 2 * M_PI),
      dist_pos(0.0, 1.0),
      dist_sign(0, 1),
      num_p(num_p),
      active_p(num_p),
      num_v(num_v),
      E_i(E_i),
      active(num_p, true),
      particle_type(particle_type)
{
    if (particle_type == PARTICLE_TYPE::BUMPY && num_v < 2) {
        throw std::runtime_error("System CONSTRUCTOR: Bumpy particles require at least 2 vertices.");
    }

    adj_contacts.setZero();

    // --- STEP 1: Create particles ---
    particles.reserve(num_p);
    for (size_t i = 0; i < num_p; ++i) {
        switch (particle_type) {
            case PARTICLE_TYPE::SMOOTH:
                particles.push_back(std::make_unique<SmoothParticle>(this, i));
                break;
            case PARTICLE_TYPE::BUMPY:
                particles.push_back(std::make_unique<BumpyParticle>(this, num_v, mu, i));
                break;
            default:
                throw std::runtime_error("System CONSTRUCTOR: Unknown particle type.");
        }
    }

    // --- STEP 2: Now that particles exist, calculate sigma and setup cell list ---
    r_cut = 2.5 * (1.0 + get_sigma());
    r_skin = 0.5 * get_sigma();
    num_cells_per_dim = static_cast<size_t>(std::floor(L / (r_cut + r_skin)));
    cell_size = L / num_cells_per_dim;

    cell_bins = std::vector<std::vector<size_t>>(num_cells_per_dim * num_cells_per_dim);
    neighbor_list = std::vector<std::vector<size_t>>(num_p);
    verlet_displacements = std::vector<Eigen::Vector2d>(num_p, Eigen::Vector2d::Zero());

    // --- STEP 3: Generate random non-overlapping config ---
    random_nonoverlap_config();
    grow_to_phi(phi);
    set_ke(E_i - get_pe());
}

void System::random_nonoverlap_config () {
    for (size_t attempt = 0; attempt < MAX_ATTEMPTS; ++attempt) {
        for (auto& p : particles) {
            p->randomize_position();
        }

        rescale(INIT_PHI);

        if (get_pe() == 0.0) return;
    }

    throw std::runtime_error("Failed to generate initial non-overlapping configuration.");
}

double System::get_sigma () {
    for (size_t i = 0; i < num_p; i++) {
        if (active[i]) return particles[i]->sigma;
    }
    return 0;
}

inline Eigen::Vector2d minimum_image (const Eigen::Vector2d& vec, double L) {
    return vec - L * (vec / L).array().round().matrix();
}

size_t System::get_flat_cell_index (const Eigen::Vector2d& pos) const {
    Eigen::Vector2d wrapped = pos;
    wrapped[0] = fmod(wrapped[0], L);
    wrapped[1] = fmod(wrapped[1], L);
    if (wrapped[0] < 0) wrapped[0] += L;
    if (wrapped[1] < 0) wrapped[1] += L;

    size_t x = static_cast<size_t>(wrapped[0] / cell_size);
    size_t y = static_cast<size_t>(wrapped[1] / cell_size);

    if (x >= num_cells_per_dim) x = num_cells_per_dim - 1;
    if (y >= num_cells_per_dim) y = num_cells_per_dim - 1;

    return y * num_cells_per_dim + x;
}

std::vector<size_t> System::get_neighboring_cells (size_t flat_index) const {
    std::vector<size_t> neighbors;
    size_t cx = flat_index % num_cells_per_dim;
    size_t cy = flat_index / num_cells_per_dim;
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            size_t nx = (cx + dx + num_cells_per_dim) % num_cells_per_dim;
            size_t ny = (cy + dy + num_cells_per_dim) % num_cells_per_dim;
            neighbors.push_back(ny * num_cells_per_dim + nx);
        }
    }

    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

    return neighbors;
}

void System::rebuild_cell_list () {
    for (auto& bin : cell_bins) bin.clear();

    for (size_t i = 0; i < num_p; ++i) {
        if (!active[i]) continue;
        size_t cell_idx = get_flat_cell_index(particles[i]->com);
        cell_bins[cell_idx].push_back(i);
    }
}

void System::rebuild_neighbor_list () {
    rebuild_cell_list();

    for (size_t i = 0; i < num_p; ++i) {
        neighbor_list[i].clear();
        verlet_displacements[i].setZero();
    }

    for (size_t i = 0; i < num_p; ++i) {
        if (!active[i]) continue;

        size_t my_cell = get_flat_cell_index(particles[i]->com);
        auto neighbor_cells = get_neighboring_cells(my_cell);

        for (size_t cid : neighbor_cells) {
            for (size_t j : cell_bins[cid]) {
                if (j <= i) continue;
                if (!active[j]) continue;

                Eigen::Vector2d dx = minimum_image(particles[j]->com - particles[i]->com, L);

                if (dx.squaredNorm() < (r_cut + r_skin) * (r_cut + r_skin)) {
                    neighbor_list[i].push_back(j);
                    // neighbor_list[j].push_back(i);
                }
            }
        }
    }
}
bool System::neighbor_list_invalid () const {
    double max_disp_sq = (0.5 * r_skin) * (0.5 * r_skin);
    for (size_t i = 0; i < num_p; ++i) {
        if (!active[i]) continue;
        if (verlet_displacements[i].squaredNorm() > max_disp_sq) {
            return true;
        }
    }
    return false;
}
double System::rescale (double _phi) {
    phi = _phi;
    double old_L = L;
    double p_area = get_particle_area();

    L = std::sqrt(p_area / phi);

    double ratio = L / old_L;
    for (auto& p : particles) {
        p->rescale_ratio(ratio);
    }

    num_cells_per_dim = static_cast<size_t>(std::floor(L / (r_cut + r_skin)));
    cell_size = L / num_cells_per_dim;
    if (num_cells_per_dim == 0) {
        num_cells_per_dim = 1;
    }
    cell_bins.clear();
    cell_bins.resize(num_cells_per_dim * num_cells_per_dim);
    rebuild_neighbor_list();

    return ratio;
}

void System::set_ke (double ke) {
    for (auto& p : particles) {
        p->set_ke(ke);
    }
}

inline bool System::is_relaxed () {
    return get_ke() < KE_tol;
}

void System::relax () {
    relaxing = true;
    for (int i = 0; i < MIN_STEPS; i++) {
        update();
    }
    while (!is_relaxed()) {
        update();
    }
    relaxing = false;
}

void System::fire_minimize (size_t max_steps, double dt_max, double alpha_start, size_t n_min, double force_tol) {
    // to check
    double alpha = alpha_start;
    double dt_fire = 0.1 * dt_max;
    size_t N_pos = 0;

    std::vector<Eigen::Vector2d> vels(num_p, Eigen::Vector2d::Zero());

    // --- Initial force computation ---
    N_c = 0;
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        particles[i]->reset_dynamics();
    }

    if (neighbor_list_invalid())
        rebuild_neighbor_list();

    std::vector<Eigen::Triplet<double>> new_contacts;
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        for (size_t j : neighbor_list[i]) if (active[j] && j > i) {
            particles[i]->interact(particles[j].get(), new_contacts);
        }
    }
    adj_contacts.setFromTriplets(new_contacts.begin(), new_contacts.end());

    // --- Main FIRE loop ---
    for (size_t step = 0; step < max_steps; ++step) {
        double P = 0.0;
        double force_norm_sq = 0.0;
        double velocity_norm_sq = 0.0;

        for (size_t i = 0; i < num_p; i++) if (active[i]) {
            P += vels[i].dot(particles[i]->force);
            force_norm_sq += particles[i]->force.squaredNorm();
            velocity_norm_sq += vels[i].squaredNorm();
        }

        double force_norm = std::sqrt(force_norm_sq);

        if (force_norm < force_tol) break;

        if (P <= 0) {
            N_pos = 0;
            dt_fire *= 0.5;
            alpha = alpha_start;
            for (size_t i = 0; i < num_p; ++i) if (active[i]) {
                vels[i].setZero();
            }
        } else {
            N_pos++;
            if (N_pos > n_min) {
                dt_fire = std::min(dt_fire * 1.1, dt_max);
                alpha *= 0.99;
            }

            double total_vel_mag = std::sqrt(velocity_norm_sq);
            for (size_t i = 0; i < num_p; ++i) if (active[i]) {
                if (force_norm > 1e-12) {
                    vels[i] = (1.0 - alpha) * vels[i] +
                              alpha * (particles[i]->force / force_norm) * total_vel_mag;
                } else {
                    vels[i].setZero();
                }
            }
        }

        // --- Position update ---
        for (size_t i = 0; i < num_p; i++) if (active[i]) {
            Eigen::Vector2d acc_i = particles[i]->force / particles[i]->get_mass();
            Eigen::Vector2d d_trans = vels[i] * dt_fire + 0.5 * acc_i * dt_fire * dt_fire;
            verlet_displacements[i] += d_trans;  // track for neighbor list rebuild
            particles[i]->move(d_trans);
            vels[i] += 0.5 * acc_i * dt_fire;
        }

        // --- Recompute forces ---
        N_c = 0;
        for (size_t i = 0; i < num_p; i++) if (active[i]) {
            particles[i]->reset_dynamics();
        }

        if (neighbor_list_invalid()) {
            rebuild_neighbor_list();
        }

        std::vector<Eigen::Triplet<double>> new_contacts;
        for (size_t i = 0; i < num_p; i++) if (active[i]) {
            for (size_t j : neighbor_list[i]) if (active[j] && j > i) {
                particles[i]->interact(particles[j].get(), new_contacts);
            }
        }
        adj_contacts.setFromTriplets(new_contacts.begin(), new_contacts.end());

        // --- Velocity update ---
        for (size_t i = 0; i < num_p; ++i) if (active[i]) {
            Eigen::Vector2d new_acc = particles[i]->force / particles[i]->get_mass();
            vels[i] += 0.5 * new_acc * dt_fire;
        }
    }
}

void System::grow_to_phi (double _phi) {
    for (double curr_phi = phi; curr_phi <= _phi + EPS; curr_phi += i_d_phi) {
        rescale(curr_phi);
        fire_minimize();
    }
    rescale(_phi);
    fire_minimize();
}

void System::send_to_jamming () {
    while (true) {
        fire_minimize();

        double d_phi = i_d_phi;
        double pe = get_pe();
        double last_pe = pe;

        while (pe < PE_tol || pe > 1.01 * PE_tol) {
            if (
                (pe < PE_tol && last_pe > 1.01 * PE_tol) ||
                (pe > 1.01 * PE_tol && last_pe < PE_tol)
            ) {
                if (d_phi * 0.5 < d_phi_min) {
                    if (pe < PE_tol) {
                        rescale(phi + d_phi);
                        fire_minimize();
                    }
                    break;
                }
                else {
                    d_phi *= 0.5;
                }
            }

            if (pe < PE_tol)
                rescale(phi + d_phi);
            else
                rescale(phi - d_phi);

            fire_minimize();
            last_pe = pe;
            pe = get_pe();
        }

        size_t rattlers = remove_rattlers();
        if (rattlers == 0) {
            break;
        }
    }
}
void System::update () {
    N_c = 0;
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        auto old_pos = particles[i]->com;
        particles[i]->reset_dynamics();
        particles[i]->update();
        verlet_displacements[i] += particles[i]->com - old_pos;
    }

    if (neighbor_list_invalid()) {
        rebuild_neighbor_list();
    }

    std::vector<Eigen::Triplet<double>> new_contacts;
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        for (size_t j : neighbor_list[i]) if (active[j] && j > i) {
            particles[i]->interact(particles[j].get(), new_contacts);
        }
    }
    adj_contacts.setFromTriplets(new_contacts.begin(), new_contacts.end());

    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        if (relaxing) particles[i]->apply_drag(RELAX_KD);
        particles[i]->integrate();
    }

    time += dt;
}

double System::get_ke () {
    double ke = 0.0;
    for (auto& p : particles)  {
        ke += p->get_ke();
    }
    return ke / num_p;
}

double System::get_pe () {
    double pe = 0.0;

    if (neighbor_list_invalid()) {
        rebuild_neighbor_list();
    }

    for (size_t i = 0; i < num_p; ++i) {
        if (!active[i]) continue;
        for (size_t j : neighbor_list[i]) {
            if (!active[j] || j <= i) continue;
            pe += particles[i]->get_energy_interaction(particles[j].get());
        }
    }

    return pe / num_p;
}

size_t System::remove_rattlers () {
    size_t rattlers_found = 0;
    for (size_t i = 0; i < num_p; i++) {
        if (active[i] && particles[i]->rattles()) {
            active[i] = false;
            rattlers_found++;
            active_p--;
        }
    }
    phi = get_particle_area() / L / L;
    return rattlers_found;
}

double System::get_particle_area () {
    double area = 0.0;
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        area += particles[i]->area;
    }
    return area;
}

std::vector<double> System::get_verts () const {
    std::vector<double> vertex_data;

    size_t verts_per_particle = (particle_type == PARTICLE_TYPE::SMOOTH) ? 1 : num_v;
    vertex_data.reserve(active_p * verts_per_particle * 2);

    for (size_t i = 0; i < num_p; i++) {
        if (!active[i]) continue; // Skip inactive particles

        if (auto* smooth = dynamic_cast<SmoothParticle*>(particles[i].get())) {
            vertex_data.push_back(smooth->com[0]);
            vertex_data.push_back(smooth->com[1]);
        }
        else if (auto* bumpy = dynamic_cast<BumpyParticle*>(particles[i].get())) {
            for (size_t v = 0; v < num_v; v++) {
                vertex_data.push_back(bumpy->verts[v][0]);
                vertex_data.push_back(bumpy->verts[v][1]);
            }
        }
    }

    return vertex_data;
}
size_t System::purge_rattlers () {
    // Determine how many particles need to be removed.
    const size_t particles_to_remove_count = num_p - active_p;

    if (particles_to_remove_count == 0) {
        return 0; // Nothing to do, all particles are active.
    }

    // --- STEP 1: Create a map from old particle indices to new ones ---
    // This map is essential for rebuilding the contact matrix correctly.
    std::vector<size_t> old_to_new_map(num_p);
    size_t new_idx = 0;
    for (size_t i = 0; i < num_p; ++i) {
        if (active[i]) {
            old_to_new_map[i] = new_idx++;
        }
    }

    // --- STEP 2: Preserve the contacts between particles that are being kept ---
    std::vector<Eigen::Triplet<double>> kept_contacts;
    kept_contacts.reserve(adj_contacts.nonZeros());

    // Iterate through all non-zero elements of the old sparse matrix.
    for (int k = 0; k < adj_contacts.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(adj_contacts, k); it; ++it) {
            size_t old_p_row = it.row() / num_v;
            size_t old_p_col = it.col() / num_v;

            // Keep a contact only if both particles involved were active.
            if (active[old_p_row] && active[old_p_col]) {
                size_t new_p_row = old_to_new_map[old_p_row];
                size_t new_p_col = old_to_new_map[old_p_col];

                size_t new_row = new_p_row * num_v + (it.row() % num_v);
                size_t new_col = new_p_col * num_v + (it.col() % num_v);

                kept_contacts.emplace_back(new_row, new_col, it.value());
            }
        }
    }

    // --- STEP 3: Create the new, smaller vector of particles ---
    std::vector<std::unique_ptr<Particle>> kept_particles;
    kept_particles.reserve(active_p); // Reserve space for the particles we are keeping.
    for (size_t i = 0; i < num_p; ++i) {
        if (active[i]) {
            kept_particles.push_back(std::move(particles[i]));
        }
    }
    particles = std::move(kept_particles);

    // --- STEP 4: Update all system counters and data structures ---
    num_p = active_p; // The new total number of particles is the old active count.

    // Resize and rebuild the contact matrix from the preserved contacts.
    adj_contacts.resize(num_v * num_p, num_v * num_p);
    adj_contacts.setFromTriplets(kept_contacts.begin(), kept_contacts.end());

    // Re-initialize lists for the new size. All remaining particles are now active.
    active = std::vector<bool>(num_p, true);
    neighbor_list = std::vector<std::vector<size_t>>(num_p);
    verlet_displacements = std::vector<Eigen::Vector2d>(num_p, Eigen::Vector2d::Zero());

    // Re-assign particle IDs to be contiguous from 0 to num_p-1.
    for (size_t i = 0; i < num_p; ++i) {
        particles[i]->id = i;
    }

    // The neighbor list must be rebuilt from scratch with the new indices.
    rebuild_neighbor_list();

    phi = get_particle_area() / L / L;

    return particles_to_remove_count;
}
