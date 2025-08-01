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

double System::calc_phi () {
    phi = get_particle_area() / L / L;
    return phi;
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

void System::set_energy (double E) {
    double pe = get_pe();
    if (E < pe) {
        LOG_WARNING << "Potential energy is too high; cannot set total energy that low";
        return;
    }
    set_ke(E - pe);
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

void System::set_last_state () {
    last_state_coms.resize(num_p);
    last_state_thetas.resize(num_p);

    for (size_t i = 0; i < num_p; ++i) {
        last_state_coms[i] = particles[i]->com;
        last_state_thetas[i] = particles[i]->theta;
    }
    last_state_L = L;
}

void System::revert_to_last_state () {
    L = last_state_L;
    for (size_t i = 0; i < num_p; ++i) {
        particles[i]->com = last_state_coms[i];
        particles[i]->theta = last_state_thetas[i];
    }

    phi = get_particle_area() / (L * L);

    num_cells_per_dim = static_cast<size_t>(std::floor(L / (r_cut + r_skin)));
    cell_size = L / num_cells_per_dim;
    if (num_cells_per_dim == 0) {
        num_cells_per_dim = 1;
    }
    cell_bins.clear();
    cell_bins.resize(num_cells_per_dim * num_cells_per_dim);

    rebuild_neighbor_list();
}

double System::get_fire_power() {
    double power = 0.0;
    for (size_t i = 0; i < num_p; ++i) if (active[i]) {
        power += particles[i]->force.dot(particles[i]->com_v);
        if (auto* bumpy = dynamic_cast<BumpyParticle*>(particles[i].get())) {
            power += bumpy->torque * bumpy->ang_v;
        }
    }
    return power;
}

void System::set_velocities_to_zero() {
    for (size_t i = 0; i < num_p; ++i) if (active[i]) {
        particles[i]->com_v.setZero();
        if (auto* bumpy = dynamic_cast<BumpyParticle*>(particles[i].get())) {
            bumpy->ang_v = 0.0;
        }
    }
}

void System::fire_velocity_update(double dt_half) {
    for (size_t i = 0; i < num_p; ++i) if (active[i]) {
        particles[i]->com_v += (particles[i]->force / particles[i]->get_mass()) * dt_half;
        if (auto* bumpy = dynamic_cast<BumpyParticle*>(particles[i].get())) {
            // Make sure BumpyParticle has a public 'moi' member or getter
            bumpy->ang_v += (bumpy->torque / bumpy->moi) * dt_half;
        }
    }
}

void System::fire_mix_velocity_and_force(double alpha) {
    // --- Step 1: Calculate global norms, keeping translational and rotational DoFs separate ---

    // Translational norms
    double v_trans_global_sq_norm = 0.0;
    double f_trans_global_sq_norm = 0.0;

    // Rotational norms
    double v_rot_global_sq_norm = 0.0;
    double f_rot_global_sq_norm = 0.0;
    bool has_rotational_dof = false;

    for (size_t i = 0; i < num_p; ++i) if (active[i]) {
        // Accumulate translational components
        v_trans_global_sq_norm += particles[i]->com_v.squaredNorm();
        f_trans_global_sq_norm += particles[i]->force.squaredNorm();

        // Accumulate rotational components if they exist
        if (auto* bumpy = dynamic_cast<BumpyParticle*>(particles[i].get())) {
            has_rotational_dof = true;
            v_rot_global_sq_norm += bumpy->ang_v * bumpy->ang_v;
            f_rot_global_sq_norm += bumpy->torque * bumpy->torque;
        }
    }

    // --- Step 2: Calculate separate mixing ratios ---

    // Translational mixing ratio
    double trans_mixing_ratio = 0.0;
    if (f_trans_global_sq_norm > 1e-40) {
        double v_trans_global_norm = std::sqrt(v_trans_global_sq_norm);
        double f_trans_global_norm = std::sqrt(f_trans_global_sq_norm);
        trans_mixing_ratio = v_trans_global_norm / f_trans_global_norm * alpha;
    }

    // Rotational mixing ratio
    double rot_mixing_ratio = 0.0;
    if (has_rotational_dof && f_rot_global_sq_norm > 1e-40) {
        double v_rot_global_norm = std::sqrt(v_rot_global_sq_norm);
        double f_rot_global_norm = std::sqrt(f_rot_global_sq_norm);
        rot_mixing_ratio = v_rot_global_norm / f_rot_global_norm * alpha;
    }

    // --- Step 3: Apply the mixing rules using the corresponding ratio ---
    for (size_t i = 0; i < num_p; ++i) if (active[i]) {
        // --- Update Translational Part ---
        if (f_trans_global_sq_norm > 1e-40) {
            particles[i]->com_v = particles[i]->com_v * (1.0 - alpha) + particles[i]->force * trans_mixing_ratio;
        } else {
            particles[i]->com_v.setZero(); // Kill velocity if translational force is globally zero
        }

        // --- Update Rotational Part ---
        if (auto* bumpy = dynamic_cast<BumpyParticle*>(particles[i].get())) {
            if (f_rot_global_sq_norm > 1e-40) {
                bumpy->ang_v = bumpy->ang_v * (1.0 - alpha) + bumpy->torque * rot_mixing_ratio;
            } else {
                bumpy->ang_v = 0.0; // Kill angular velocity if torque is globally zero
            }
        }
    }
}

void System::fire_position_update(double dt) {
    for (size_t i = 0; i < num_p; ++i) if (active[i]) {
        Eigen::Vector2d d_trans = particles[i]->com_v * dt;
        verlet_displacements[i] += d_trans; // Important for neighbor list
        particles[i]->move(d_trans);

        if (auto* bumpy = dynamic_cast<BumpyParticle*>(particles[i].get())) {
            double d_theta = bumpy->ang_v * dt;
            bumpy->rotate(d_theta);
        }
    }
}

void System::recalculate_forces_and_energy() {
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
}

void System::fire_minimize (size_t max_steps, double dt_fire, double alpha_start, size_t n_min) {
    // --- FIRE Parameters ---
    double alpha = alpha_start;
    const double dt_fire_max = 10.0 * dt_fire;
    const double dt_fire_min = 1e-3 * dt_fire;
    const double f_inc = 1.1;
    const double f_dec = 0.5;
    const double f_alpha = 0.99;
    size_t N_pos = 0, N_neg = 0;
    const size_t N_neg_max = 10;

    // --- Convergence Parameters ---
    double pe_previous = std::numeric_limits<double>::max();

    // --- Initial State ---
    set_velocities_to_zero();
    recalculate_forces_and_energy();
    double pe_current = get_pe();

    size_t step = 0;
    // --- Main FIRE Loop ---
    for (step = 0; step < max_steps; ++step) {
        // --- 1. Check for Convergence ---
        if (pe_current < PE_tol || std::abs(pe_current / pe_previous - 1.0) < PE_tol) break; // Converged if PE is near zero

        // --- 2. FIRE Power-Based Logic ---
        double power = get_fire_power();
        if (power > 0) { // Moving downhill: accelerate
            N_pos++;
            N_neg = 0;
            if (N_pos > n_min) {
                dt_fire = std::min(dt_fire * f_inc, dt_fire_max);
                alpha *= f_alpha;
            }
        }
        else { // Moving uphill: stop, reverse, and reset
            N_pos = 0;
            N_neg++;
            if (N_neg > N_neg_max) {
                set_velocities_to_zero();
                return;
            }
            dt_fire = std::max(dt_fire * f_dec, dt_fire_min);
            alpha = alpha_start;

            fire_position_update(-dt_fire / 2.0);
            set_velocities_to_zero();
        }

        // --- 3. Modified Velocity-Verlet Integrator ---
        fire_velocity_update(dt_fire / 2.0);      // v(t + dt/2)
        fire_mix_velocity_and_force(alpha);       // FIRE's velocity steering
        fire_position_update(dt_fire);            // x(t + dt)
        recalculate_forces_and_energy();          // F(t + dt)
        fire_velocity_update(dt_fire / 2.0);      // v(t + dt)

        // Update energy for next iteration's convergence check
        pe_previous = pe_current;
        pe_current = get_pe();
    }
    if (step >= max_steps) {
        LOG_WARNING << "Max steps was set too small; FIRE wanted to go longer";
    }

    // --- Cleanup ---
    set_velocities_to_zero();
}

void System::grow_to_phi (double _phi) {
    for (double curr_phi = phi; curr_phi <= _phi + EPS; curr_phi += i_d_phi) {
        rescale(curr_phi);
        fire_minimize();
    }
    rescale(_phi);
    fire_minimize();
}

void System::send_to_jamming() {
    const size_t max_compression_steps = 100000;

    fire_minimize();

    set_last_state();
    double phi_low = phi;
    double phi_high = -1.0;
    double target_phi = phi;

    size_t step = 0;
    for (step = 0; step < max_compression_steps; ++step) {
        fire_minimize();
        double pe_per_particle = get_pe();

        if (pe_per_particle > PE_tol) {
            phi_high = phi;
            revert_to_last_state();
            target_phi = (phi_high + phi_low) / 2.0;
        }
        else {
            set_last_state();
            phi_low = phi;

            if (phi_high > 0) {
                target_phi = (phi_high + phi_low) / 2.0;
            }
            else {
                target_phi = phi + i_d_phi;
            }
        }

        if (phi_high > 0 && std::fabs(phi_high / phi_low - 1.0) < d_phi_min) {
            break;
            // auto is_stable = [&] () -> bool {
            //     double pe = get_pe();
            //     set_last_state();
            //     rescale(calc_phi() + i_d_phi);
            //     fire_minimize();
            //     double pe2 = get_pe();
            //     revert_to_last_state();
            //     return pe2 > 0.1 * i_d_phi * i_d_phi;
            // };
            //
            // bool stable = is_stable();
            // if (stable) break;
            //
            // phi_low = phi;
            // phi_high = -1.0;
            // target_phi = calc_phi();
            // set_last_state();
        }

        rescale(target_phi);
    }
    // LOG_INFO << "Steps took: " << step + 1;

    remove_rattlers();
    // LOG_INFO << num_ratt;
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

void System::make_contact_network () {
    std::vector<Eigen::Triplet<double>> new_contacts;
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        for (size_t j : neighbor_list[i]) if (active[j] && j > i) {
            particles[i]->interact(particles[j].get(), new_contacts);
        }
    }
    adj_contacts.setFromTriplets(new_contacts.begin(), new_contacts.end());

    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        particles[i]->reset_dynamics();
    }
}

size_t System::remove_rattlers () {
    size_t total_rattlers = 0;
    while (true) {
        rebuild_neighbor_list();
        make_contact_network();
        size_t rattlers_found = 0;
        for (size_t i = 0; i < num_p; i++) {
            if (active[i] && particles[i]->rattles()) {
                active[i] = false;
                rattlers_found++;
                active_p--;
            }
        }
        phi = get_particle_area() / L / L;
        total_rattlers += rattlers_found;
        if (rattlers_found == 0) break;
    }
    return total_rattlers;
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
    remove_rattlers();
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
