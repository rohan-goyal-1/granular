// TODO: Fix temp calculations

#include "include/system.h"
#include "utils/include/logger.h"
#include "include/cell_list.h"

std::atomic<int> System::curr_id{0};

System::System (PARTICLE_TYPE particle_type, size_t num_p, size_t num_v, double E_i, double dt, double mu, double phi) :
    adj_contacts(num_v * num_p, num_v * num_p),
    phi(phi),
    mu(mu),
    dt(dt),
    gen((std::random_device{}())),
    dist_angle(std::uniform_real_distribution<>(0.0, 2 * M_PI)),
    dist_pos(std::uniform_real_distribution<>(0.0, 1)),
    dist_sign(std::uniform_int_distribution<>(0, 1)),
    num_p(num_p),
    num_v(num_v),
    E_i(E_i),
    active(num_p, true),
    particle_type(particle_type),
    verlet_displacements(num_p)
{
    if (particle_type == PARTICLE_TYPE::BUMPY && num_v < 2) {
        throw std::runtime_error("System CONSTRUCTOR: Bumpy particles require at least 2 vertices.");
    }

    adj_contacts.setZero();
    id = get_id();

    for (size_t curr_attempt = 0; curr_attempt < MAX_ATTEMPTS; curr_attempt++) {
        particles.clear();
        particles.reserve(num_p);

        for (size_t i = 0; i < num_p; i++) {
            switch (particle_type) {
                case PARTICLE_TYPE::SMOOTH:
                    particles.push_back(std::make_unique<SmoothParticle>(this, 1.0, i));
                    break;
                case PARTICLE_TYPE::BUMPY:
                    particles.push_back(std::make_unique<BumpyParticle>(this, num_v, mu, INIT_PHI / num_p, i));
                    break;
                default:
                    throw std::runtime_error("System CONSTRUCTOR: Unknown particle type.");
            }
        }

        rescale(INIT_PHI);
        if (get_pe() == 0.0) break;
    }

    grow_to_phi(phi);

    double interaction_radius = 2.0 * get_radius();
    verlet_skin_radius = 0.4 * interaction_radius;
    double verlet_cutoff = interaction_radius + verlet_skin_radius;
    verlet_cutoff_sq = verlet_cutoff * verlet_cutoff;

    cell_list = std::make_unique<CellList>(this, verlet_cutoff);
    // build_verlet_list();

    if (E_i < get_pe()) {
        throw std::runtime_error("System CONSTRUCTOR: Energy minimum is greater than requested initial energy. Try a lower packing fraction.");
    }

    set_temp(E_i - get_pe());
}

System::~System() = default;

void System::apply_pbc_to_vector (Eigen::Vector2d& vec) const {
    vec.x() -= 1 * std::round(vec.x() / 1);
    vec.y() -= 1 * std::round(vec.y() / 1);
}

Eigen::Vector2d System::get_min_dist_vec (const Particle* p1, const Particle* p2) const {
    Eigen::Vector2d dr = p2->com - p1->com;
    apply_pbc_to_vector(dr);
    return dr;
}

void System::build_verlet_list () {
    cell_list->build();
    verlet_list.clear();

    std::vector<int> neighbor_cells;
    for (size_t i = 0; i < num_p; i++) {
        if (!active[i]) continue;

        int cell_idx = cell_list->get_cell_index(particles[i]->com);
        cell_list->get_neighbor_cells(cell_idx, neighbor_cells);

        for (int neighbor_cell_idx : neighbor_cells) {
            for (size_t j : cell_list->get_particles_in_cell(neighbor_cell_idx)) {
                if (i >= j || !active[j]) continue;

                Eigen::Vector2d dr = get_min_dist_vec(particles[i].get(), particles[j].get());
                if (dr.squaredNorm() < verlet_cutoff_sq) {
                    verlet_list.emplace_back(i, j);
                }
            }
        }
    }

    std::fill(verlet_displacements.begin(), verlet_displacements.end(), Eigen::Vector2d::Zero());
    rebuild_verlet_list = false;
}

void System::check_verlet_rebuild () {
    if (rebuild_verlet_list) return;

    double max_disp_sq = 0.0;
    for (size_t i = 0; i < num_p; i++) {
        if (active[i]) {
            max_disp_sq = std::max(max_disp_sq, verlet_displacements[i].squaredNorm());
        }
    }

    if (max_disp_sq * 4.0 > verlet_skin_radius * verlet_skin_radius) {
        rebuild_verlet_list = true;
    }
}

double System::get_radius () {
    return particles[0]->radius;
}

double System::get_sigma () {
    return particles[0]->sigma;
}

void System::rescale (double _phi) {
    // TODO: Dont rescale based on phi but on radius increments
    phi = _phi;
    double area_part = phi / num_p;
    for (auto& p : particles) {
        p->rescale(area_part);
    }
    fire_minimize();

    // double interaction_radius = 2.5 * get_sigma();
    // verlet_skin_radius = 0.4 * interaction_radius;
    // double verlet_cutoff = interaction_radius + verlet_skin_radius;
    // verlet_cutoff_sq = verlet_cutoff * verlet_cutoff;
    //
    // cell_list = std::make_unique<CellList>(this, verlet_cutoff);
    // rebuild_verlet_list = true;
}

void System::set_temp (double temp) {
    for (auto& p : particles) {
        p->set_ke(temp);
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
    double alpha = alpha_start;
    double dt_fire = 0.1 * dt_max;
    size_t N_pos = 0;

    std::vector<Eigen::Vector2d> vels(num_p);
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        vels[i] = Eigen::Vector2d(0, 0);
    }

    N_c = 0;
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        particles[i]->reset_dynamics();
    }
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        auto &p1 = particles[i];
        for (size_t j = i + 1; j < num_p; j++) if (active[j]) {
            auto &p2 = particles[j];
            p1->interact(p2.get());
        }
    }

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

        if (force_norm < force_tol) {
            break;
        }

        if (P <= 0) {
            N_pos = 0;
            dt_fire *= 0.5;
            alpha = alpha_start;
            for (size_t i = 0; i < num_p; ++i) if (active[i]) {
                vels[i].setZero();
            }
        }
        else {
            N_pos++;
            if (N_pos > n_min) {
                dt_fire = std::min(dt_fire * 1.1, dt_max);
                alpha *= 0.99;
            }

            double total_system_velocity_magnitude = std::sqrt(velocity_norm_sq);
            for (size_t i = 0; i < num_p; ++i) if (active[i]) {
                if (force_norm > 1e-12) {
                     vels[i] = (1.0 - alpha) * vels[i] + alpha * (particles[i]->force / force_norm) * total_system_velocity_magnitude;
                }
                else {
                     vels[i] = Eigen::Vector2d(0, 0);
                }
            }
        }

        for (size_t i = 0; i < num_p; i++) if (active[i]) {
            Eigen::Vector2d acc_i = particles[i]->force / num_v;
            auto d_trans = vels[i] * dt_fire + 0.5 * acc_i * dt_fire * dt_fire;
            particles[i]->move(d_trans);
            vels[i] += 0.5 * acc_i * dt_fire;
        }

        N_c = 0;
        for (size_t i = 0; i < num_p; i++) if (active[i]) {
            particles[i]->reset_dynamics();
        }
        for (size_t i = 0; i < num_p; i++) if (active[i]) {
            auto &p1 = particles[i];
            for (size_t j = i + 1; j < num_p; j++) if (active[j]) {
                auto &p2 = particles[j];
                p1->interact(p2.get());
            }
        }

        for (size_t i = 0; i < num_p; ++i) if (active[i]) {
            Eigen::Vector2d new_acc_i = particles[i]->force / num_v;
            vels[i] += 0.5 * new_acc_i * dt_fire;
        }
    }
}

void System::grow_to_phi (double _phi) {
    for (double curr_phi = phi; curr_phi <= _phi + EPS; curr_phi += i_d_phi) {
        rescale(curr_phi);
    }
    rescale(_phi);
}

void System::send_to_jamming () {
    fire_minimize();

    double d_phi = i_d_phi;

    double pe = get_pe();
    double last_pe = pe;

    while (pe != PE_tol) {
        if (
            (pe < PE_tol && last_pe > PE_tol) ||
            (pe > PE_tol && last_pe < PE_tol)
        ) {
            if (d_phi * 0.5 < d_phi_min) {
                if (pe < PE_tol) {
                    rescale(phi + d_phi);
                }
                return;
            }
            else {
                d_phi *= 0.5;
            }
        }

        if (pe < PE_tol)
            rescale(phi + d_phi);
        else
            rescale(phi - d_phi);

        last_pe = pe;
        pe = get_pe();
    }
}

void System::update () {
    N_c = 0;
    for (auto &p : particles) {
        p->reset_dynamics();
        p->update();
    }

    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        auto &p1 = particles[i];
        for (size_t j = i + 1; j < num_p; j++) if (active[j]) {
            auto &p2 = particles[j];
            p1->interact(p2.get());
        }
    }

    for (auto &p : particles) {
        if (relaxing) p->apply_drag(RELAX_KD);
        p->integrate();
    }

    time += dt;

    // check_verlet_rebuild();
    // if (rebuild_verlet_list) {
    //     build_verlet_list();
    // }
    //
    // N_c = 0;
    // for (auto &p : particles) {
    //     p->reset_dynamics();
    //     p->update();
    // }
    //
    // for (const auto& pair : verlet_list) {
    //     auto &p1 = particles[pair.first];
    //     auto &p2 = particles[pair.second];
    //     if (active[p1->id] && active[p2->id]) {
    //         p1->interact(p2.get());
    //     }
    // }
    //
    // for (size_t i = 0; i < num_p; i++) if (active[i]) {
    //     auto& p = particles[i];
    //     Eigen::Vector2d old_pos = p->com;
    //     if (relaxing) p->apply_drag(RELAX_KD);
    //     p->integrate();
    //     verlet_displacements[i] += p->com - old_pos;
    // }
    //
    // time += dt;
}

double System::get_ke () {
    double ke = 0.0;
    for (auto& p : particles)  {
        ke += p->get_ke();
    }
    return ke / num_p;
}

double System::get_pe () {
    double re = 0.0;
    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        auto &p1 = particles[i];
        for (size_t j = i + 1; j < num_p; j++) if (active[j]) {
            auto &p2 = particles[j];
            re += p1->get_energy_interaction(p2.get());
        }
    }

    return re / num_p;
}

void System::remove_rattlers () {
    for (size_t i = 0; i < num_p; i++) {
        if (particles[i]->rattles()) {
            active[i] = false;
        }
    }
}

size_t System::get_id () {
    curr_id.fetch_add(1);
    return curr_id.load() - 1;
}


std::vector<double> System::get_verts () const {
    // TODO: Fix me
    LOG_WARNING << "Using get_verts() is currently not advised because it has troubles when removing rattlers";
    std::vector<double> vertex_data(num_p * num_v * 2);

    for (size_t i = 0; i < num_p; i++) if (active[i]) {
        if (auto* smooth = dynamic_cast<SmoothParticle*>(particles[i].get())) {
            size_t idx = i * 2;
            vertex_data[idx + 0] = smooth->com[0];
            vertex_data[idx + 1] = smooth->com[1];
        }
        else if (auto* bumpy = dynamic_cast<BumpyParticle*>(particles[i].get())) {
            for (int v = 0; v < num_v; v++) {
                size_t idx = (i * num_v + v) * 2;
                vertex_data[idx + 0] = bumpy->verts[v][0];
                vertex_data[idx + 1] = bumpy->verts[v][1];
            }
        }
    }

    return vertex_data;
}
