#include "include/smooth_particle.h"

SmoothParticle::SmoothParticle (System* system, size_t id) :
    Particle(system, 2.0, id)
{
    theta = 0.0;
    force = Eigen::Vector2d(0, 0);
    com = Eigen::Vector2d(0, 0);
    for (size_t d = 0; d < 2; d++) {
        com[d] = system->dist_pos(system->gen);
    }
    set_ke(0.0);
}

void SmoothParticle::rescale_ratio (double ratio) {
    com *= ratio;
}

double SmoothParticle::get_area () {
    return M_PI * sigma * sigma / 4;
}

void SmoothParticle::update (void) {
    Eigen::Vector2d d_trans = com_v * system->dt + com_a * 0.5 * system->dt * system->dt;
    move(d_trans);
}

void SmoothParticle::move (Eigen::Vector2d translation) {
    com += translation;
}

void SmoothParticle::randomize_position () {
    for (size_t d = 0; d < 2; d++) {
        com[d] = system->dist_pos(system->gen);
    }
}

void SmoothParticle::set_ke (double ke) {
    double dir = system->dist_angle(system->gen);
    double trans_vel = std::sqrt(ke * 2);
    com_v = Eigen::Vector2d(trans_vel * cos(dir), trans_vel * sin(dir));
    com_a = Eigen::Vector2d(0, 0);
}

double SmoothParticle::get_ke () {
    return 0.5 * com_v.squaredNorm();
}

double SmoothParticle::get_energy_interaction (Particle* other) {
    if (auto* smooth = dynamic_cast<SmoothParticle*>(other)) {
        double sigma_l = (sigma + smooth->sigma) / 2;
        Eigen::Vector2d r = minimum_image(smooth->com - com, system->L);
        double dist = r.norm();
        double pe = 0.5 * pow((sigma_l - dist), 2) * heaviside(sigma_l - dist);
        return pe;
    }
    else {
        throw std::runtime_error("SmoothParticle INTERACTION: Didn't get another SmoothParticle to interact with");
    }
}

void SmoothParticle::interact (Particle* other) {
    if (auto* smooth = dynamic_cast<SmoothParticle*>(other)) {
        double sigma_l = (sigma + smooth->sigma) / 2;

        Eigen::Vector2d r = minimum_image(smooth->com - com, system->L);
        double dist = r.norm();

        Eigen::Vector2d r_u;
        if (dist != 0.0) {
            r_u = r / dist;
        }
        else {
            r_u = Eigen::Vector2d(0, 0);
        }

        Eigen::Vector2d l_force = -r_u * (sigma_l - dist) * heaviside(sigma_l - dist);
        force += l_force;
        smooth->force -= l_force;

        // Update contact matrix
        if (l_force.norm() > 0.0) {
            system->N_c++;
            if (!system->adj_contacts(id, smooth->id)) system->num_collisions++;
            system->adj_contacts(id, smooth->id) = 1;
            system->adj_contacts(smooth->id, id) = 1;
        }
        else {
            system->adj_contacts(id, smooth->id) = 0;
            system->adj_contacts(smooth->id, id) = 0;
        }
    }
    else {
        throw std::runtime_error("SmoothParticle INTERACTION: Didn't get another SmoothParticle to interact with");
    }
}

bool SmoothParticle::rattles () {
    bool ratt = true;
    for (size_t i = 0; i < system->num_p; i++) {
        ratt &= (system->adj_contacts(id, i) == 0);
    }
    return ratt;
}

void SmoothParticle::apply_drag (double kd) {
    force += com_v * -kd;
}

void SmoothParticle::integrate () {
    Eigen::Vector2d n_com_a = force;
    com_v += (com_a + n_com_a) * (system->dt / 2.0);
    com_a = n_com_a;
}

void SmoothParticle::reset_dynamics () {
    force.setZero();
}
