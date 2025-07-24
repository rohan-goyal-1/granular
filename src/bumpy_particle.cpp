#include "include/bumpy_particle.h"

BumpyParticle::BumpyParticle (System* system, size_t num_v, double mu_eff, double area, size_t id) :
    Particle(system, calc_sigma(mu_eff, area), id),
    ang_v(0.0),
    ang_a(0.0),
    theta(0.0),
    num_v(num_v),
    mu_eff(mu_eff),
    area(area),
    torque(0.0)
{
    force = Eigen::Vector2d(0, 0);
    com = Eigen::Vector2d(0, 0);
    verts.reserve(num_v);
    for (size_t i = 0; i < num_v; i++) {
        verts.push_back(Eigen::Vector2d(0, 0));
    }
    for (size_t d = 0; d < 2; d++) {
        com[d] = system->dist_pos(system->gen);
    }
    theta = system->dist_angle(system->gen);

    generate_polygon();
    set_ke(0.0);
}

double BumpyParticle::calc_sigma (double mu, double target_area) {
    if (mu <= 0.0 || target_area <= 0.0) {
        throw std::invalid_argument("mu and target_area must be positive");
    }

    const double tol = 1e-8;
    const int max_iter = 1000;
    const double R_min = 1e-8;
    const double R_max = 1e8;

    double ratio = calc_ratio(mu);

    auto compute_area = [&] (double R) -> double {
        double Nv = static_cast<double>(num_v);
        double sigma_local = ratio * R;
        double r = sigma_local / 2.0;
        double theta = 2.0 * M_PI / Nv;
        double phi = M_PI - theta;
        double d = 2.0 * R * std::sin(M_PI / Nv);

        double A_polygon = 0.5 * Nv * R * R * std::sin(theta);

        double A_asperities = Nv * r * r * (M_PI - phi / 2.0);

        double A_overlap = 0.0;
        if (d < 2 * r) {
            double part1 = 2 * r * r * std::acos(d / (2 * r));
            double part2 = 0.5 * d * std::sqrt(4 * r * r - d * d);
            A_overlap = Nv * (part1 - part2);
        }

        return A_polygon + A_asperities - A_overlap;
    };

    double R_low = R_min;
    double R_high = R_max;
    double R_mid = 0.0;

    for (int iter = 0; iter < max_iter; ++iter) {
        R_mid = 0.5 * (R_low + R_high);
        double area = compute_area(R_mid);

        if (std::abs(area - target_area) < tol) {
            break;
        }

        if (area < target_area) {
            R_low = R_mid;
        }
        else {
            R_high = R_mid;
        }
    }

    return ratio * R_mid;
}

double BumpyParticle::calc_ratio (double mu) {
    return (std::sin(M_PI/num_v) / 2) * std::sqrt(1/(mu * mu) + 1);
}

void BumpyParticle::set_ke (double ke) {
    double dir = system->dist_angle(system->gen);
    double trans_vel = sqrt((ke / 2) * 2 / num_v);
    com_v = Eigen::Vector2d(trans_vel * cos(dir), trans_vel * sin(dir));
    com_a = Eigen::Vector2d(0, 0);

    double ang_vel = sqrt((ke / 2) * 2 / moi);
    short spin = system->dist_sign(system->gen) == 0 ? -1 : 1;
    ang_v = spin * ang_vel;
    ang_a = 0.0;
}


double BumpyParticle::get_ke () {
    return 0.5 * num_v * com_v.squaredNorm() + 0.5 * moi * (ang_v * ang_v);
}

void BumpyParticle::rescale (double a) {
    area = a;
    sigma = calc_sigma(mu_eff, a);

    generate_polygon();
}

void BumpyParticle::move (Eigen::Vector2d translation) {
    com += translation;
    for (auto& c : verts) {
        c += translation;
    }
}

void BumpyParticle::update () {
    Eigen::Vector2d d_trans = com_v * system->dt + com_a * 0.5 * system->dt * system->dt;
    move(d_trans);
    double d_theta = ang_v * system->dt + 0.5 * ang_a * system->dt * system->dt;
    rotate(d_theta);
}

void BumpyParticle::rotate (double angle) {
    double cos_theta = cos(angle);
    double sin_theta = sin(angle);

    for (auto& c : verts) {
        Eigen::Vector2d r = c - com;
        Eigen::Vector2d r_rotated(cos_theta * r[0] - sin_theta * r[1], sin_theta * r[0] + cos_theta * r[1]);
        c = com + r_rotated;
    }
    theta += angle;
}

void BumpyParticle::generate_polygon () {
    double angle = 2 * M_PI / num_v;
    radius = sigma / calc_ratio(mu_eff);

    for (size_t i = 0; i < num_v; i++) {
        Eigen::Vector2d displacement(radius * cos(i * angle + theta), radius * sin(i * angle + theta));
        verts[i] = com + displacement;
    }

    moi = _moi();
}

double BumpyParticle::_moi () {
    double ret = 0.0;
    for (auto& c : verts) {
        Eigen::Vector2d r = c - com;
        ret += r.squaredNorm();
    }
    return ret;
}

double BumpyParticle::get_energy_interaction (Particle* other) {
    if (auto* bumpy = dynamic_cast<BumpyParticle*>(other)) {
        double sigma_l = (sigma + bumpy->sigma) / 2;
        double pe = 0.0;
        for (auto &c1 : verts) {
            for (auto &c2 : bumpy->verts) {
                Eigen::Vector2d r = minimum_image(c2 - c1);

                double dist = r.norm();
                pe += 0.5 * pow((sigma_l - dist), 2) * heaviside(sigma_l - dist);
            }
        }
        return pe;
    }
    else {
        throw std::runtime_error("BumpyParticle INTERACTION: Didn't get another BumpyParticle to interact with");
    }
}

bool BumpyParticle::rattles () {
    bool ratt = true;
    for (size_t i = 0; i < num_v; i++) {
        for (size_t j = 0; j < system->num_p * num_v; j++) {
            ratt &= (system->adj_contacts(num_v * id + i, j) == 0);
        }
    }
    return ratt;
}

void BumpyParticle::interact (Particle* other) {
    // TODO: Apply some optimizations to collision handling for bumpy particles
    if (auto* bumpy = dynamic_cast<BumpyParticle*>(other)) {
        bool collides = false;
        double sigma_l = (sigma + bumpy->sigma) / 2;
        for (size_t i = 0; i < num_v; i++) {
            auto& c1 = verts[i];
            for (size_t j = 0; j < bumpy->num_v; j++) {
                auto& c2 = bumpy->verts[j];

                Eigen::Vector2d r = minimum_image(c2 - c1);

                double dist = r.norm();
                Eigen::Vector2d r_u;
                if (dist != 0.0) {
                    r_u = r / dist;
                }
                else {
                    r_u = Eigen::Vector2d(0, 0); // We pray this never happens
                }

                Eigen::Vector2d l_force = -r_u * (sigma_l - dist) * heaviside(sigma_l - dist);

                force += l_force;
                bumpy->force -= l_force;

                Eigen::Vector2d ra1 = c1 - com;
                Eigen::Vector2d ra2 = c2 - bumpy->com;

                torque += ra1(0) * l_force(1) - ra1(1) * l_force(0);
                bumpy->torque -= ra2(0) * l_force(1) - ra2(1) * l_force(0);

                if (l_force.norm() > 0.0) {
                    collides = true;
                    if (!system->adj_contacts(num_v * id + i, bumpy->num_v * bumpy->id + j)) system->num_collisions++;
                    system->adj_contacts(num_v * id + i, bumpy->num_v * bumpy->id + j) = 1;
                    system->adj_contacts(bumpy->num_v * bumpy->id + j, num_v * id + i) = 1;
                }
                else {
                    system->adj_contacts(num_v * id + i, bumpy->num_v * bumpy->id + j) = 0;
                    system->adj_contacts(bumpy->num_v * bumpy->id + j, num_v * id + i) = 0;
                }
            }
        }
        system->N_c += collides;
    }
    else {
        throw std::runtime_error("BumpyParticle INTERACTION: Didn't get another BumpyParticle to interact with");
    }
}

void BumpyParticle::apply_drag (double kd) {
    force += com_v * -kd;
    torque += ang_v * -kd * (sigma * sigma);
}

void BumpyParticle::integrate () {
    Eigen::Vector2d n_com_a = force / num_v;
    com_v += (com_a + n_com_a) * (system->dt / 2.0);
    com_a = n_com_a;

    double n_ang_a = torque / moi;
    ang_v += (ang_a + n_ang_a) * (system->dt / 2.0);
    ang_a = n_ang_a;
}

void BumpyParticle::reset_dynamics () {
    force.setZero();
    torque = 0.0;
}
