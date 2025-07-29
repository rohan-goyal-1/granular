#include "include/bumpy_particle.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/algorithms/union.hpp>
#include <boost/geometry/algorithms/area.hpp>

namespace bg = boost::geometry;
typedef bg::model::d2::point_xy<double> BoostPoint;
typedef bg::model::polygon<BoostPoint> BoostPolygon;
typedef bg::model::multi_polygon<BoostPolygon> BoostMultiPolygon;

BumpyParticle::BumpyParticle (System* system, size_t num_v, double mu_eff, size_t id) :
    Particle(system, calc_ratio(mu_eff, num_v), id),
    ang_v(0.0),
    ang_a(0.0),
    theta(0.0),
    num_v(num_v),
    mu_eff(mu_eff),
    torque(0.0)
{
    force = Eigen::Vector2d(0, 0);
    com = Eigen::Vector2d(0, 0);
    verts.reserve(num_v);
    for (size_t i = 0; i < num_v; i++) {
        verts.push_back(Eigen::Vector2d(0, 0));
    }
    randomize_position();
    set_random_ke(0.0);
    area = get_area();
}

double BumpyParticle::calc_ratio (double mu, size_t n) {
    return (std::sin(M_PI/n)) * std::sqrt(1/(mu * mu) + 1);
}

void BumpyParticle::set_random_ke (double ke) {
    double dir = system->dist_angle(system->gen);
    double trans_vel = sqrt((ke / 2) * 2 / get_mass());
    com_v = Eigen::Vector2d(trans_vel * cos(dir), trans_vel * sin(dir));
    com_a = Eigen::Vector2d(0, 0);

    double ang_vel = sqrt((ke / 2) * 2 / moi);
    short spin = system->dist_sign(system->gen) == 0 ? -1 : 1;
    ang_v = spin * ang_vel;
    ang_a = 0.0;
}

void BumpyParticle::set_ke (double ke) {
    double current_ke = get_ke();

    if (current_ke <= 1E-16) {
        set_random_ke(ke);
        return;
    }

    double scale = sqrt(ke / current_ke);

    com_v *= scale;
    ang_v *= scale;
}

double BumpyParticle::get_ke () {
    return 0.5 * num_v * com_v.squaredNorm() + 0.5 * moi * (ang_v * ang_v);
}

void BumpyParticle::rescale_ratio (double ratio) {
    auto old_com = com;
    com *= ratio;
    for (auto& c : verts) {
        c += com - old_com;
    }
}

void BumpyParticle::move (Eigen::Vector2d translation) {
    com += translation;
    for (auto& c : verts) {
        c += translation;
    }
}

void BumpyParticle::randomize_position () {
    for (size_t d = 0; d < 2; d++) {
        com[d] = system->dist_pos(system->gen);
    }
    theta = system->dist_angle(system->gen);
    generate_polygon();
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

    for (size_t i = 0; i < num_v; i++) {
        Eigen::Vector2d displacement(cos(i * angle + theta), sin(i * angle + theta));
        verts[i] = com + displacement;
    }

    moi = _moi();
}

double BumpyParticle::_moi () {
    return num_v;
    // double ret = 0.0;
    // for (auto& c : verts) {
    //     Eigen::Vector2d r = c - com;
    //     ret += r.squaredNorm();
    // }
    // return ret;
}

double BumpyParticle::get_area () {
    if (verts.empty()) return 0.0;

    BoostPolygon base;
    for (const auto& v : verts) {
        bg::append(base.outer(), BoostPoint(v.x(), v.y()));
    }
    bg::append(base.outer(), BoostPoint(verts[0].x(), verts[0].y()));
    bg::correct(base);

    std::vector<BoostPolygon> circles;
    bg::strategy::buffer::distance_symmetric<double> distance_strategy(sigma / 2.0);
    bg::strategy::buffer::join_round join_strategy(10000);
    bg::strategy::buffer::end_round end_strategy;
    bg::strategy::buffer::point_circle circle_strategy(10000);
    bg::strategy::buffer::side_straight side_strategy;

    for (const auto& v : verts) {
        BoostPoint center(v.x(), v.y());
        BoostMultiPolygon circle_result;
        bg::buffer(center, circle_result,
            distance_strategy, side_strategy, join_strategy, end_strategy, circle_strategy);

        for (const auto& p : circle_result) {
            circles.push_back(p);
        }
    }

    BoostMultiPolygon union_result;
    bg::union_(base, circles[0], union_result);

    for (size_t i = 1; i < circles.size(); ++i) {
        BoostMultiPolygon temp_result;
        bg::union_(union_result, circles[i], temp_result);
        union_result = std::move(temp_result);
    }

    double total_area = 0.0;
    for (const auto& poly : union_result) {
        total_area += bg::area(poly);
    }

    return total_area;
}

double BumpyParticle::get_mass () const {
    return num_v;
}

double BumpyParticle::get_energy_interaction (Particle* other) {
    if (auto* bumpy = dynamic_cast<BumpyParticle*>(other)) {
        double sigma_l = (sigma + bumpy->sigma) / 2;
        double pe = 0.0;
        for (auto &c1 : verts) {
            for (auto &c2 : bumpy->verts) {
                Eigen::Vector2d r = minimum_image(c2 - c1, system->L);

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
    size_t count = 0;
    for (size_t i = 0; i < num_v; i++) {
        for (size_t j = 0; j < system->num_p * num_v; j++) {
            count += system->adj_contacts.coeffRef(num_v * id + i, j);
        }
    }
    return count < 3;
}

void BumpyParticle::interact (Particle* other, std::vector<Eigen::Triplet<double>>& new_contacts) {
    if (auto* bumpy = dynamic_cast<BumpyParticle*>(other)) {
        double sigma_l = (sigma + bumpy->sigma) / 2;
        for (size_t i = 0; i < num_v; i++) {
            auto& c1 = verts[i];
            for (size_t j = 0; j < bumpy->num_v; j++) {
                auto& c2 = bumpy->verts[j];

                Eigen::Vector2d r = minimum_image(c2 - c1, system->L);

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

                size_t row = num_v * id + i;
                size_t col = bumpy->num_v * bumpy->id + j;
                if (l_force.norm() > 0.0) {
                    new_contacts.emplace_back(row, col, 1.0);
                    new_contacts.emplace_back(col, row, 1.0);
                }
                else {
                    new_contacts.emplace_back(row, col, 0.0);
                    new_contacts.emplace_back(col, row, 0.0);
                }
            }
        }
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
    Eigen::Vector2d n_com_a = force / get_mass();
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
