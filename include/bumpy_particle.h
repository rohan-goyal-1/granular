#ifndef BUMPY_PARTICLE_H
#define BUMPY_PARTICLE_H

#include <vector>

#include <Eigen/Dense>

#include "include/system.h"
#include "include/particle.h"

class System;

class BumpyParticle : public Particle {

public:
    double ang_v;
    double ang_a;
    double theta;

    size_t num_v;
    double mu_eff;
    double area;
    std::vector<Eigen::Vector2d> verts;

    double torque;
    double moi;

    BumpyParticle(System* system, size_t num_v, double mu_eff, double area, size_t id);

    void rescale(double area) override;

    void update(void) override;
    void move(Eigen::Vector2d translation) override;
    void rotate(double angle) override;
    void apply_drag(double kd) override;

    void set_ke(double ke) override;
    double get_ke(void) override;
    double _moi(void);

    double get_energy_interaction(Particle* other) override;
    void interact(Particle* other) override;
    void integrate(void) override;
    void reset_dynamics(void) override;

    bool rattles(void) override;

    double calc_sigma(double mu, double target_area);
    double calc_ratio(double mu);

private:
    void generate_polygon(void);

};

#endif // BUMPY_PARTICLE_H
