#ifndef SMOOTH_PARTICLE_H
#define SMOOTH_PARTICLE_H

#include <vector>

#include <Eigen/Dense>

#include "include/system.h"
#include "include/particle.h"

class System;

class SmoothParticle : public Particle {

public:
    SmoothParticle(System* system, double sigma, size_t id);

    void rescale(double area) override;

    void update(void) override;
    void move(Eigen::Vector2d translation) override;
    void rotate (double angle) override {

    }
    void apply_drag(double kd) override;

    void set_ke(double ke) override;
    double get_ke(void) override;

    double get_energy_interaction(Particle* other) override;
    void interact(Particle* other) override;
    void integrate(void) override;

    void reset_dynamics(void) override;

    bool rattles(void) override;

private:
};

#endif // SMOOTH_PARTICLE_H
