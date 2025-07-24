#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

#include <Eigen/Dense>

#include "include/system.h"

class System;

class Particle {

public:
    System* system;

    Eigen::Vector2d com_v;
    Eigen::Vector2d com_a;

    Eigen::Vector2d com;
    double theta;

    double sigma;
    double radius;

    Eigen::Vector2d force;

    size_t id;

    Particle(System* system, double sigma, size_t id)
        : system(system), sigma(sigma), radius(sigma), id(id) {}
    virtual ~Particle() = default;

    virtual void rescale(double area) = 0;

    virtual void update(void) = 0;
    virtual void move(Eigen::Vector2d translation) = 0;
    virtual void rotate(double angle) = 0;
    virtual void apply_drag(double kd) = 0;

    virtual void set_ke(double ke) = 0;
    virtual double get_ke(void) = 0;

    virtual double get_energy_interaction(Particle* other) = 0;
    virtual void interact(Particle* other) = 0;
    virtual void integrate (void) = 0;
    virtual void reset_dynamics(void) = 0;

    virtual bool rattles(void) = 0;

    inline Eigen::Vector2d minimum_image (const Eigen::Vector2d& vec) {
        return vec - vec.array().round().matrix();
    }

    inline double heaviside (double x) {
        return (x < 0.0) ? 0.0 : 1.0;
    }

private:

};

#endif // PARTICLE_H
