#ifndef SYSTEM_H
#define SYSTEM_H

#include <Eigen/Dense>

#include "include/particle.h"
#include "include/smooth_particle.h"
#include "include/bumpy_particle.h"

#include <vector>
#include <random>
#include <string>
#include <math.h>
#include <cstddef>
#include <memory>

class CellList;

enum class PARTICLE_TYPE {
    BUMPY,
    SMOOTH,
};

class BumpyParticle;
class SmoothParticle;

class System {

public:
    System(PARTICLE_TYPE particle_type, size_t num_p, size_t num_v, double E_i, double dt, double mu, double phi = INIT_PHI);
    ~System();

    double get_sigma(void);
    void update(void);
    double get_particle_area(void);
    double get_ke(void);
    double get_pe(void);
    void set_temp(double temp);
    void rescale(double _phi);
    void grow_to_phi(double _phi);
    inline bool is_relaxed(void);
    void relax(void);
    void fire_minimize(size_t max_steps = 10000, double dt_max = 1e-2, double alpha_start = 0.1, size_t n_min = 5, double force_tol = 1e-10);
    void send_to_jamming(void);
    void remove_rattlers(void);
    std::vector<double> get_verts(void) const;
    void build_verlet_list();
    void check_verlet_rebuild();

    void apply_pbc_to_vector(Eigen::Vector2d& vec) const;
    Eigen::Vector2d get_min_dist_vec(const Particle* p1, const Particle* p2) const;

    bool relaxing = false;

    size_t N_c = 0;
    Eigen::MatrixXd adj_contacts;
    size_t num_collisions = 0;
    double time = 0.0;
    double phi;
    double mu;
    double L;
    const double dt;
    std::mt19937 gen;
    std::uniform_real_distribution<> dist_angle;
    std::uniform_real_distribution<> dist_pos;
    std::uniform_int_distribution<> dist_sign;
    size_t id;
    static size_t get_id(void);
    static std::atomic<int> curr_id;
    const size_t num_p;
    const size_t num_v;
    const double E_i;
    std::vector<std::unique_ptr<Particle>> particles;
    std::vector<bool> active;
    const size_t MAX_ATTEMPTS = 1E3;
    const double i_d_phi = 1E-2;
    const double d_phi_min = 1E-10;
    static constexpr double INIT_PHI = 0.01;
    const double KE_tol = 1E-20;
    const double PE_tol = 1E-10;
    const size_t MIN_STEPS = 10;
    const PARTICLE_TYPE particle_type;
    const double RELAX_KD = 10.0;

    std::unique_ptr<CellList> cell_list;
    std::vector<std::pair<size_t, size_t>> verlet_list;
    std::vector<Eigen::Vector2d> verlet_displacements;
    double verlet_cutoff_sq;
    double verlet_skin_radius;
    bool rebuild_verlet_list;

    double EPS = 1e-9;

private:

};

#endif // SYSTEM_H
