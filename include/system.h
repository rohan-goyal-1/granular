#ifndef SYSTEM_H
#define SYSTEM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

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
    ~System() = default;

    void random_nonoverlap_config(void);
    double get_sigma(void);
    void update(void);
    double get_particle_area(void);
    double get_ke(void);
    double get_pe(void);
    void set_ke(double ke);
    double rescale(double _phi);
    void grow_to_phi(double _phi);
    inline bool is_relaxed(void);
    void relax(void);
    void fire_minimize(size_t max_steps = 10000, double dt_max = 1e-2, double alpha_start = 0.1, size_t n_min = 5, double force_tol = 1e-10);
    void send_to_jamming(void);
    size_t remove_rattlers(void);
    size_t purge_rattlers(void);
    std::vector<double> get_verts(void) const;
    void build_verlet_list();
    void check_verlet_rebuild();

    void apply_pbc_to_vector(Eigen::Vector2d& vec) const;
    Eigen::Vector2d get_min_dist_vec(const Particle* p1, const Particle* p2) const;

    bool relaxing = false;

    size_t N_c = 0;
    Eigen::SparseMatrix<double> adj_contacts;
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
    size_t num_p;
    size_t active_p;
    const size_t num_v;
    const double E_i;
    std::vector<std::unique_ptr<Particle>> particles;
    std::vector<bool> active;
    const size_t MAX_ATTEMPTS = 1E4;
    const double i_d_phi = 1E-2;
    const double d_phi_min = 1E-8;
    static constexpr double INIT_PHI = 1E-4;
    const double KE_tol = 1E-20;
    const double PE_tol = 1E-10;
    const size_t MIN_STEPS = 100;
    const PARTICLE_TYPE particle_type;
    const double RELAX_KD = 10.0;

    double EPS = 1e-9;

    // --- Cell + neighbor list ---
    void rebuild_cell_list();
    void rebuild_neighbor_list();
    bool neighbor_list_invalid() const;

    std::vector<std::vector<size_t>> cell_bins;
    std::vector<std::vector<size_t>> neighbor_list;
    std::vector<Eigen::Vector2d> verlet_displacements;
    double r_cut = 0.0;
    double r_skin = 0.0;
    size_t num_cells_per_dim = 0;
    double cell_size = 0.0;

    size_t get_flat_cell_index(const Eigen::Vector2d& pos) const;
    std::vector<size_t> get_neighboring_cells(size_t flat_index) const;

private:

};

#endif // SYSTEM_H
