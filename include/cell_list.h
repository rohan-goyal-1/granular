#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <Eigen/Dense>

#include <vector>

class System;

class CellList {
public:
    CellList(System* system, double cutoff);
    void build();
    void get_neighbor_cells(int particle_cell_index, std::vector<int>& neighbor_cells) const;
    const std::vector<size_t>& get_particles_in_cell(int cell_index) const { return cells[cell_index]; }
    int get_cell_index(const Eigen::Vector2d& position) const;

private:
    System* system;
    double cell_size;
    int nx, ny;
    std::vector<std::vector<size_t>> cells;
};

#endif // CELL_LIST_H
