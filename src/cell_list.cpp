#include "include/cell_list.h"
#include "include/system.h"

CellList::CellList (System* sys, double cutoff) :
    system(sys),
    cell_size(cutoff)
{
    nx = static_cast<int>(1.0 / cell_size);
    ny = static_cast<int>(1.0 / cell_size);
    if (nx == 0) nx = 1;
    if (ny == 0) ny = 1;
    cells.resize(nx * ny);
}

void CellList::build () {
    for (auto& cell : cells) {
        cell.clear();
    }
    for (size_t i = 0; i < system->num_p; i++) {
        if (system->active[i]) {
            int index = get_cell_index(system->particles[i]->com);
            if (index >= 0 && index < cells.size()) {
                cells[index].push_back(i);
            }
        }
    }
}

int CellList::get_cell_index (const Eigen::Vector2d& position) const {
    int ix = static_cast<int>(position.x() / cell_size);
    int iy = static_cast<int>(position.y() / cell_size);
    ix = std::max(0, std::min(ix, nx - 1));
    iy = std::max(0, std::min(iy, ny - 1));
    return iy * nx + ix;
}

void CellList::get_neighbor_cells (int particle_cell_index, std::vector<int>& neighbor_cells) const {
    neighbor_cells.clear();
    int ix_center = particle_cell_index % nx;
    int iy_center = particle_cell_index / nx;

    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            int ix = ix_center + dx;
            int iy = iy_center + dy;

            ix = (ix % nx + nx) % nx;
            iy = (iy % ny + ny) % ny;

            neighbor_cells.push_back(iy * nx + ix);
        }
    }
}
