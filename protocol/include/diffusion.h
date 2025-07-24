#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "include/system.h"
#include "utils/include/logger.h"

#include <vector>
#include <utility>

class System;

std::vector<std::vector<std::pair<double, double>>> find_diffusion(System* sys, double start, double end, double time, size_t count, std::vector<double>& phis);

#endif // DIFFUSION_H
