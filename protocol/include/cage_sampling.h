#ifndef CAGE_SAMPLING_H
#define CAGE_SAMPLING_H

#include "include/system.h"

#include <numeric>
#include <set>

#include <Eigen/Dense>

// === DATA STRUCTURES ===
struct Triangle {
    Eigen::Vector2d a, b, c;
    double area;
};

// --- HELPER FOR CONVEX HULL ---
// Computes the Z-component of the cross product of vectors (p1-p0) and (p2-p0).
// > 0 for counter-clockwise turn, < 0 for clockwise, = 0 for collinear.
double cross_product_z (const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    return (p1.x() - p0.x()) * (p2.y() - p0.y()) - (p1.y() - p0.y()) * (p2.x() - p0.x());
}

// === GEOMETRY UTILS ===
std::vector<Triangle> triangulate (const std::vector<Eigen::Vector2d>& poly) {
    std::vector<Triangle> tris;
    for (size_t i = 1; i + 1 < poly.size(); i++) {
        Triangle tri{poly[0], poly[i], poly[i + 1]};
        tri.area = 0.5 * std::abs((tri.b - tri.a).x() * (tri.c - tri.a).y() -
                                  (tri.b - tri.a).y() * (tri.c - tri.a).x());
        tris.push_back(tri);
    }
    return tris;
}

Eigen::Vector2d sample_point (const Triangle& tri, std::mt19937& gen, std::uniform_real_distribution<>& dist) {
    double u = dist(gen), v = dist(gen);
    if (u + v > 1.0) { u = 1.0 - u; v = 1.0 - v; }
    return tri.a + u * (tri.b - tri.a) + v * (tri.c - tri.a);
}

struct CageResult {
    std::vector<Eigen::Vector2d> cage_polygon;
    std::vector<size_t> neighbor_indices;
};

std::vector<Eigen::Vector2d> generate_cage (const std::vector<Eigen::Vector2d>& particles, size_t center_idx) {
    const std::vector<Eigen::Vector2d> shifts = {
        {-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 0}, {0, 1}, {1, -1}, {1, 0}, {1, 1}
    };

    const Eigen::Vector2d& center = particles[center_idx];
    std::vector<Eigen::Vector2d> cage = {
        center + Eigen::Vector2d(-1.0, -1.0),
        center + Eigen::Vector2d(1.0, -1.0),
        center + Eigen::Vector2d(1.0, 1.0),
        center + Eigen::Vector2d(-1.0, 1.0)
    };

    for (size_t pi = 0; pi < particles.size(); ++pi) {
        if (pi == center_idx) continue;
        for (const auto& shift : shifts) {
            Eigen::Vector2d img = particles[pi] + shift;
            Eigen::Vector2d dr = center - img;
            Eigen::Vector2d midpoint = center - 0.5 * dr;
            Eigen::Vector2d normal = dr;

            std::vector<Eigen::Vector2d> new_cage;
            if (cage.empty()) return {};

            Eigen::Vector2d prev = cage.back();
            bool prev_in = normal.dot(prev - midpoint) >= 0;
            for (const auto& curr : cage) {
                bool curr_in = normal.dot(curr - midpoint) >= 0;
                if (prev_in != curr_in) {
                    Eigen::Vector2d edge = curr - prev;
                    double t = normal.dot(midpoint - prev) / normal.dot(edge);
                    new_cage.push_back(prev + t * edge);
                }
                if (curr_in) new_cage.push_back(curr);
                prev = curr;
                prev_in = curr_in;
            }
            cage = std::move(new_cage);
        }
    }
    return cage;
}

// --- CONVEX HULL ALGORITHM (Monotone Chain) ---
std::vector<Eigen::Vector2d> compute_convex_hull(std::vector<Eigen::Vector2d>& points) {
    if (points.size() <= 3) return points;

    // Sort points lexicographically
    std::sort(points.begin(), points.end(), [](const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
        return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y());
    });

    std::vector<Eigen::Vector2d> upper_hull, lower_hull;
    for (const auto& p : points) {
        while (lower_hull.size() >= 2 && cross_product_z(lower_hull[lower_hull.size()-2], lower_hull.back(), p) <= 0) {
            lower_hull.pop_back();
        }
        lower_hull.push_back(p);
    }

    for (int i = points.size() - 1; i >= 0; --i) {
        const auto& p = points[i];
        while (upper_hull.size() >= 2 && cross_product_z(upper_hull[upper_hull.size()-2], upper_hull.back(), p) <= 0) {
            upper_hull.pop_back();
        }
        upper_hull.push_back(p);
    }

    lower_hull.pop_back();
    upper_hull.pop_back();
    lower_hull.insert(lower_hull.end(), upper_hull.begin(), upper_hull.end());
    return lower_hull;
}

// --- THE FULL OPTIMIZED PROTOCOL ---
// This function replaces the old generate_cage call in main.cpp
std::vector<Eigen::Vector2d> generate_minkowski_cage(
    size_t center_local_idx,
    const std::vector<Eigen::Vector2d>& active_particle_COMs,
    const std::unordered_map<size_t, size_t>& local_to_global_map,
    const std::vector<std::unique_ptr<Particle>>& all_particles
) {
    std::vector<Eigen::Vector2d> cage = generate_cage(active_particle_COMs, center_local_idx);
    if (cage.empty()) return {};

    size_t center_global_idx = local_to_global_map.at(center_local_idx);
    double center_radius = all_particles[center_global_idx]->sigma;

    const int points_per_circle = 32;
    for (size_t neighbor_local_idx = 0; neighbor_local_idx < active_particle_COMs.size(); ++neighbor_local_idx) {
        if (neighbor_local_idx == center_local_idx) continue;

        size_t neighbor_global_idx = local_to_global_map.at(neighbor_local_idx);
        const auto& neighbor = all_particles[neighbor_global_idx];
        double forbidden_radius = center_radius + neighbor->sigma;

        std::vector<Eigen::Vector2d> circle_pts;
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                Eigen::Vector2d neighbor_pos = neighbor->com + Eigen::Vector2d(dx, dy);
                for (int k = 0; k < points_per_circle; ++k) {
                    double angle = 2.0 * M_PI * k / points_per_circle;
                    circle_pts.push_back(neighbor_pos + forbidden_radius * Eigen::Vector2d(std::cos(angle), std::sin(angle)));
                }
            }
        }

        std::vector<Eigen::Vector2d> hull = compute_convex_hull(circle_pts);
        if (hull.size() < 3) continue;

        for (size_t i = 0; i < hull.size(); ++i) {
            Eigen::Vector2d p1 = hull[i];
            Eigen::Vector2d p2 = hull[(i + 1) % hull.size()];
            Eigen::Vector2d edge = p2 - p1;
            Eigen::Vector2d normal(-edge.y(), edge.x());

            std::vector<Eigen::Vector2d> new_cage;
            Eigen::Vector2d prev = cage.back();
            bool prev_in = normal.dot(prev - p1) >= 0;

            for (const auto& curr : cage) {
                bool curr_in = normal.dot(curr - p1) >= 0;
                if (prev_in != curr_in) {
                    Eigen::Vector2d clip_edge = curr - prev;
                    double denominator = normal.dot(clip_edge);
                    if (std::abs(denominator) > 1e-9) {
                        double t = normal.dot(p1 - prev) / denominator;
                        new_cage.push_back(prev + t * clip_edge);
                    }
                }
                if (curr_in) new_cage.push_back(curr);
                prev = curr;
                prev_in = curr_in;
            }

            cage = std::move(new_cage);
            if (cage.empty()) break;
        }
        if (cage.empty()) break;
    }

    return cage;
}
// === MAIN SAMPLING FUNCTION ===
std::vector<Eigen::Vector2d> generate_sample_points (const std::vector<Eigen::Vector2d>& cage, size_t N) {
    if (cage.size() < 3) {
        return {};
    }

    auto tris = triangulate(cage);
    if (tris.empty()) {
        return {};
    }

    std::vector<double> cumulative(tris.size());
    for (size_t i = 0; i < tris.size(); i++) cumulative[i] = tris[i].area;
    std::partial_sum(cumulative.begin(), cumulative.end(), cumulative.begin());
    double total_area = cumulative.back();

    std::vector<Eigen::Vector2d> samples;
    samples.reserve(N);

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    for (size_t i = 0; i < N; i++) {
        double r = dist(gen) * total_area;
        size_t idx = std::lower_bound(cumulative.begin(), cumulative.end(), r) - cumulative.begin();
        samples.push_back(sample_point(tris[idx], gen, dist));
    }
    return samples;
}

#endif // CAGE_SAMPLING_H
