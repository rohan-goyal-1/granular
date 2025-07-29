#ifndef CAGE_SAMPLING_H
#define CAGE_SAMPLING_H

#include "include/system.h"

#include <numeric>
#include <set>

#include <Eigen/Dense>

using Vector2d = Eigen::Vector2d;
using Polygon = std::vector<Vector2d>;

struct Line {
    Vector2d normal;
    double offset; // n · x = offset
};

// Minimum image distance (wrapped vector)
Vector2d minimum_image (const Vector2d& vec, double L) {
    return vec - (vec / L).array().round().matrix() * L;
}

// Bisector function
Line bisector (const Vector2d& a, const Vector2d& b) {
    Vector2d mid = 0.5 * (a + b);
    Vector2d normal = b - a;

    if (normal.squaredNorm() < 1e-12) {
        return {Vector2d(0,0), 0};
    }

    normal.normalize();
    double offset = normal.dot(mid);
    return {normal, offset};
}
// Clip polygon with a half-plane defined by line (n · x <= offset)
Polygon clip_polygon (const Polygon& poly, const Line& line) {
    Polygon new_poly;
    size_t n = poly.size();

    for (size_t i = 0; i < n; ++i) {
        const Vector2d& A = poly[i];
        const Vector2d& B = poly[(i + 1) % n];

        double da = line.normal.dot(A) - line.offset;
        double db = line.normal.dot(B) - line.offset;

        if (da <= 0) new_poly.push_back(A); // A is inside

        if (da * db < 0) { // crosses boundary
            double t = da / (da - db);
            Vector2d intersect = A + t * (B - A);
            new_poly.push_back(intersect);
        }
    }

    return new_poly;
}

// Main cage function
Polygon make_cage (size_t center, const std::vector<Vector2d>& points, double L) {
    const double INF = 1e9;
    Polygon cage = {
        Vector2d(-INF, -INF),
        Vector2d(INF, -INF),
        Vector2d(INF, INF),
        Vector2d(-INF, INF)
    };

    size_t n = points.size();
    for (size_t i = 0; i < n; i++) if (center != i) {
        Vector2d true_neighbor = points[center] + minimum_image(points[i] - points[center], L);
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                Vector2d shift(dx * L, dy * L);
                Vector2d shifted = true_neighbor + shift;
                Line bis = bisector(points[center], shifted);
                cage = clip_polygon(cage, bis);
                if (cage.empty()) {
                    std::cout << "degenerate cage\n";
                    return {}; // degenerate cage
                }
            }
        }
    }

    return cage;
}

struct Triangle {
    Eigen::Vector2d a, b, c;
    double area;
};

std::vector<Triangle> triangulate (const Polygon& poly) {
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

std::vector<Eigen::Vector2d> generate_sample_points (const Polygon& cage, size_t N) {
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
