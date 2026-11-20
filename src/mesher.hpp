#pragma once
// Closed triangle mesh of a TPMS lattice clipped to a box.

#include "tpms.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace tpms {

struct Vec3 {
    double x = 0, y = 0, z = 0;
    double& operator[](int i) { return i == 0 ? x : i == 1 ? y : z; }
    double operator[](int i) const { return i == 0 ? x : i == 1 ? y : z; }
};

inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(Vec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(Vec3 a) { return std::sqrt(dot(a, a)); }

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<std::array<uint32_t, 3>> triangles;  // counter-clockwise seen from outside
};

// Box [0, cells[a] * cell[a]] per axis, sampled with res grid steps per unit cell.
Mesh meshLattice(const Lattice& lat, const int cells[3], int res,
                 size_t* droppedDegenerate = nullptr);

struct ClosedCheck {
    size_t openEdges = 0;         // directed edges without a reverse partner
    size_t nonManifoldEdges = 0;  // directed edges used more than once
    bool closed() const { return openEdges == 0 && nonManifoldEdges == 0; }
};
ClosedCheck checkClosed(const Mesh& m);

// Connected pieces, each with its own vertex list.
std::vector<Mesh> splitComponents(const Mesh& m);

// Enclosed volume (positive for an outward-oriented closed mesh).
double volume(const Mesh& m);

}  // namespace tpms
