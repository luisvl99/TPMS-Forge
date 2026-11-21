#pragma once
// Closed triangle mesh of an implicit field in a box.

#include "tpms.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace tpms {

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<std::array<uint32_t, 3>> triangles;  // counter-clockwise seen from outside
};

// The box [origin, origin + size] sampled with steps[a] grid steps per axis.
// Axes set in `periodic` wrap around: the field matches at both ends, so
// those faces get no caps and the two ends are welded together.
struct Grid {
    double origin[3] = {0, 0, 0}, size[3] = {1, 1, 1};
    int steps[3] = {1, 1, 1};
    int periodic = 0;
};

// The region field < 0 inside the grid box, closed with planar caps on the
// non-periodic box faces.
using Field = std::function<double(double, double, double)>;
Mesh meshField(const Field& field, const Grid& grid, size_t* droppedDegenerate = nullptr);

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
