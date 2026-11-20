// Marching tetrahedra on a regular grid over the box, plus planar caps on the
// six box faces. Every cube is split into 6 tetrahedra along its main diagonal
// (Freudenthal split). The split is the same in every cube, so neighbor cubes
// and the caps share edges and vertices exactly. The result is watertight.

#include "mesher.hpp"

#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <utility>

namespace tpms {
namespace {

// Corner bitmask m of a cube: bit 0 = +x, bit 1 = +y, bit 2 = +z.
// A grid edge is keyed as (lower node, bitmask of its direction).
// Key mask 0 means the grid node itself.
constexpr double kEdgeMargin = 0.02;  // fraction of a grid edge

class Builder {
public:
    Builder(const Lattice& lat, const int cells[3], int res) {
        for (int a = 0; a < 3; ++a) {
            n_[a] = cells[a] * res;
            h_[a] = lat.cell[a] / res;
        }
        const double hmin = std::min({h_[0], h_[1], h_[2]});
        minArea2_ = 1e-12 * hmin * hmin;
        stride_[0] = 1;
        stride_[1] = n_[0] + 1;
        stride_[2] = int64_t(n_[0] + 1) * (n_[1] + 1);
        for (int m = 0; m < 8; ++m)
            off_[m] = (m & 1) * stride_[0] + ((m >> 1) & 1) * stride_[1] + ((m >> 2) & 1) * stride_[2];
        evaluate(lat);
    }

    Mesh run(size_t* dropped) {
        static const int kPerm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                                        {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
        for (int k = 0; k < n_[2]; ++k)
            for (int j = 0; j < n_[1]; ++j)
                for (int i = 0; i < n_[0]; ++i) {
                    const int64_t base = i + j * stride_[1] + k * stride_[2];
                    int nin = 0;
                    for (int m = 0; m < 8; ++m) nin += inside(base + off_[m]);
                    if (nin == 0 || nin == 8) continue;
                    for (const auto& p : kPerm) {
                        const int m[4] = {0, 1 << p[0], (1 << p[0]) | (1 << p[1]), 7};
                        tet(base, m);
                    }
                }
        for (int axis = 0; axis < 3; ++axis)
            for (int side = 0; side < 2; ++side) cap(axis, side);
        if (dropped) *dropped = dropped_;
        return std::move(mesh_);
    }

private:
    int n_[3];
    double h_[3];
    double minArea2_;
    int64_t stride_[3];
    int64_t off_[8];
    std::vector<float> g_;
    Mesh mesh_;
    std::unordered_map<uint64_t, uint32_t> vertexOf_;
    std::vector<uint32_t> poly_;
    size_t dropped_ = 0;

    void evaluate(const Lattice& lat) {
        g_.resize(size_t(stride_[2] * (n_[2] + 1)));
#pragma omp parallel for schedule(static)
        for (int k = 0; k <= n_[2]; ++k)
            for (int j = 0; j <= n_[1]; ++j)
                for (int i = 0; i <= n_[0]; ++i)
                    g_[size_t(i + j * stride_[1] + k * stride_[2])] =
                        float(lat(i * h_[0], j * h_[1], k * h_[2]));
    }

    bool inside(int64_t node) const { return g_[size_t(node)] < 0.0f; }

    Vec3 nodePos(int64_t node) const {
        const int64_t k = node / stride_[2], r = node % stride_[2];
        return {double(r % stride_[1]) * h_[0], double(r / stride_[1]) * h_[1], double(k) * h_[2]};
    }

    Vec3 cornerOffset(int m) const {
        return {(m & 1) * h_[0], ((m >> 1) & 1) * h_[1], ((m >> 2) & 1) * h_[2]};
    }

    const Vec3& pos(uint32_t v) const { return mesh_.vertices[v]; }

    uint32_t addVertex(uint64_t key, Vec3 p) {
        const auto v = uint32_t(mesh_.vertices.size());
        mesh_.vertices.push_back(p);
        vertexOf_.emplace(key, v);
        return v;
    }

    uint32_t nodeVertex(int64_t node) {
        const uint64_t key = uint64_t(node) * 8;
        const auto it = vertexOf_.find(key);
        return it != vertexOf_.end() ? it->second : addVertex(key, nodePos(node));
    }

    // Surface vertex on the grid edge from `lower` along direction `mask`.
    // The crossing is kept at least kEdgeMargin from both grid nodes. This
    // prevents near-zero edges. The vertex stays on its own edge, so the
    // mesh topology does not change.
    uint32_t edgeVertex(int64_t lower, int mask) {
        const uint64_t key = uint64_t(lower) * 8 + uint64_t(mask);
        const auto it = vertexOf_.find(key);
        if (it != vertexOf_.end()) return it->second;
        const double ga = g_[size_t(lower)], gb = g_[size_t(lower + off_[mask])];
        const double t = std::clamp(ga / (ga - gb), kEdgeMargin, 1.0 - kEdgeMargin);
        return addVertex(key, nodePos(lower) + cornerOffset(mask) * t);
    }

    void tet(int64_t base, const int m[4]) {
        int64_t node[4];
        bool in[4];
        int nin = 0;
        for (int q = 0; q < 4; ++q) {
            node[q] = base + off_[m[q]];
            in[q] = inside(node[q]);
            nin += in[q];
        }
        if (nin == 0 || nin == 4) return;

        // Corners lie on a monotone path, so for p < q the mask of q contains the mask of p.
        auto cut = [&](int p, int q) {
            if (p > q) std::swap(p, q);
            return edgeVertex(node[p], m[q] & ~m[p]);
        };

        Vec3 cin, cout;
        for (int q = 0; q < 4; ++q) (in[q] ? cin : cout) = (in[q] ? cin : cout) + cornerOffset(m[q]);
        const Vec3 outward = cout * (1.0 / (4 - nin)) - cin * (1.0 / nin);

        poly_.clear();
        if (nin == 1 || nin == 3) {
            const bool loneIn = nin == 1;
            int lone = 0;
            while (in[lone] != loneIn) ++lone;
            for (int q = 0; q < 4; ++q)
                if (q != lone) poly_.push_back(cut(lone, q));
        } else {
            int i[2], o[2], ni = 0, no = 0;
            for (int q = 0; q < 4; ++q) (in[q] ? i[ni++] : o[no++]) = q;
            poly_ = {cut(i[0], o[0]), cut(i[0], o[1]), cut(i[1], o[1]), cut(i[1], o[0])};
        }
        polygon(outward);
    }

    // Cap on box face `axis` at min (side 0) or max (side 1).
    void cap(int axis, int side) {
        int u = (axis + 1) % 3, v = (axis + 2) % 3;
        if (u > v) std::swap(u, v);
        Vec3 outward;
        outward[axis] = side ? 1.0 : -1.0;
        const int mu = 1 << u, mv = 1 << v;
        // Same face diagonal (min corner to max corner) as the tetrahedra use.
        const int tris[2][3] = {{0, mu, mu | mv}, {0, mu | mv, mv}};
        int idx[3];
        idx[axis] = side ? n_[axis] : 0;
        for (idx[v] = 0; idx[v] < n_[v]; ++idx[v])
            for (idx[u] = 0; idx[u] < n_[u]; ++idx[u]) {
                const int64_t base = idx[0] + idx[1] * stride_[1] + idx[2] * stride_[2];
                for (const auto& t : tris) clipTriangle(base, t, outward);
            }
    }

    // Part of a face triangle that lies inside, as a polygon.
    void clipTriangle(int64_t base, const int m[3], const Vec3& outward) {
        int64_t node[3];
        bool in[3];
        int nin = 0;
        for (int q = 0; q < 3; ++q) {
            node[q] = base + off_[m[q]];
            in[q] = inside(node[q]);
            nin += in[q];
        }
        if (nin == 0) return;
        poly_.clear();
        for (int p = 0; p < 3; ++p) {
            const int q = (p + 1) % 3;
            if (in[p]) poly_.push_back(nodeVertex(node[p]));
            if (in[p] != in[q]) {
                const bool pLower = (m[p] & m[q]) == m[p];
                const int lo = pLower ? p : q, hi = pLower ? q : p;
                poly_.push_back(edgeVertex(node[lo], m[hi] & ~m[lo]));
            }
        }
        polygon(outward);
    }

    void polygon(const Vec3& outward) {
        const auto& p = poly_;
        if (p.size() == 4) {  // split along the shorter diagonal
            if (norm(pos(p[0]) - pos(p[2])) <= norm(pos(p[1]) - pos(p[3]))) {
                tri(p[0], p[1], p[2], outward);
                tri(p[0], p[2], p[3], outward);
            } else {
                tri(p[0], p[1], p[3], outward);
                tri(p[1], p[2], p[3], outward);
            }
            return;
        }
        for (size_t q = 1; q + 1 < p.size(); ++q) tri(p[0], p[q], p[q + 1], outward);
    }

    void tri(uint32_t a, uint32_t b, uint32_t c, const Vec3& outward) {
        const Vec3 n = cross(pos(b) - pos(a), pos(c) - pos(a));
        if (norm(n) <= minArea2_) {
            ++dropped_;
            return;
        }
        if (dot(n, outward) < 0) std::swap(b, c);
        mesh_.triangles.push_back({a, b, c});
    }
};

}  // namespace

Mesh meshLattice(const Lattice& lat, const int cells[3], int res, size_t* droppedDegenerate) {
    return Builder(lat, cells, res).run(droppedDegenerate);
}

ClosedCheck checkClosed(const Mesh& m) {
    std::unordered_map<uint64_t, uint32_t> count;
    count.reserve(m.triangles.size() * 3);
    for (const auto& t : m.triangles)
        for (int e = 0; e < 3; ++e) ++count[(uint64_t(t[e]) << 32) | t[(e + 1) % 3]];
    ClosedCheck c;
    for (const auto& [key, n] : count) {
        if (n > 1) ++c.nonManifoldEdges;
        if (count.find((key << 32) | (key >> 32)) == count.end()) ++c.openEdges;
    }
    return c;
}

std::vector<Mesh> splitComponents(const Mesh& m) {
    std::vector<uint32_t> parent(m.vertices.size());
    std::iota(parent.begin(), parent.end(), 0u);
    auto find = [&](uint32_t a) {
        while (parent[a] != a) a = parent[a] = parent[parent[a]];
        return a;
    };
    for (const auto& t : m.triangles) {
        parent[find(t[1])] = find(t[0]);
        parent[find(t[2])] = find(t[0]);
    }

    std::vector<Mesh> parts;
    std::unordered_map<uint32_t, size_t> partOf;
    std::vector<uint32_t> local(m.vertices.size(), UINT32_MAX);
    for (const auto& t : m.triangles) {
        const auto [it, isNew] = partOf.emplace(find(t[0]), parts.size());
        if (isNew) parts.emplace_back();
        Mesh& part = parts[it->second];
        std::array<uint32_t, 3> lt;
        for (int e = 0; e < 3; ++e) {
            uint32_t& id = local[t[e]];
            if (id == UINT32_MAX) {
                id = uint32_t(part.vertices.size());
                part.vertices.push_back(m.vertices[t[e]]);
            }
            lt[e] = id;
        }
        part.triangles.push_back(lt);
    }
    return parts;
}

double volume(const Mesh& m) {
    double v = 0;
    for (const auto& t : m.triangles)
        v += dot(m.vertices[t[0]], cross(m.vertices[t[1]], m.vertices[t[2]]));
    return v / 6.0;
}

}  // namespace tpms
