// Smooth B-rep of a TPMS lattice in a box.
//
// Each wall (a level set f = c) gets a quad mesh by dual contouring on a coarse
// grid. Every quad becomes one bicubic B-spline face that interpolates a grid
// of points projected onto the exact surface. Neighbor faces share the same
// boundary samples, so their boundary curves are identical and the shell is
// watertight without sewing. Planar faces close the solid on the box faces.

#include "export.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <utility>

namespace tpms {
namespace {

constexpr int kDeg = 3;

int bitCount(int x) {
    int c = 0;
    for (; x; x &= x - 1) ++c;
    return c;
}

// Cubic interpolation of n+1 points at parameters i/n. Clamped knots by
// averaging (The NURBS Book, eq. 9.8).
struct Spline {
    int n;
    std::vector<double> knots;             // n + 5 values
    std::vector<std::vector<double>> inv;  // inverse collocation matrix
    std::vector<std::vector<double>> mid;   // basis values at (i + 0.5) / n
    std::vector<std::vector<double>> fine;  // basis values at i / (4 n), for volumes

    explicit Spline(int n_) : n(n_) {
        knots.assign(kDeg + 1, 0.0);
        for (int j = 1; j <= n - kDeg; ++j) knots.push_back((3.0 * j + 3.0) / (3.0 * n));
        knots.insert(knots.end(), kDeg + 1, 1.0);
        std::vector<std::vector<double>> a;
        for (int i = 0; i <= n; ++i) a.push_back(basis(double(i) / n));
        inv = invert(a);
        for (int i = 0; i < n; ++i) mid.push_back(basis((i + 0.5) / n));
        for (int i = 0; i <= 4 * n; ++i) fine.push_back(basis(i / (4.0 * n)));
    }

    std::vector<double> basis(double t) const {
        const size_t m = knots.size() - 1;
        std::vector<double> N(m, 0.0);
        if (t >= 1.0) N[size_t(n)] = 1.0;
        else
            for (size_t i = 0; i < m; ++i) N[i] = knots[i] <= t && t < knots[i + 1] ? 1.0 : 0.0;
        for (size_t p = 1; p <= kDeg; ++p)
            for (size_t i = 0; i + p < m; ++i) {
                const double d1 = knots[i + p] - knots[i], d2 = knots[i + p + 1] - knots[i + 1];
                N[i] = (d1 > 0 ? (t - knots[i]) / d1 * N[i] : 0.0) +
                       (d2 > 0 ? (knots[i + p + 1] - t) / d2 * N[i + 1] : 0.0);
            }
        N.resize(size_t(n) + 1);
        return N;
    }

    static std::vector<std::vector<double>> invert(std::vector<std::vector<double>> a) {
        const size_t n = a.size();
        std::vector<std::vector<double>> b(n, std::vector<double>(n, 0.0));
        for (size_t i = 0; i < n; ++i) b[i][i] = 1.0;
        for (size_t c = 0; c < n; ++c) {
            size_t p = c;
            for (size_t r = c + 1; r < n; ++r)
                if (std::fabs(a[r][c]) > std::fabs(a[p][c])) p = r;
            std::swap(a[c], a[p]);
            std::swap(b[c], b[p]);
            const double d = a[c][c];
            for (size_t k = 0; k < n; ++k) a[c][k] /= d, b[c][k] /= d;
            for (size_t r = 0; r < n; ++r) {
                const double f = a[r][c];
                if (r == c || f == 0) continue;
                for (size_t k = 0; k < n; ++k) a[r][k] -= f * a[c][k], b[r][k] -= f * b[c][k];
            }
        }
        return b;
    }

    std::vector<Vec3> fit(const std::vector<Vec3>& pts) const {
        std::vector<Vec3> c(size_t(n) + 1);
        for (int i = 0; i <= n; ++i)
            for (int k = 0; k <= n; ++k) c[size_t(i)] = c[size_t(i)] + pts[size_t(k)] * inv[size_t(i)][size_t(k)];
        return c;
    }

    // STEP lists: knot multiplicities and distinct knots.
    std::string mults() const {
        std::string s = "(4";
        for (int i = 0; i < n - kDeg; ++i) s += ",1";
        return s + ",4)";
    }
    std::string distinct() const {
        std::string s = "(";
        for (size_t i = kDeg; i + kDeg < knots.size(); ++i) s += (i > kDeg ? "," : "") + real(knots[i]);
        return s + ")";
    }
};

struct Vert {
    Vec3 p;
    int faces;  // bitmask of box faces the vertex lies on, bit 2 * axis + side
};

struct Edge {
    int a, b;                // vertex ids, a < b; points run from a to b
    int wall;                // -1: straight segment on a box edge
    std::vector<Vec3> pts;   // samples on the surface, or the two end points
    std::vector<Vec3> ctrl;  // B-spline control points
};

struct Patch {
    int wall;
    std::array<int, 4> v;                       // corners, counter-clockwise seen from outside
    std::array<int, 4> e;                       // edge from v[k] to v[k + 1]
    std::vector<std::vector<Vec3>> pts, ctrl;   // [i][j], i along v0 -> v1, j along v0 -> v3
};

using Loop = std::vector<std::pair<int, bool>>;  // (edge, runs from a to b)

struct Face {
    int patch = -1;           // B-spline patch, or -1 for a planar cap
    int box = -1;             // box face of a cap
    std::vector<Loop> loops;  // first loop is the outer one
    double volume = 0;        // contribution to the enclosed volume
};

class Builder {
public:
    Builder(const Lattice& lat, const int cells[3], int res) : lat_(lat) {
        for (int a = 0; a < 3; ++a) {
            N_[a] = cells[a] * res;
            H_[a] = lat.cell[a] / res;
            L_[a] = cells[a] * lat.cell[a];
        }
        hmin_ = std::min({H_[0], H_[1], H_[2]});
        sx_ = N_[0] + 1;
        sxy_ = int64_t(sx_) * (N_[1] + 1);
        // phi = sign * (f - level) equals the lattice field near each wall.
        const double s = lat.part == Part::Solid ? 1.0 : -1.0;
        walls_.push_back({lat.level, s});
        if (lat.form == Form::Sheet) walls_.push_back({-lat.level, -s});
        for (size_t w = 0; w < walls_.size(); ++w) {
            phi_.emplace_back(size_t(sxy_ * (N_[2] + 1)));
            int i[3];
            for (i[2] = 0; i[2] <= N_[2]; ++i[2])
                for (i[1] = 0; i[1] <= N_[1]; ++i[1])
                    for (i[0] = 0; i[0] <= N_[0]; ++i[0]) phi_[w][size_t(node(i))] = phi(int(w), nodePos(i));
        }
    }

    // Quads of all walls by dual contouring. Around each grid edge with a sign
    // change: the cubes, plus boundary squares and the edge itself at the box.
    void buildPatches() {
        static const int D[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        for (int w = 0; w < int(walls_.size()); ++w)
            for (int a = 0; a < 3; ++a) {
                const int u = (a + 1) % 3, v = (a + 2) % 3;
                int i[3];
                for (i[2] = 0; i[2] <= N_[2]; ++i[2])
                    for (i[1] = 0; i[1] <= N_[1]; ++i[1])
                        for (i[0] = 0; i[0] <= N_[0]; ++i[0]) {
                            if (i[a] == N_[a]) continue;
                            int j[3] = {i[0], i[1], i[2]};
                            ++j[a];
                            const double v0 = nodePhi(w, i), v1 = nodePhi(w, j);
                            if ((v0 < 0) == (v1 < 0)) continue;
                            bool valid[4];
                            int nValid = 0;
                            for (int q = 0; q < 4; ++q) {
                                const int cu = i[u] - D[q][0], cv = i[v] - D[q][1];
                                valid[q] = cu >= 0 && cu < N_[u] && cv >= 0 && cv < N_[v];
                                nValid += valid[q];
                            }
                            std::vector<int> poly;
                            for (int q = 0; q < 4; ++q) {
                                const int r = (q + 1) % 4;
                                if (valid[q]) {
                                    int lo[3] = {i[0], i[1], i[2]};
                                    lo[u] -= D[q][0];
                                    lo[v] -= D[q][1];
                                    poly.push_back(cellVertex(w, lo, 7, i, a));
                                }
                                if (valid[q] == valid[r]) continue;
                                int sq[3] = {i[0], i[1], i[2]}, mask;
                                if (D[q][1] == D[r][1]) {
                                    sq[v] -= D[q][1];
                                    mask = (1 << a) | (1 << v);
                                } else {
                                    sq[u] -= D[q][0];
                                    mask = (1 << a) | (1 << u);
                                }
                                poly.push_back(cellVertex(w, sq, mask, i, a));
                                if (valid[q] && nValid == 1) poly.push_back(cellVertex(w, i, 1 << a, i, a));
                            }
                            Vec3 out;
                            out[a] = v0 < 0 ? 1.0 : -1.0;
                            const Vec3 nrm = cross(verts_[size_t(poly[2])].p - verts_[size_t(poly[0])].p,
                                                   verts_[size_t(poly[3])].p - verts_[size_t(poly[1])].p);
                            if (dot(nrm, out) < 0) std::reverse(poly.begin(), poly.end());
                            Patch p;
                            p.wall = w;
                            for (int k = 0; k < 4; ++k) p.v[size_t(k)] = poly[size_t(k)];
                            for (int k = 0; k < 4; ++k) p.e[size_t(k)] = edge(p.v[size_t(k)], p.v[size_t((k + 1) % 4)], w);
                            patches_.push_back(p);
                        }
            }
    }

    // Samples, fits and caps for n samples per patch edge. Returns the max
    // fit error, or a negative value with `error` set.
    double run(int n, std::string& error) {
        if (offSurface_) {
            error = "a patch corner did not converge onto the surface";
            return -1;
        }
        spline_ = Spline(n);
        const Spline& sp = spline_;
        for (Edge& e : edges_)
            if (e.wall >= 0) sampleEdge(e, sp);
        faces_.clear();
        double maxErr = 0;
        for (size_t k = 0; k < patches_.size(); ++k) {
            maxErr = std::max(maxErr, fitPatch(patches_[k], sp));
            Face f;
            f.patch = int(k);
            Loop loop;
            for (int c = 0; c < 4; ++c)
                loop.push_back({patches_[k].e[size_t(c)], edges_[size_t(patches_[k].e[size_t(c)])].a == patches_[k].v[size_t(c)]});
            f.loops.push_back(loop);
            f.volume = patchVolume(patches_[k]);
            faces_.push_back(f);
        }
        for (int F = 0; F < 6; ++F)
            if (!buildCaps(F, error)) return -1;
        return maxErr;
    }

    // Connected closed shells: (faces, volume).
    bool shells(std::vector<std::pair<std::vector<int>, double>>& out, std::string& error) const {
        std::vector<int> parent(faces_.size());
        std::iota(parent.begin(), parent.end(), 0);
        auto find = [&](int x) {
            while (parent[size_t(x)] != x) x = parent[size_t(x)] = parent[size_t(parent[size_t(x)])];
            return x;
        };
        std::vector<std::vector<int>> use(edges_.size());
        for (size_t f = 0; f < faces_.size(); ++f)
            for (const Loop& l : faces_[f].loops)
                for (const auto& [e, same] : l) use[size_t(e)].push_back(int(f) * 2 + same);
        for (size_t e = 0; e < use.size(); ++e) {
            if (use[e].empty()) continue;
            if (use[e].size() != 2 || (use[e][0] & 1) == (use[e][1] & 1)) {
                error = "edge not shared by exactly two faces in opposite directions";
                return false;
            }
            parent[size_t(find(use[e][0] / 2))] = find(use[e][1] / 2);
        }
        std::unordered_map<int, size_t> shellOf;
        for (size_t f = 0; f < faces_.size(); ++f) {
            const auto [it, isNew] = shellOf.emplace(find(int(f)), out.size());
            if (isNew) out.push_back({});
            out[it->second].first.push_back(int(f));
            out[it->second].second += faces_[f].volume;
        }
        return true;
    }

    std::vector<long> write(StepFile& st, const std::vector<std::vector<int>>& shells) const {
        const std::string mults = spline_.mults(), knots = spline_.distinct();
        std::unordered_map<int, long> vertexId, edgeId;
        auto point = [&](Vec3 p) { return st.add("CARTESIAN_POINT(''," + triple(p) + ")"); };
        auto vertex = [&](int v) {
            auto it = vertexId.find(v);
            if (it != vertexId.end()) return it->second;
            const long id = st.add("VERTEX_POINT(''," + ref(point(verts_[size_t(v)].p)) + ")");
            vertexId.emplace(v, id);
            return id;
        };
        auto edgeCurve = [&](int e) {
            auto it = edgeId.find(e);
            if (it != edgeId.end()) return it->second;
            const Edge& E = edges_[size_t(e)];
            long curve;
            if (E.wall < 0) {
                const Vec3 d = E.pts[1] - E.pts[0];
                const long dir = st.add("DIRECTION(''," + triple(d * (1.0 / norm(d))) + ")");
                const long vec = st.add("VECTOR(''," + ref(dir) + "," + real(norm(d)) + ")");
                curve = st.add("LINE(''," + ref(point(E.pts[0])) + "," + ref(vec) + ")");
            } else {
                std::vector<long> cps;
                for (const Vec3& c : E.ctrl) cps.push_back(point(c));
                curve = st.add("B_SPLINE_CURVE_WITH_KNOTS('',3," + refList(cps) + ",.UNSPECIFIED.,.F.,.F.," +
                               mults + "," + knots + ",.UNSPECIFIED.)");
            }
            const long id = st.add("EDGE_CURVE(''," + ref(vertex(E.a)) + "," + ref(vertex(E.b)) + "," + ref(curve) + ",.T.)");
            edgeId.emplace(e, id);
            return id;
        };
        auto loop = [&](const Loop& l) {
            std::vector<long> oes;
            for (const auto& [e, same] : l)
                oes.push_back(st.add("ORIENTED_EDGE('',*,*," + ref(edgeCurve(e)) + (same ? ",.T.)" : ",.F.)")));
            return st.add("EDGE_LOOP(''," + refList(oes) + ")");
        };

        std::vector<long> solids;
        for (const auto& shell : shells) {
            std::vector<long> faces;
            for (int fi : shell) {
                const Face& f = faces_[size_t(fi)];
                long surface;
                if (f.patch >= 0) {
                    std::string rows = "(";
                    for (const auto& row : patches_[size_t(f.patch)].ctrl) {
                        std::vector<long> ids;
                        for (const Vec3& c : row) ids.push_back(point(c));
                        rows += (rows.size() > 1 ? ",\n" : "") + refList(ids);
                    }
                    surface = st.add("B_SPLINE_SURFACE_WITH_KNOTS('',3,3," + rows + "),.UNSPECIFIED.,.F.,.F.,.F.," +
                                     mults + "," + mults + "," + knots + "," + knots + ",.UNSPECIFIED.)");
                } else {
                    const int a = f.box / 2, side = f.box % 2, s = side ? (a + 1) % 3 : (a + 2) % 3;
                    Vec3 origin, nrm, ref_dir;
                    origin[a] = side ? L_[a] : 0.0;
                    nrm[a] = side ? 1.0 : -1.0;
                    ref_dir[s] = 1.0;
                    const long ax = st.add("AXIS2_PLACEMENT_3D(''," + ref(point(origin)) + "," +
                                           ref(st.add("DIRECTION(''," + triple(nrm) + ")")) + "," +
                                           ref(st.add("DIRECTION(''," + triple(ref_dir) + ")")) + ")");
                    surface = st.add("PLANE(''," + ref(ax) + ")");
                }
                std::vector<long> bounds;
                for (size_t k = 0; k < f.loops.size(); ++k)
                    bounds.push_back(st.add(std::string(k ? "FACE_BOUND" : "FACE_OUTER_BOUND") + "(''," +
                                            ref(loop(f.loops[k])) + ",.T.)"));
                faces.push_back(st.add("ADVANCED_FACE(''," + refList(bounds) + "," + ref(surface) + ",.T.)"));
            }
            solids.push_back(st.add("MANIFOLD_SOLID_BREP(''," + ref(st.add("CLOSED_SHELL(''," + refList(faces) + ")")) + ")"));
        }
        return solids;
    }

private:
    struct Wall {
        double level, sign;  // phi = sign * (f - level), inside where phi < 0
    };

    const Lattice& lat_;
    int N_[3];
    double H_[3], L_[3], hmin_;
    int sx_;
    int64_t sxy_;
    std::vector<Wall> walls_;
    std::vector<std::vector<double>> phi_;
    std::vector<Vert> verts_;
    std::unordered_map<uint64_t, int> vertOf_;
    std::vector<Edge> edges_;
    std::unordered_map<uint64_t, int> edgeOf_;
    std::vector<Patch> patches_;
    std::vector<Face> faces_;
    Spline spline_{4};
    bool offSurface_ = false;

    int64_t node(const int i[3]) const { return i[0] + sx_ * i[1] + sxy_ * i[2]; }
    double nodePhi(int w, const int i[3]) const { return phi_[size_t(w)][size_t(node(i))]; }
    Vec3 nodePos(const int i[3]) const { return {i[0] * H_[0], i[1] * H_[1], i[2] * H_[2]}; }

    double phi(int w, Vec3 p) const {
        return walls_[size_t(w)].sign * (lat_.value(p.x, p.y, p.z) - walls_[size_t(w)].level);
    }

    Vec3 grad(int w, Vec3 p) const {
        const double d = 1e-6 * hmin_;
        Vec3 g;
        for (int a = 0; a < 3; ++a) {
            Vec3 q = p, r = p;
            q[a] += d;
            r[a] -= d;
            g[a] = (phi(w, q) - phi(w, r)) / (2 * d);
        }
        return g;
    }

    // Newton steps onto phi = 0. Axes set in `fixed` keep their coordinate.
    Vec3 project(int w, Vec3 p, int fixed) const {
        const double maxStep = 0.25 * hmin_;
        for (int it = 0; it < 50; ++it) {
            const double v = phi(w, p);
            if (std::fabs(v) < 1e-13) break;
            Vec3 g = grad(w, p);
            for (int a = 0; a < 3; ++a)
                if (fixed >> a & 1) g[a] = 0;
            const double gg = dot(g, g);
            if (gg == 0) {  // critical point of f: step off it
                for (int a = 0; a < 3; ++a)
                    if (!(fixed >> a & 1)) p[a] += 1e-3 * hmin_;
                continue;
            }
            Vec3 step = g * (v / gg);
            if (norm(step) > maxStep) step = step * (maxStep / norm(step));
            p = p - step;
        }
        return p;
    }

    // Dual vertex of the surface piece through grid edge (e0, ea) in a grid
    // cell: cube (mask 7), boundary square (2 bits) or box edge (1 bit). A cell
    // can hold several pieces. Crossing edges are grouped through the cell's
    // faces; a face with 4 crossings is split by the asymptotic decider, which
    // gives the same answer in both cells that share the face. The vertex
    // starts at the centroid of the piece's crossings and is projected onto
    // the surface inside the cell's box faces.
    int cellVertex(int w, const int lo[3], int mask, const int e0[3], int ea) {
        struct Crossing {
            std::array<int, 3> p;
            int axis;
        };
        std::vector<Crossing> xs;
        for (int a = 0; a < 3; ++a) {
            if (!(mask >> a & 1)) continue;
            for (int s = 0; s < 8; ++s) {
                if ((s & ~mask) || (s >> a & 1)) continue;
                Crossing x{{lo[0] + (s & 1), lo[1] + (s >> 1 & 1), lo[2] + (s >> 2 & 1)}, a};
                std::array<int, 3> q = x.p;
                ++q[size_t(a)];
                if ((nodePhi(w, x.p.data()) < 0) != (nodePhi(w, q.data()) < 0)) xs.push_back(x);
            }
        }
        auto indexOf = [&](std::array<int, 3> p, int axis) {
            for (size_t k = 0; k < xs.size(); ++k)
                if (xs[k].p == p && xs[k].axis == axis) return int(k);
            return -1;
        };
        std::vector<int> parent(xs.size());
        std::iota(parent.begin(), parent.end(), 0);
        auto root = [&](int x) {
            while (parent[size_t(x)] != x) x = parent[size_t(x)];
            return x;
        };
        std::vector<std::pair<std::array<int, 3>, int>> cellFaces;  // (origin, normal axis)
        if (mask == 7) {
            for (int b = 0; b < 3; ++b)
                for (int o = 0; o < 2; ++o) {
                    std::array<int, 3> f0{lo[0], lo[1], lo[2]};
                    f0[size_t(b)] += o;
                    cellFaces.push_back({f0, b});
                }
        } else if (bitCount(mask) == 2) {
            int nb = 0;
            while (mask >> nb & 1) ++nb;
            cellFaces.push_back({{lo[0], lo[1], lo[2]}, nb});
        }
        for (const auto& [f0, nb] : cellFaces) {
            const int s1 = (nb + 1) % 3, s2 = (nb + 2) % 3;
            std::array<int, 3> c10 = f0, c01 = f0, c11 = f0;
            ++c10[size_t(s1)], ++c01[size_t(s2)], ++c11[size_t(s1)], ++c11[size_t(s2)];
            const int id[4] = {indexOf(f0, s1), indexOf(c01, s1), indexOf(f0, s2), indexOf(c10, s2)};
            std::vector<int> on;
            for (int k : id)
                if (k >= 0) on.push_back(k);
            if (on.size() == 2) parent[size_t(root(on[0]))] = root(on[1]);
            if (on.size() != 4) continue;
            const double A = nodePhi(w, f0.data()), B = nodePhi(w, c10.data()), C = nodePhi(w, c01.data()),
                         D = nodePhi(w, c11.data());
            const bool saddleIn = (A * D - B * C) / (A + D - B - C) < 0;
            const double val[4] = {A, B, C, D};      // corners 00, 10, 01, 11
            const int inc[4][2] = {{0, 2}, {0, 3}, {1, 2}, {1, 3}};  // their two face edges
            for (int k = 0; k < 4; ++k)
                if (saddleIn ? val[k] >= 0 : val[k] < 0) parent[size_t(root(id[inc[k][0]]))] = root(id[inc[k][1]]);
        }

        const int target = root(indexOf({e0[0], e0[1], e0[2]}, ea));
        int comp = 0;
        std::vector<int> seen;
        Vec3 sum;
        std::vector<Vec3> crossings;
        for (size_t k = 0; k < xs.size(); ++k) {
            const int r = root(int(k));
            if (std::find(seen.begin(), seen.end(), r) == seen.end()) {
                if (r == target) comp = int(seen.size());
                seen.push_back(r);
            }
            if (r != target) continue;
            std::array<int, 3> q = xs[k].p;
            ++q[size_t(xs[k].axis)];
            const double vp = nodePhi(w, xs[k].p.data()), vq = nodePhi(w, q.data());
            const Vec3 P = nodePos(xs[k].p.data());
            crossings.push_back(P + (nodePos(q.data()) - P) * (vp / (vp - vq)));
            sum = sum + crossings.back();
        }
        const uint64_t key = (uint64_t(w) << 56) | (uint64_t(node(lo)) << 5) | uint64_t(mask << 2) | uint64_t(comp);
        if (auto it = vertOf_.find(key); it != vertOf_.end()) return it->second;
        int faces = 0, fixed = 0;
        for (int b = 0; b < 3; ++b) {
            if (mask >> b & 1) continue;  // the other axes lie on the box boundary
            fixed |= 1 << b;
            faces |= 1 << (2 * b + (lo[b] == N_[b]));
        }
        // Start at the centroid. If Newton stalls there, start at a crossing:
        // crossings lie next to the surface.
        Vec3 p = project(w, sum * (1.0 / double(crossings.size())), fixed);
        for (size_t k = 0; k < crossings.size() && std::fabs(phi(w, p)) > 1e-9; ++k) p = project(w, crossings[k], fixed);
        offSurface_ |= std::fabs(phi(w, p)) > 1e-9;
        verts_.push_back({p, faces});
        vertOf_.emplace(key, int(verts_.size()) - 1);
        return int(verts_.size()) - 1;
    }

    int cornerVertex(int c) {  // c: bit b set = max side on axis b
        const uint64_t key = (uint64_t(15) << 56) | uint64_t(c);
        if (auto it = vertOf_.find(key); it != vertOf_.end()) return it->second;
        Vec3 p;
        int faces = 0;
        for (int b = 0; b < 3; ++b) {
            p[b] = (c >> b & 1) ? L_[b] : 0.0;
            faces |= 1 << (2 * b + (c >> b & 1));
        }
        verts_.push_back({p, faces});
        vertOf_.emplace(key, int(verts_.size()) - 1);
        return int(verts_.size()) - 1;
    }

    int edge(int a, int b, int wall) {
        const uint64_t key = (uint64_t(std::min(a, b)) << 32) | uint64_t(std::max(a, b));
        if (auto it = edgeOf_.find(key); it != edgeOf_.end()) return it->second;
        edges_.push_back({std::min(a, b), std::max(a, b), wall, {}, {}});
        edgeOf_.emplace(key, int(edges_.size()) - 1);
        return int(edges_.size()) - 1;
    }

    void sampleEdge(Edge& e, const Spline& sp) {
        const Vec3 A = verts_[size_t(e.a)].p, B = verts_[size_t(e.b)].p;
        const int common = verts_[size_t(e.a)].faces & verts_[size_t(e.b)].faces;
        int fixed = 0;
        for (int F = 0; F < 6; ++F)
            if (common >> F & 1) fixed |= 1 << (F / 2);
        e.pts.assign(size_t(sp.n) + 1, A);
        e.pts.back() = B;
        for (int k = 1; k < sp.n; ++k) e.pts[size_t(k)] = project(e.wall, A + (B - A) * (double(k) / sp.n), fixed);
        e.ctrl = sp.fit(e.pts);
    }

    std::vector<Vec3> side(const Patch& p, int k) const {  // samples from v[k] to v[k + 1]
        const Edge& e = edges_[size_t(p.e[size_t(k)])];
        std::vector<Vec3> s = e.pts;
        if (e.a != p.v[size_t(k)]) std::reverse(s.begin(), s.end());
        return s;
    }

    // Coons grid from the four boundary curves, projected onto the surface,
    // then interpolated. Returns the max distance at the cell mid points.
    double fitPatch(Patch& p, const Spline& sp) {
        const int n = sp.n;
        const auto B = side(p, 0), R = side(p, 1);
        auto T = side(p, 2), Lf = side(p, 3);
        std::reverse(T.begin(), T.end());    // v3 -> v2
        std::reverse(Lf.begin(), Lf.end());  // v0 -> v3
        const size_t m = size_t(n) + 1;
        p.pts.assign(m, std::vector<Vec3>(m));
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < m; ++j) {
                Vec3& P = p.pts[i][j];
                if (j == 0) P = B[i];
                else if (j == m - 1) P = T[i];
                else if (i == 0) P = Lf[j];
                else if (i == m - 1) P = R[j];
                else {
                    const double u = double(i) / n, v = double(j) / n;
                    P = B[i] * (1 - v) + T[i] * v + Lf[j] * (1 - u) + R[j] * u -
                        (B[0] * ((1 - u) * (1 - v)) + B[m - 1] * (u * (1 - v)) + T[m - 1] * (u * v) + T[0] * ((1 - u) * v));
                    P = project(p.wall, P, 0);
                }
            }
        std::vector<std::vector<Vec3>> tmp(m, std::vector<Vec3>(m));
        p.ctrl.assign(m, std::vector<Vec3>(m));
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < m; ++j)
                for (size_t k = 0; k < m; ++k) tmp[i][j] = tmp[i][j] + p.pts[k][j] * sp.inv[i][k];
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < m; ++j)
                for (size_t k = 0; k < m; ++k) p.ctrl[i][j] = p.ctrl[i][j] + tmp[i][k] * sp.inv[j][k];

        double err = 0;
        for (size_t i = 0; i < m - 1; ++i)
            for (size_t j = 0; j < m - 1; ++j) {
                Vec3 q;
                for (size_t k = 0; k < m; ++k)
                    for (size_t l = 0; l < m; ++l) q = q + p.ctrl[k][l] * (sp.mid[i][k] * sp.mid[j][l]);
                err = std::max(err, std::fabs(phi(p.wall, q)) / norm(grad(p.wall, q)));
            }
        return err;
    }

    // Points of an edge curve at the fine parameters (a line: its end points).
    std::vector<Vec3> finePoints(const Edge& e) const {
        if (e.wall < 0) return e.pts;
        std::vector<Vec3> out(spline_.fine.size());
        for (size_t i = 0; i < out.size(); ++i)
            for (size_t k = 0; k < e.ctrl.size(); ++k) out[i] = out[i] + e.ctrl[k] * spline_.fine[i][k];
        return out;
    }

    // Volume term of a patch, from the B-spline at the fine parameters.
    // Its boundary points equal finePoints() of its edges, so the cap
    // polygons close the same polyhedron.
    double patchVolume(const Patch& p) const {
        const auto& F = spline_.fine;
        const size_t f = F.size(), m = p.ctrl.size();
        std::vector<std::vector<Vec3>> tmp(f, std::vector<Vec3>(m)), s(f, std::vector<Vec3>(f));
        for (size_t i = 0; i < f; ++i)
            for (size_t l = 0; l < m; ++l)
                for (size_t k = 0; k < m; ++k) tmp[i][l] = tmp[i][l] + p.ctrl[k][l] * F[i][k];
        for (size_t i = 0; i < f; ++i)
            for (size_t j = 0; j < f; ++j)
                for (size_t l = 0; l < m; ++l) s[i][j] = s[i][j] + tmp[i][l] * F[j][l];
        double v = 0;
        for (size_t i = 0; i + 1 < s.size(); ++i)
            for (size_t j = 0; j + 1 < s.size(); ++j)
                v += dot(s[i][j], cross(s[i + 1][j], s[i + 1][j + 1])) +
                     dot(s[i][j], cross(s[i + 1][j + 1], s[i][j + 1]));
        return v / 6.0;
    }

    // Planar faces on box face F: the region inside the lattice, bounded by
    // the patch boundary curves (reversed) and inside parts of the perimeter.
    bool buildCaps(int F, std::string& error) {
        const int a = F / 2, side = F % 2;
        const int s = side ? (a + 1) % 3 : (a + 2) % 3, t = side ? (a + 2) % 3 : (a + 1) % 3;  // e_s x e_t = outward
        struct Directed {
            int from, to, edge;
        };
        std::vector<Directed> de;
        for (const Patch& p : patches_)
            for (size_t k = 0; k < 4; ++k) {
                const int x = p.v[k], y = p.v[(k + 1) % 4];
                if (verts_[size_t(x)].faces & verts_[size_t(y)].faces & (1 << F)) de.push_back({y, x, p.e[k]});
            }

        for (int c = 0; c < 8; ++c)
            if ((c >> a & 1) == side) cornerVertex(c);
        const double S = L_[s], T = L_[t], eps = 1e-9 * (S + T);
        std::vector<std::pair<double, int>> perim;  // counter-clockwise position, vertex
        for (size_t id = 0; id < verts_.size(); ++id) {
            if (!(verts_[id].faces >> F & 1) || bitCount(verts_[id].faces) < 2) continue;
            const double ps = verts_[id].p[s], pt = verts_[id].p[t];
            const double pos = pt < eps ? ps : ps > S - eps ? S + pt : pt > T - eps ? 2 * S + T - ps : 2 * S + 2 * T - pt;
            perim.push_back({pos, int(id)});
        }
        std::sort(perim.begin(), perim.end());
        for (size_t k = 0; k < perim.size(); ++k) {
            const int x = perim[k].second, y = perim[(k + 1) % perim.size()].second;
            const Vec3 m = (verts_[size_t(x)].p + verts_[size_t(y)].p) * 0.5;
            if (lat_(m.x, m.y, m.z) >= 0) continue;
            const int e = edge(x, y, -1);
            edges_[size_t(e)].pts = {verts_[size_t(edges_[size_t(e)].a)].p, verts_[size_t(edges_[size_t(e)].b)].p};
            de.push_back({x, y, e});
        }

        std::unordered_map<int, size_t> next;
        for (size_t k = 0; k < de.size(); ++k)
            if (!next.emplace(de[k].from, k).second) {
                error = "box face cap: vertex with two outgoing edges";
                return false;
            }
        std::vector<bool> used(de.size(), false);
        std::vector<Loop> loops;
        std::vector<std::vector<std::array<double, 2>>> polys;
        std::vector<double> areas;
        for (size_t k0 = 0; k0 < de.size(); ++k0) {
            if (used[k0]) continue;
            Loop loop;
            std::vector<std::array<double, 2>> poly;
            for (size_t k = k0; !used[k];) {
                used[k] = true;
                const Edge& e = edges_[size_t(de[k].edge)];
                const bool same = de[k].from == e.a;
                loop.push_back({de[k].edge, same});
                const std::vector<Vec3> pts = finePoints(e);
                for (size_t q = 0; q + 1 < pts.size(); ++q) {
                    const Vec3& P = pts[same ? q : pts.size() - 1 - q];
                    poly.push_back({P[s], P[t]});
                }
                const auto it = next.find(de[k].to);
                if (it == next.end()) {
                    error = "box face cap: open loop";
                    return false;
                }
                k = it->second;
            }
            double area = 0;
            for (size_t q = 0; q < poly.size(); ++q) {
                const auto& A = poly[q];
                const auto& B = poly[(q + 1) % poly.size()];
                area += A[0] * B[1] - B[0] * A[1];
            }
            loops.push_back(loop);
            polys.push_back(poly);
            areas.push_back(area / 2);
        }

        // Counter-clockwise loops are outer bounds. Each clockwise loop is a
        // hole in the smallest outer loop around it.
        auto inside = [&](const std::array<double, 2>& p, const std::vector<std::array<double, 2>>& poly) {
            bool in = false;
            for (size_t q = 0, r = poly.size() - 1; q < poly.size(); r = q++)
                if ((poly[q][1] > p[1]) != (poly[r][1] > p[1]) &&
                    p[0] < (poly[r][0] - poly[q][0]) * (p[1] - poly[q][1]) / (poly[r][1] - poly[q][1]) + poly[q][0])
                    in = !in;
            return in;
        };
        const double d = side ? L_[a] : 0.0;  // plane: outward normal . p = d
        std::unordered_map<size_t, size_t> faceOf;
        for (size_t k = 0; k < loops.size(); ++k) {
            if (areas[k] <= 0) continue;
            faceOf[k] = faces_.size();
            Face f;
            f.box = F;
            f.loops.push_back(loops[k]);
            f.volume = d * areas[k] / 3;
            faces_.push_back(f);
        }
        for (size_t k = 0; k < loops.size(); ++k) {
            if (areas[k] > 0) continue;
            size_t best = SIZE_MAX;
            for (const auto& [o, fi] : faceOf)
                if (inside(polys[k][0], polys[o]) && (best == SIZE_MAX || areas[o] < areas[best])) best = o;
            if (best == SIZE_MAX) {
                error = "box face cap: hole outside every outer loop";
                return false;
            }
            Face& f = faces_[faceOf[best]];
            f.loops.push_back(loops[k]);
            f.volume += d * areas[k] / 3;
        }
        return true;
    }
};

}  // namespace

SmoothReport writeSmoothStep(const std::string& path, const std::string& name, const Lattice& lat,
                             const int cells[3], int patches, double tol, double minFragment) {
    // A coarse patch grid can miss small features: a tunnel inside two cells,
    // or a corner that does not converge. The lattice repeats per cell, so a
    // patch count that is valid for one cell is valid for all. Step it up
    // until the topology is closed and the fit is within tolerance.
    // Folded patches can still close up, so also compare the enclosed volume
    // with a sample of the field. Quasi-random R3 points: a regular grid can
    // line up with the lattice symmetry and give a biased count.
    const double L[3] = {cells[0] * lat.cell[0], cells[1] * lat.cell[1], cells[2] * lat.cell[2]};
    const double boxVol = L[0] * L[1] * L[2], g = 1.2207440846057595;  // plastic number
    const double alpha[3] = {1 / g, 1 / (g * g), 1 / (g * g * g)};
    const int samples = 1 << 18;
    int in = 0;
    for (int k = 0; k < samples; ++k) {
        double x[3];
        for (int a = 0; a < 3; ++a) x[a] = std::fmod(0.5 + k * alpha[a], 1.0) * L[a];
        in += lat(x[0], x[1], x[2]) < 0;
    }
    const double sampled = boxVol * in / samples;

    SmoothReport r;
    for (int p = patches; p <= std::max(patches + 7, 16); ++p) {
        r = SmoothReport{};
        r.patches = p;
        Builder b(lat, cells, p);
        b.buildPatches();
        for (int n = 4; r.error.empty(); n *= 2) {
            r.maxError = b.run(n, r.error);
            r.samples = n;
            // A large error is folded geometry, which more samples do not fix.
            if (r.maxError <= tol || r.maxError > 20 * tol || n >= 32) break;
        }
        if (r.error.empty() && r.maxError > tol) r.error = "fit error above tolerance";
        std::vector<std::pair<std::vector<int>, double>> shells;
        if (!r.error.empty() || !b.shells(shells, r.error)) continue;
        double total = 0, maxVol = 0;
        for (const auto& s : shells) total += s.second, maxVol = std::max(maxVol, s.second);
        if (std::fabs(total - sampled) > 0.02 * boxVol) {
            r.error = "enclosed volume does not match the field (folded patches)";
            continue;
        }
        std::vector<std::vector<int>> kept;
        for (const auto& [faces, vol] : shells) {
            if (vol <= 0) ++r.cavities;
            else if (vol < minFragment * maxVol) ++r.fragments;
            else {
                kept.push_back(faces);
                r.volume += vol;
                r.faces += faces.size();
            }
        }
        r.solids = kept.size();
        if (kept.empty()) r.error = "nothing to export";
        else if (!writeStepFile(path, name, [&](StepFile& st) { return b.write(st, kept); }))
            r.error = "cannot write '" + path + "'";
        return r;
    }
    return r;
}

}  // namespace tpms
