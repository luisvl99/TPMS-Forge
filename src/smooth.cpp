// Smooth B-rep of a TPMS lattice in a box, or bent into a tube or torus shell.
//
// Each wall (a level set f = c) gets a quad mesh by dual contouring on a coarse
// grid. Every quad becomes one bicubic B-spline face that interpolates a grid
// of points projected onto the exact surface. Neighbor faces share the same
// boundary samples, so their boundary curves are identical and the shell is
// watertight without sewing. Caps close the solid on the box faces.
//
// Curved shapes run the same steps in their parameter box: samples are
// projected there and mapped to mm, where the B-splines are fitted. The
// faces where the parameter box wraps around (the seam) get no caps: their
// two copies are welded together at the end. Volumes are fluxes in the
// parameter box (see Shape::flux), so they stay exact for every shape.

#include "export.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <queue>
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
};

// STEP lists of a clamped knot vector: multiplicities and distinct values.
std::pair<std::string, std::string> knotLists(const std::vector<double>& U) {
    std::string m = "(", k = "(";
    for (size_t i = 0, j; i < U.size(); i = j) {
        for (j = i; j < U.size() && U[j] == U[i]; ++j) {
        }
        m += (i ? "," : "") + std::to_string(j - i);
        k += (i ? "," : "") + real(U[i]);
    }
    return {m + ")", k + ")"};
}

// Clamped cubic knot vector on [0, 1] with the given internal knots.
std::vector<double> clamped(const std::vector<double>& inner) {
    std::vector<double> U(4, 0.0);
    U.insert(U.end(), inner.begin(), inner.end());
    U.insert(U.end(), 4, 1.0);
    return U;
}

// The 4 nonzero cubic basis values at t (The NURBS Book, A2.2). Returns the
// index of the first one.
int cubicBasis(const std::vector<double>& U, double t, double N[4]) {
    const int k = int(U.size()) - 4;
    t = std::clamp(t, 0.0, 1.0);
    const int span = std::clamp(int(std::upper_bound(U.begin(), U.end(), t) - U.begin()) - 1, 3, k - 1);
    double left[4], right[4];
    N[0] = 1;
    for (int j = 1; j <= 3; ++j) {
        left[j] = t - U[size_t(span + 1 - j)];
        right[j] = U[size_t(span + j)] - t;
        double saved = 0;
        for (int r = 0; r < j; ++r) {
            const double tmp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * tmp;
            saved = left[j - r] * tmp;
        }
        N[j] = saved;
    }
    return span - 3;
}

// Solves A x = b (in place in b) for a symmetric positive definite band
// matrix: A[i * (w + 1) + (j - i)] holds A(i, j) for i <= j <= i + w.
bool bandSolve(std::vector<double> A, size_t n, size_t w, std::vector<double>& b) {
    auto at = [&](size_t i, size_t j) -> double& { return A[i * (w + 1) + (j - i)]; };
    // A = U^T U, row by row: each row updates the rows below it in its band,
    // so all accesses run along rows (several times faster than columns).
    for (size_t i = 0; i < n; ++i) {
        const size_t end = std::min(n, i + w + 1);
        if (at(i, i) <= 0) return false;
        const double d = at(i, i) = std::sqrt(at(i, i));
        for (size_t j = i + 1; j < end; ++j) at(i, j) /= d;
        for (size_t r = i + 1; r < end; ++r) {
            const double f = at(i, r);
            if (f == 0) continue;
            double* row = &at(r, r);
            const double* src = &at(i, r);
            for (size_t c = 0; c < end - r; ++c) row[c] -= f * src[c];
        }
    }
    for (size_t i = 0; i < n; ++i) {
        b[i] /= at(i, i);
        for (size_t j = i + 1; j < std::min(n, i + w + 1); ++j) b[j] -= at(i, j) * b[i];
    }
    for (size_t i = n; i-- > 0;) {
        for (size_t j = i + 1; j < std::min(n, i + w + 1); ++j) b[i] -= at(i, j) * b[j];
        b[i] /= at(i, i);
    }
    return true;
}

// Least-squares cubic B-spline with knots U through pts at parameters t,
// with fixed end points. False if the normal equations are singular.
bool lsqCurve(const std::vector<Vec3>& pts, const std::vector<double>& t, const std::vector<double>& U,
              std::vector<Vec3>& ctrl) {
    const size_t m = pts.size(), M = U.size() - 4, n = M - 2;  // n unknown interior control points
    std::vector<double> A(n * 4, 0.0), rhs[3];
    for (auto& r : rhs) r.assign(n, 0.0);
    for (size_t k = 0; k < m; ++k) {
        double N[4];
        const int i0 = cubicBasis(U, t[k], N);
        Vec3 P = pts[k];
        for (int p = 0; p < 4; ++p) {
            const size_t c = size_t(i0 + p);
            if (c == 0) P = P - pts.front() * N[p];
            if (c == M - 1) P = P - pts.back() * N[p];
        }
        for (int p = 0; p < 4; ++p) {
            const size_t c = size_t(i0 + p);
            if (c == 0 || c == M - 1) continue;
            for (int q = p; q < 4; ++q) {
                const size_t c2 = size_t(i0 + q);
                if (c2 != 0 && c2 != M - 1) A[(c - 1) * 4 + (c2 - c)] += N[p] * N[q];
            }
            for (int a = 0; a < 3; ++a) rhs[a][c - 1] += N[p] * P[a];
        }
    }
    for (auto& r : rhs)
        if (!bandSolve(A, n, 3, r)) return false;
    ctrl.assign(1, pts.front());
    for (size_t i = 0; i < n; ++i) ctrl.push_back({rhs[0][i], rhs[1][i], rhs[2][i]});
    ctrl.push_back(pts.back());
    return true;
}

// Least-squares cubic B-spline through pts with fixed end points: chord
// length parameters t, knots by The NURBS Book eq. 9.69, the fewest control
// points that keep every point within tol. Least squares stops at m / 2
// control points: with more, a knot span can hold a single sample at its
// end, and the curve swings far off between samples. Past that it
// interpolates (averaged knots, eq. 9.8), which is stable.
void fitCurve(std::vector<Vec3> pts, double tol, std::vector<Vec3>& ctrl, std::vector<double>& U) {
    // Repeated points (ends of near zero-length edges) would repeat a
    // parameter and make interpolation singular.
    double length = 0;
    for (size_t k = 1; k < pts.size(); ++k) length += norm(pts[k] - pts[k - 1]);
    std::vector<Vec3> kept(1, pts.front());
    for (size_t k = 1; k < pts.size(); ++k)
        if (norm(pts[k] - kept.back()) > 1e-9 * length) kept.push_back(pts[k]);
        else if (k + 1 == pts.size()) kept.back() = pts[k];
    pts.swap(kept);
    const size_t m = pts.size();
    std::vector<double> t(m, 0.0);
    for (size_t k = 1; k < m; ++k) t[k] = t[k - 1] + norm(pts[k] - pts[k - 1]);
    for (size_t k = 0; k < m; ++k) t[k] = t.back() > 0 ? t[k] / t.back() : double(k) / double(m - 1);
    for (size_t M = std::min<size_t>(m, 4);;) {
        std::vector<double> inner;
        const double d = double(m) / double(M - 3);
        for (size_t j = 1; j + 3 < M; ++j) {
            if (M == m) {
                inner.push_back((t[j] + t[j + 1] + t[j + 2]) / 3);
                continue;
            }
            const size_t i = size_t(j * d);
            const double a = j * d - double(i);
            inner.push_back((1 - a) * t[i - 1] + a * t[i]);
        }
        U = clamped(inner);
        double err = 0;
        const bool solved = lsqCurve(pts, t, U, ctrl);
        for (size_t k = 0; k < m && solved; ++k) {
            double N[4];
            const int i0 = cubicBasis(U, t[k], N);
            Vec3 q;
            for (int p = 0; p < 4; ++p) q = q + ctrl[size_t(i0 + p)] * N[p];
            err = std::max(err, norm(q - pts[k]));
        }
        if ((solved && err <= tol) || M == m) return;
        M = M + M / 2 + 1 > m / 2 ? m : M + M / 2 + 1;
    }
}

struct Vert {
    Vec3 p;        // in the parameter box
    int faces;     // bitmask of box faces the vertex lies on, bit 2 * axis + side
    uint64_t key;  // in Builder::vertOf_
};

// Points are kept in the parameter box, where the lattice is sampled, and in
// mm, where the B-splines are fitted. For the box both are the same.
struct Edge {
    int a, b;                   // vertex ids; points run from a to b
    int wall;                   // -1: on a box edge (straight in the parameter box)
    std::vector<Vec3> pts, mm;  // samples in the parameter box and in mm
    std::vector<Vec3> ctrl;     // B-spline control points in mm; empty: a line
    std::vector<double> knots;  // empty: the patch spline's knots
    std::vector<Vec3> pctrl;    // wall edge: control points in the parameter box
};

struct Patch {
    int wall;
    std::array<int, 4> v;  // corners, counter-clockwise seen from outside
    std::array<int, 4> e;  // edge from v[k] to v[k + 1]
    // [i][j], i along v0 -> v1, j along v0 -> v3: samples and control points
    // in the parameter box (pts, pctrl) and in mm (mm, ctrl).
    std::vector<std::vector<Vec3>> pts, pctrl, mm, ctrl;
    std::vector<Vec3> mids;  // the exact surface next to the spline's cell mid points, mm
};

struct Directed {
    int from, to, edge;
    bool same;  // runs from the edge's a to b (a closed edge has a == b)
};

using Loop = std::vector<std::pair<int, bool>>;  // (edge, runs from a to b)

struct Face {
    int patch = -1;                          // B-spline patch
    int box = -1;                            // box face of a planar cap
    std::vector<std::vector<Vec3>> ctrl;     // trimmed B-spline of a merged region
    std::vector<double> knots[2];            // its knots along u and v
    std::vector<Loop> loops;                 // first loop is the outer one
    double volume = 0;                       // contribution to the enclosed volume
};

class Builder {
public:
    Builder(const Lattice& lat, const int cells[3], int res, double tol)
        : lat_(lat), shape_(lat.shape), mapped_(lat.shape.mapped()), periodic_(lat.shape.periodic()),
          tol_(tol) {
        for (int a = 0; a < 3; ++a) {
            N_[a] = cells[a] * res;
            H_[a] = lat.cell[a] / res;
            L_[a] = cells[a] * lat.cell[a];
        }
        hmin_ = std::min({H_[0], H_[1], H_[2]});
        sx_ = N_[0] + 1;
        sxy_ = int64_t(sx_) * (N_[1] + 1);
        // phi = sign * (f - level) equals the lattice field near each wall:
        // f = c, and for a sheet f = -c.
        const double s = lat.part == Part::Solid ? 1.0 : -1.0;
        walls_.push_back({1, s});
        if (lat.form == Form::Sheet) walls_.push_back({-1, -s});
        // The far face of a periodic axis copies the near one, so both get
        // the same patches and their vertices can be welded.
        for (size_t w = 0; w < walls_.size(); ++w) {
            phi_.emplace_back(size_t(sxy_ * (N_[2] + 1)));
#pragma omp parallel for schedule(static)
            for (int z = 0; z <= N_[2]; ++z) {
                int i[3] = {0, 0, z};
                for (i[1] = 0; i[1] <= N_[1]; ++i[1])
                    for (i[0] = 0; i[0] <= N_[0]; ++i[0]) {
                        int k[3] = {i[0], i[1], i[2]};
                        for (int a = 0; a < 3; ++a)
                            if ((periodic_ >> a & 1) && k[a] == N_[a]) k[a] = 0;
                        phi_[w][size_t(node(i))] = phi(int(w), nodePos(k));
                    }
            }
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
    // fit error in mm, or a negative value with `error` set.
    double run(int n, std::string& error) {
        if (offSurface_) {
            error = "a patch corner did not converge onto the surface";
            return -1;
        }
        spline_ = Spline(n);
        const Spline& sp = spline_;
        // Each edge and each patch writes only its own samples: in parallel.
        const int ne = int(edges_.size()), np = int(patches_.size());
#pragma omp parallel for schedule(dynamic, 64)
        for (int e = 0; e < ne; ++e)
            if (edges_[size_t(e)].wall >= 0) sampleEdge(edges_[size_t(e)], sp);
        std::vector<double> errs(patches_.size()), vols(patches_.size());
#pragma omp parallel for schedule(dynamic, 16)
        for (int k = 0; k < np; ++k) {
            errs[size_t(k)] = fitPatch(patches_[size_t(k)], sp);
            vols[size_t(k)] = patchVolume(patches_[size_t(k)]);
        }
        faces_.clear();
        double maxErr = 0;
        for (size_t k = 0; k < patches_.size(); ++k) {
            maxErr = std::max(maxErr, errs[k]);
            Face f;
            f.patch = int(k);
            Loop loop;
            for (int c = 0; c < 4; ++c)
                loop.push_back({patches_[k].e[size_t(c)], edges_[size_t(patches_[k].e[size_t(c)])].a == patches_[k].v[size_t(c)]});
            f.loops.push_back(loop);
            f.volume = vols[k];
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

    // Curved shapes: welds the two copies of each periodic face (the seam).
    // The far copy of every vertex and edge becomes the near one. A cap on a
    // curved surface keeps the seam as an edge it uses twice, as CAD seams
    // are; a flat cap that the seam cuts (the tube's ends) is joined across
    // it into faces with holes. Call after the last run().
    bool weld(std::string& error) {
        if (!periodic_) return true;
        // The copy of a vertex on the near face of the given axes.
        auto nearCopy = [&](int v, int axes) {
            const auto it = vertOf_.find(nearKey(verts_[size_t(v)], axes));
            return it == vertOf_.end() ? -1 : it->second;
        };
        std::vector<int> rep(verts_.size());
        for (size_t v = 0; v < verts_.size(); ++v)
            if ((rep[v] = nearCopy(int(v), farFaces(int(v)))) < 0) return error = "seam vertices do not match", false;
        // An edge on far faces is a copy of the edge between the near copies
        // of its ends. (Matching by welded ends would not do: on a torus cap
        // both seam circles run from the one welded corner to itself.)
        std::vector<std::pair<int, bool>> canon(edges_.size());  // (edge, same direction)
        for (size_t e = 0; e < edges_.size(); ++e) {
            const Edge& E = edges_[e];
            const int axes = farFaces(E.a) & farFaces(E.b), a = nearCopy(E.a, axes), b = nearCopy(E.b, axes);
            const auto it = a < 0 || b < 0 ? edgeOf_.end()
                                           : edgeOf_.find((uint64_t(std::min(a, b)) << 32) | uint64_t(std::max(a, b)));
            if (it == edgeOf_.end()) return error = "seam edges do not match", false;
            canon[e] = {it->second, edges_[size_t(it->second)].a == a};
        }
        for (Edge& E : edges_) E.a = rep[size_t(E.a)], E.b = rep[size_t(E.b)];
        for (Face& f : faces_)
            for (Loop& l : f.loops)
                for (auto& [e, same] : l) {
                    same = same == canon[size_t(e)].second;
                    e = canon[size_t(e)].first;
                }
        for (Patch& p : patches_)
            for (size_t k = 0; k < 4; ++k) p.v[k] = rep[size_t(p.v[k])], p.e[k] = canon[size_t(p.e[k])].first;

        for (int F = 0; F < 6; ++F) {
            if (!flatCap(F) || !(periodic_ & ~(1 << F / 2))) continue;
            std::vector<Directed> de;
            std::vector<Face> rest;
            for (const Face& f : faces_) {
                if (f.box != F) {
                    rest.push_back(f);
                    continue;
                }
                for (const Loop& l : f.loops)
                    for (const auto& [e, same] : l) {
                        const Edge& E = edges_[size_t(e)];
                        de.push_back(same ? Directed{E.a, E.b, e, true} : Directed{E.b, E.a, e, false});
                    }
            }
            std::unordered_map<int, int> uses;
            for (const Directed& d : de) ++uses[d.edge];
            de.erase(std::remove_if(de.begin(), de.end(), [&](const Directed& d) { return uses[d.edge] == 2; }),
                     de.end());  // seam edges run both ways
            faces_ = rest;
            if (!closeCaps(F, de, true, error)) return false;
        }
        return true;
    }

    // Merge patches into large trimmed faces. Each region is a height field
    // over a plane, so one B-spline surface holds it within tol. Chains of
    // patch edges between the same two faces become one curve within tol / 4.
    // Edges are then tolerant: they lie within tol of their faces, not on
    // them. On failure nothing changes.
    bool merge(double& maxErr, std::string& error) {
        const double tol = tol_;
        // Cone half-angle of a region's normals: 82 degrees gave the fewest
        // faces (wider cones fail the fit and get split again).
        const double kPi = 3.141592653589793, kStart = 82 * kPi / 180, kStep = 15 * kPi / 180;
        const size_t np = patches_.size();
        // Zero-length edges (two corners on one point, e.g. on a mirror plane)
        // cannot be projected onto a merged face: collapse them.
        rep_.resize(verts_.size());
        std::iota(rep_.begin(), rep_.end(), 0);
        auto root = [&](int v) {
            while (rep_[size_t(v)] != v) v = rep_[size_t(v)] = rep_[size_t(rep_[size_t(v)])];
            return v;
        };
        tiny_.assign(edges_.size(), 0);
        for (size_t e = 0; e < edges_.size(); ++e) {
            const auto& p = edges_[e].pts;  // closed edges have equal ends: measure the path
            double length = 0;
            for (size_t k = 0; k + 1 < p.size(); ++k) length += norm(p[k + 1] - p[k]);
            if (length < 1e-6 * hmin_) {
                tiny_[e] = 1;
                rep_[size_t(root(edges_[e].a))] = root(edges_[e].b);
            }
        }
        for (size_t v = 0; v < rep_.size(); ++v) rep_[v] = root(int(v));
        edgeFaces_.assign(edges_.size(), {});
        for (size_t f = 0; f < faces_.size(); ++f)
            for (const Loop& l : faces_[f].loops)
                for (const auto& [e, same] : l) edgeFaces_[size_t(e)].push_back({int(f), same});
        pn_.assign(np, {});
        pnbr_.assign(np, {});
#pragma omp parallel for schedule(dynamic, 64)
        for (int k = 0; k < int(np); ++k) {
            const Patch& p = patches_[size_t(k)];
            for (const auto& row : p.pts)
                for (const Vec3& q : row) pn_[size_t(k)].push_back(normal(p.wall, q));
            for (int e : p.e)
                for (const auto& [f, same] : edgeFaces_[size_t(e)])
                    if (f != k && f < int(np)) pnbr_[size_t(k)].push_back(f);
        }

        // Regions; a region that does not fit is clustered again with a
        // smaller cone, down to single patches.
        std::vector<std::pair<std::vector<int>, double>> work;
        for (int w = 0; w < int(walls_.size()); ++w) {
            std::vector<int> members;
            for (size_t k = 0; k < np; ++k)
                if (patches_[k].wall == w) members.push_back(int(k));
            for (auto& c : cluster(members, std::cos(kStart))) work.push_back({c, kStart});
        }
        // In rounds: the regions of a round are fitted in parallel, then the
        // ones that did not fit are split for the next round.
        std::vector<Face> regions;
        std::vector<int> regionOf(faces_.size(), -1);
        maxErr = 0;
        while (!work.empty()) {
            const int nw = int(work.size());
            std::vector<Face> fit(work.size());
            std::vector<double> err(work.size(), 0.0);
            std::vector<char> fitted(work.size());
#pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < nw; ++i)
                fitted[size_t(i)] = fitRegion(work[size_t(i)].first, tol, fit[size_t(i)], err[size_t(i)]);
            std::vector<std::pair<std::vector<int>, double>> next;
            for (size_t i = 0; i < work.size(); ++i) {
                const auto& [members, angle] = work[i];
                maxErr = std::max(maxErr, err[i]);
                if (!fitted[i] && members.size() > 1) {
                    const double smaller = angle - kStep;
                    if (smaller > 0)
                        for (auto& c : cluster(members, std::cos(smaller))) next.push_back({c, smaller});
                    else
                        for (int p : members) next.push_back({{p}, smaller});
                    continue;
                }
                // A single patch too curved for a height field keeps its own surface.
                for (int p : members) regionOf[size_t(p)] = int(regions.size());
                regions.push_back(fitted[i] ? fit[i] : faces_[size_t(members[0])]);
            }
            work.swap(next);
        }
        for (size_t f = np; f < faces_.size(); ++f) regionOf[f] = int(regions.size() + f);

        // Chains: maximal runs of boundary edges with the same face pair.
        auto pairOf = [&](int e) {
            const auto& u = edgeFaces_[size_t(e)];
            const int r0 = regionOf[size_t(u[0].first)], r1 = regionOf[size_t(u[1].first)];
            return std::make_pair(std::min(r0, r1), std::max(r0, r1));
        };
        std::unordered_map<int, std::vector<int>> atVertex;
        for (size_t e = 0; e < edges_.size(); ++e) {
            if (edgeFaces_[e].size() != 2 || tiny_[e]) continue;
            // An edge inside a region goes; a seam edge used twice by a cap stays.
            const auto pr = pairOf(int(e));
            if (pr.first == pr.second && pr.first < int(regions.size())) continue;
            atVertex[rep_[size_t(edges_[e].a)]].push_back(int(e));
            atVertex[rep_[size_t(edges_[e].b)]].push_back(int(e));
        }
        auto away = [&](int e, int v) {  // unit tangent of edge e leaving vertex v
            const auto& p = edges_[size_t(e)].mm;
            const Vec3 d = rep_[size_t(edges_[size_t(e)].a)] == v ? p[1] - p[0] : p[p.size() - 2] - p.back();
            return d * (1.0 / norm(d));
        };
        // v is inside a chain of pair pr: two edges of that pair, and no kink
        // sharper than 30 degrees (one curve cannot follow a corner).
        // A patch that kept its own surface keeps its own edges too: they lie
        // exactly on it, and a merged curve would need projecting onto it.
        auto kept = [&](std::pair<int, int> pr) {
            for (int r : {pr.first, pr.second})
                if (r < int(regions.size()) && regions[size_t(r)].patch >= 0) return true;
            return false;
        };
        auto inner = [&](int v, std::pair<int, int> pr) {
            if (kept(pr)) return false;
            const auto& es = atVertex.at(v);
            return es.size() == 2 && pairOf(es[0]) == pr && pairOf(es[1]) == pr &&
                   dot(away(es[0], v), away(es[1], v)) < -0.866;
        };
        auto end = [&](int e, bool b) { return rep_[size_t(b ? edges_[size_t(e)].b : edges_[size_t(e)].a)]; };
        auto across = [&](int e, int v) { return end(e, false) == v ? end(e, true) : end(e, false); };
        auto otherEdge = [&](int v, int e) {
            const auto& es = atVertex.at(v);
            return es[0] == e ? es[1] : es[0];
        };
        std::vector<Edge> chains;
        std::vector<int> chainOf(edges_.size(), -1);
        std::vector<bool> chainSame(edges_.size(), false);
        for (const auto& [v0, es0] : atVertex)
            for (int e0 : es0) {
                if (chainOf[size_t(e0)] >= 0) continue;
                const auto pr = pairOf(e0);
                int v = end(e0, false), e = e0;  // rewind to the chain start
                bool closed = false;
                while (inner(v, pr)) {
                    e = otherEdge(v, e);
                    if (e == e0) {
                        closed = true;
                        break;
                    }
                    v = across(e, v);
                }
                if (closed) v = end(e0, false), e = e0;
                Edge c{v, v, edges_[size_t(e)].wall, {}, {}, {}, {}, {}};
                const int id = int(chains.size());
                for (;;) {
                    const bool fwd = end(e, false) == v;
                    chainOf[size_t(e)] = id;
                    chainSame[size_t(e)] = fwd;
                    std::vector<Vec3> s = edges_[size_t(e)].mm;
                    if (!fwd) std::reverse(s.begin(), s.end());
                    c.mm.insert(c.mm.end(), s.begin() + (c.mm.empty() ? 0 : 1), s.end());
                    v = across(e, v);
                    if (closed ? v == c.a : !inner(v, pr)) break;
                    e = otherEdge(v, e);
                }
                c.b = v;
                if (c.wall < 0) fitLine(c);
                else if (kept(pr)) {
                    c.ctrl = edges_[size_t(e)].ctrl;
                    if (!chainSame[size_t(e)]) std::reverse(c.ctrl.begin(), c.ctrl.end());
                    c.knots = spline_.knots;
                }
                else
                    fitCurve(c.mm, 0.25 * tol, c.ctrl, c.knots);
                chains.push_back(c);
            }

        auto convert = [&](const Loop& l) {
            Loop out;
            for (const auto& [e, same] : l) {
                const std::pair<int, bool> step{chainOf[size_t(e)], same == chainSame[size_t(e)]};
                if (step.first < 0 || (!out.empty() && out.back() == step)) continue;  // collapsed edge
                out.push_back(step);
            }
            if (out.size() > 1 && out.front() == out.back()) out.pop_back();
            return out;
        };
        std::vector<Face> faces;
        for (Face& r : regions) {
            for (Loop& l : r.loops) l = convert(l);
            faces.push_back(r);
        }
        for (size_t f = np; f < faces_.size(); ++f) {
            Face cap = faces_[f];
            for (Loop& l : cap.loops) l = convert(l);
            faces.push_back(cap);
        }
        std::swap(faces, faces_);
        std::swap(chains, edges_);
        std::vector<std::pair<std::vector<int>, double>> check;
        if (!shells(check, error)) {
            std::swap(faces, faces_);
            std::swap(chains, edges_);
            return false;
        }
        return true;
    }

    std::vector<long> write(StepFile& st, const std::vector<std::vector<int>>& shells) const {
        std::unordered_map<int, long> vertexId, edgeId, surfaceId;

        auto point = [&](Vec3 p) { return st.add("CARTESIAN_POINT(''," + triple(p) + ")"); };
        auto bspline = [&](const std::vector<Vec3>& ctrl, const std::vector<double>& U) {
            std::vector<long> cps;
            for (const Vec3& c : ctrl) cps.push_back(point(c));
            const auto [m, k] = knotLists(U);
            return st.add("B_SPLINE_CURVE_WITH_KNOTS('',3," + refList(cps) + ",.UNSPECIFIED.,.F.,.F.," + m + "," + k +
                          ",.UNSPECIFIED.)");
        };
        auto vertex = [&](int v) {
            auto it = vertexId.find(v);
            if (it != vertexId.end()) return it->second;
            const long id = st.add("VERTEX_POINT(''," + ref(point(mm(verts_[size_t(v)].p))) + ")");
            vertexId.emplace(v, id);
            return id;
        };
        auto surfaceOf = [&](int fi) {
            auto it = surfaceId.find(fi);
            if (it != surfaceId.end()) return it->second;
            const Face& f = faces_[size_t(fi)];
            long surface;
            if (f.patch >= 0 || !f.ctrl.empty()) {
                const auto& net = f.patch >= 0 ? patches_[size_t(f.patch)].ctrl : f.ctrl;
                const auto [mu, ku] = knotLists(f.patch >= 0 ? spline_.knots : f.knots[0]);
                const auto [mv, kv] = knotLists(f.patch >= 0 ? spline_.knots : f.knots[1]);
                std::string rows = "(";
                for (const auto& row : net) {
                    std::vector<long> ids;
                    for (const Vec3& c : row) ids.push_back(point(c));
                    rows += (rows.size() > 1 ? ",\n" : "") + refList(ids);
                }
                surface = st.add("B_SPLINE_SURFACE_WITH_KNOTS('',3,3," + rows + "),.UNSPECIFIED.,.F.,.F.,.F.," + mu +
                                 "," + mv + "," + ku + "," + kv + ",.UNSPECIFIED.)");
            } else {
                const int a = f.box / 2, side = f.box % 2;
                auto placement = [&](Vec3 origin, Vec3 z, Vec3 x) {
                    return ref(st.add("AXIS2_PLACEMENT_3D(''," + ref(point(origin)) + "," +
                                      ref(st.add("DIRECTION(''," + triple(z) + ")")) + "," +
                                      ref(st.add("DIRECTION(''," + triple(x) + ")")) + ")"));
                };
                if (!flatCap(f.box)) {  // radial cap of a curved shape, around z
                    const std::string ax = placement({}, {0, 0, 1}, {1, 0, 0}), r = real(side ? shape_.r1 : shape_.r0);
                    surface = shape_.kind == ShapeKind::Tube
                                  ? st.add("CYLINDRICAL_SURFACE(''," + ax + "," + r + ")")
                                  : st.add("TOROIDAL_SURFACE(''," + ax + "," + real(shape_.major) + "," + r + ")");
                } else {  // plane through a corner of the face, normal out of the solid
                    Vec3 corner, nrm, x;  // flat caps lie across a parameter axis that maps to the same mm axis
                    corner[a] = side ? L_[a] : 0.0;
                    nrm[a] = side ? 1.0 : -1.0;
                    x[side ? (a + 1) % 3 : (a + 2) % 3] = 1;
                    surface = st.add("PLANE(''," + placement(mm(corner), nrm, x) + ")");
                }
            }
            surfaceId.emplace(fi, surface);
            return surface;
        };
        auto edgeCurve = [&](int e) {
            auto it = edgeId.find(e);
            if (it != edgeId.end()) return it->second;
            const Edge& E = edges_[size_t(e)];
            long curve;
            if (E.ctrl.empty()) {
                const Vec3 d = E.mm.back() - E.mm.front();
                const long dir = st.add("DIRECTION(''," + triple(d * (1.0 / norm(d))) + ")");
                const long vec = st.add("VECTOR(''," + ref(dir) + "," + real(norm(d)) + ")");
                curve = st.add("LINE(''," + ref(point(E.mm.front())) + "," + ref(vec) + ")");
            } else {
                // No pcurves: written ones made OpenCASCADE's SameParameter
                // repair run away (50 min, 22 GB) on Fischer-Koch S.
                curve = bspline(E.ctrl, E.knots.empty() ? spline_.knots : E.knots);
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
                const long surface = surfaceOf(fi);
                std::vector<long> bounds;
                for (size_t k = 0; k < f.loops.size(); ++k)
                    bounds.push_back(st.add(std::string(k ? "FACE_BOUND" : "FACE_OUTER_BOUND") + "(''," +
                                            ref(loop(f.loops[k])) + ",.T.)"));
                // The inner radial cap faces against its surface's normal.
                const bool against = f.box == 0 && !flatCap(0);
                faces.push_back(st.add("ADVANCED_FACE(''," + refList(bounds) + "," + ref(surface) +
                                       (against ? ",.F.)" : ",.T.)")));
            }
            solids.push_back(st.add("MANIFOLD_SOLID_BREP(''," + ref(st.add("CLOSED_SHELL(''," + refList(faces) + ")")) + ")"));
        }
        return solids;
    }

private:
    struct Wall {
        double levelSign, sign;  // phi = sign * (f - levelSign c), inside where phi < 0
    };

    const Lattice& lat_;
    const Shape& shape_;
    const bool mapped_;  // built in a parameter box (see Shape)
    const int periodic_;  // bit per parameter axis that wraps around
    const double tol_;
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
    std::vector<std::vector<std::pair<int, bool>>> edgeFaces_;  // (face, same) using each edge
    std::vector<std::vector<Vec3>> pn_;                          // unit normals at patch samples
    std::vector<std::vector<int>> pnbr_;                         // patches sharing an edge
    std::vector<int> rep_;                                       // vertex after collapsing zero-length edges
    std::vector<char> tiny_;                                     // zero-length edges

    // Periodic axes on whose far face vertex v lies.
    int farFaces(int v) const {
        int axes = 0;
        for (int b = 0; b < 3; ++b)
            if ((periodic_ >> b & 1) && (verts_[size_t(v)].faces >> (2 * b + 1) & 1)) axes |= 1 << b;
        return axes;
    }

    // Key of the copy of a vertex on the near face of the given axes (the
    // vertex lies on their far faces). The grid repeats there, so the copy
    // has the same crossings and component.
    uint64_t nearKey(const Vert& v, int axes) const {
        const bool corner = v.key >> 56 == 15;
        const int64_t stride[3] = {1, sx_, sxy_};
        uint64_t key = v.key;
        for (int b = 0; b < 3; ++b)
            if (axes >> b & 1) key -= corner ? uint64_t(1) << b : uint64_t(N_[b] * stride[b]) << 5;
        return key;
    }

    Vec3 mm(Vec3 p) const { return mapped_ ? shape_.map(p) : p; }
    std::vector<Vec3> mapped(std::vector<Vec3> ps) const {
        for (Vec3& p : ps) p = mm(p);
        return ps;
    }

    // Distance estimate of a point in mm from wall w: |phi| / |grad phi|.
    double distance(int w, Vec3 q) const {
        const Vec3 p = mapped_ ? shape_.unmap(q) : q, g = grad(w, p);
        return std::fabs(phi(w, p)) / norm(mapped_ ? shape_.gradient(p, g) : g);
    }

    // Unit normal in mm, out of the solid, at a point of wall w.
    Vec3 normal(int w, Vec3 p) const {
        Vec3 g = grad(w, p);
        if (mapped_) g = shape_.gradient(p, g);
        return g * (1.0 / norm(g));
    }

    // Box face F is capped by a plane (else a cylinder or torus around z).
    bool flatCap(int F) const { return shape_.kind == ShapeKind::Box || (shape_.kind == ShapeKind::Tube && F / 2 == 2); }

    // Agglomerative clustering of adjacent patches: merge the neighbor pair
    // whose union has the narrowest normal cone, while every normal stays
    // within acos(cosMax) of the union's mean normal.
    std::vector<std::vector<int>> cluster(const std::vector<int>& members, double cosMax) const {
        std::unordered_map<int, int> of;  // patch -> cluster
        std::vector<std::vector<int>> cl;
        std::vector<Vec3> sum;
        std::vector<int> version;
        for (int p : members) {
            of[p] = int(cl.size());
            cl.push_back({p});
            Vec3 s;
            for (const Vec3& n : pn_[size_t(p)]) s = s + n;
            sum.push_back(s);
            version.push_back(0);
        }
        struct Cand {
            double cosMin;
            int a, b, va, vb;
            bool operator<(const Cand& o) const { return cosMin < o.cosMin; }
        };
        std::priority_queue<Cand> pq;
        auto push = [&](int a, int b) {
            Vec3 ax = sum[size_t(a)] + sum[size_t(b)];
            ax = ax * (1.0 / norm(ax));
            double c = 1;
            for (int k : {a, b})
                for (int p : cl[size_t(k)])
                    for (const Vec3& n : pn_[size_t(p)]) c = std::min(c, dot(n, ax));
            if (c >= cosMax) pq.push({c, a, b, version[size_t(a)], version[size_t(b)]});
        };
        for (int p : members)
            for (int q : pnbr_[size_t(p)]) {
                const auto it = of.find(q);
                if (it != of.end() && p < q) push(of[p], it->second);
            }
        while (!pq.empty()) {
            const Cand c = pq.top();
            pq.pop();
            if (version[size_t(c.a)] != c.va || version[size_t(c.b)] != c.vb) continue;
            int a = c.a, b = c.b;
            if (cl[size_t(a)].size() < cl[size_t(b)].size()) std::swap(a, b);
            for (int p : cl[size_t(b)]) {
                of[p] = a;
                cl[size_t(a)].push_back(p);
            }
            cl[size_t(b)].clear();
            sum[size_t(a)] = sum[size_t(a)] + sum[size_t(b)];
            ++version[size_t(a)], ++version[size_t(b)];
            std::vector<int> nb;
            for (int p : cl[size_t(a)])
                for (int q : pnbr_[size_t(p)]) {
                    const auto it = of.find(q);
                    if (it != of.end() && it->second != a) nb.push_back(it->second);
                }
            std::sort(nb.begin(), nb.end());
            nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
            for (int o : nb) push(a, o);
        }
        std::vector<std::vector<int>> out;
        for (auto& c : cl)
            if (!c.empty()) out.push_back(c);
        return out;
    }

    // One trimmed height-field B-spline for a set of patches. Fails if the
    // region touches itself at a vertex, is not a height field over its mean
    // plane, or does not fit within tol at the data and 2 tol at check points
    // between them, with up to about 60 x 60 control points.
    bool fitRegion(const std::vector<int>& members, double tol, Face& out, double& maxErr) const {
        const int np = int(patches_.size());
        std::vector<char> in(size_t(np), 0);
        for (int p : members) in[size_t(p)] = 1;

        struct Step {
            int to, edge;
            bool same;
        };
        std::unordered_map<int, Step> next;  // boundary edges, by start vertex
        for (int p : members)
            for (size_t k = 0; k < 4; ++k) {
                const int e = patches_[size_t(p)].e[k];
                if (tiny_[size_t(e)]) continue;
                bool internal = false;
                for (const auto& [f, s] : edgeFaces_[size_t(e)]) internal |= f != p && f < np && in[size_t(f)];
                if (internal) continue;
                const int from = rep_[size_t(patches_[size_t(p)].v[k])], to = rep_[size_t(patches_[size_t(p)].v[(k + 1) % 4])];
                if (!next.emplace(from, Step{to, e, rep_[size_t(edges_[size_t(e)].a)] == from}).second) return false;
            }
        std::vector<Loop> loops;
        std::unordered_map<int, bool> used;
        for (const auto& [start, first] : next) {
            if (used.count(start)) continue;
            Loop l;
            int v = start;
            do {
                used[v] = true;
                const auto it = next.find(v);
                if (it == next.end()) return false;
                l.push_back({it->second.edge, it->second.same});
                v = it->second.to;
            } while (!used.count(v));
            if (v != start) return false;
            loops.push_back(l);
        }

        Vec3 axis;
        for (int p : members)
            for (const Vec3& n : pn_[size_t(p)]) axis = axis + n;
        axis = axis * (1.0 / norm(axis));
        for (int p : members)
            for (const Vec3& n : pn_[size_t(p)])
                if (dot(n, axis) < 0.1) return false;
        Vec3 e1 = cross(axis, std::fabs(axis.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0});
        e1 = e1 * (1.0 / norm(e1));
        Vec3 e2 = cross(axis, e1);

        // Data: points on the exact surface, at the patch samples and next to
        // the patch cell mid points. So the fit is within tol of the surface
        // there, the same check as for a patch.
        std::vector<Vec3> data;
        for (int p : members) {
            const Patch& pa = patches_[size_t(p)];
            for (const auto& row : pa.mm) data.insert(data.end(), row.begin(), row.end());
            data.insert(data.end(), pa.mids.begin(), pa.mids.end());
        }
        Vec3 o;
        for (const Vec3& d : data) o = o + d;
        o = o * (1.0 / double(data.size()));
        // The rectangle follows the principal directions of the data, so it
        // stays tight around long, thin regions; e1 the longer one.
        double cuu = 0, cvv = 0, cuv = 0;
        for (const Vec3& d : data) {
            const double u = dot(d - o, e1), v = dot(d - o, e2);
            cuu += u * u, cvv += v * v, cuv += u * v;
        }
        const double turn = 0.5 * std::atan2(2 * cuv, cuu - cvv);
        e1 = e1 * std::cos(turn) + e2 * std::sin(turn);
        e2 = cross(axis, e1);
        double u0 = 1e300, u1 = -1e300, v0 = 1e300, v1 = -1e300;
        std::vector<std::array<double, 3>> uvh;
        for (const Vec3& d : data) {
            const Vec3 q = d - o;
            uvh.push_back({dot(q, e1), dot(q, e2), dot(q, axis)});
            u0 = std::min(u0, uvh.back()[0]), u1 = std::max(u1, uvh.back()[0]);
            v0 = std::min(v0, uvh.back()[1]), v1 = std::max(v1, uvh.back()[1]);
        }
        const double margin = 0.02 * std::max(u1 - u0, v1 - v0) + 1e-9;
        u0 -= margin, u1 += margin, v0 -= margin, v1 += margin;
        for (auto& x : uvh) x[0] = (x[0] - u0) / (u1 - u0), x[1] = (x[1] - v0) / (v1 - v0);
        // Check points between the data: the patch splines at half steps,
        // (s, t, wall). The fit must stay near the exact surface there too.
        std::vector<std::array<double, 3>> checks;
        for (int pm : members) {
            const Patch& pa = patches_[size_t(pm)];
            for (const auto& row : halfSteps(pa.ctrl))
                for (const Vec3& f : row)
                    checks.push_back({(dot(f - o, e1) - u0) / (u1 - u0), (dot(f - o, e2) - v0) / (v1 - v0), double(pa.wall)});
        }

        // A Ku x Kv net with about K^2 points, spaced alike along u and v.
        constexpr int kMax = 60;
        constexpr double kDrift = 2;  // between the data, the fit may be up to 2 tol off
        const double aspect = std::sqrt(std::min(16.0, (u1 - u0) / (v1 - v0)));
        auto uniform = [](int k) {
            std::vector<double> inner;
            for (int j = 1; j <= k - 4; ++j) inner.push_back(double(j) / (k - 3));
            return clamped(inner);
        };
        for (int K = std::min(kMax, int(3 * std::sqrt(double(members.size()))) + 4);; K = std::min(kMax, K + K / 5 + 1)) {
            const int Ku = std::max(4, int(K * aspect + 0.5)), Kv = std::max(4, int(K / aspect + 0.5));
            const std::vector<double> U = uniform(Ku), V = uniform(Kv);
            const size_t n = size_t(Ku) * size_t(Kv), w = size_t(3 * Kv + 3);
            std::vector<double> A(n * (w + 1), 0.0), h(n, 0.0);
            auto add = [&](const size_t* idx, const double* c, int cnt, double weight, double rhs) {
                for (int i = 0; i < cnt; ++i) {
                    h[idx[i]] += weight * c[i] * rhs;
                    for (int j = 0; j < cnt; ++j)
                        if (idx[j] >= idx[i]) A[idx[i] * (w + 1) + (idx[j] - idx[i])] += weight * c[i] * c[j];
                }
            };
            for (const auto& x : uvh) {
                double Ns[4], Nt[4], c[16];
                size_t idx[16];
                const int is = cubicBasis(U, x[0], Ns), it = cubicBasis(V, x[1], Nt);
                for (int p = 0; p < 4; ++p)
                    for (int q = 0; q < 4; ++q) {
                        idx[p * 4 + q] = size_t(is + p) * size_t(Kv) + size_t(it + q);
                        c[p * 4 + q] = Ns[p] * Nt[q];
                    }
                add(idx, c, 16, 1.0, x[2]);
            }
            // Small second-difference penalty: fills control points that no
            // data point reaches, barely touches the rest.
            const double second[3] = {1, -2, 1};
            for (int a = 0; a < Ku; ++a)
                for (int b = 0; b < Kv; ++b) {
                    const size_t i = size_t(a * Kv + b);
                    const size_t alongV[3] = {i - 1, i, i + 1}, alongU[3] = {i - size_t(Kv), i, i + size_t(Kv)};
                    if (b > 0 && b + 1 < Kv) add(alongV, second, 3, 1e-5, 0.0);
                    if (a > 0 && a + 1 < Ku) add(alongU, second, 3, 1e-5, 0.0);
                }
            if (!bandSolve(std::move(A), n, w, h)) return false;

            // Distance of each data point from the fitted surface: the height
            // difference over sqrt(1 + |grad h|^2), the distance to the tangent
            // plane. It is measured to the data, so a fit cannot pass by
            // touching another sheet. (The height difference alone overstates
            // it up to 7 times where the face is steep, which made nets dense.)
            auto height = [&](double s, double t) {
                double Ns[4], Nt[4], z = 0;
                const int is = cubicBasis(U, s, Ns), it = cubicBasis(V, t, Nt);
                for (int p = 0; p < 4; ++p)
                    for (int q = 0; q < 4; ++q) z += Ns[p] * Nt[q] * h[size_t(is + p) * size_t(Kv) + size_t(it + q)];
                return z;
            };
            const double d = 1e-5, du = 2 * d * (u1 - u0), dv = 2 * d * (v1 - v0);
            double err = 0;
            for (size_t k = 0; k < uvh.size() && err <= tol; ++k) {
                const double s = uvh[k][0], t = uvh[k][1];
                const double hu = (height(s + d, t) - height(s - d, t)) / du, hv = (height(s, t + d) - height(s, t - d)) / dv;
                err = std::max(err, std::fabs(height(s, t) - uvh[k][2]) / std::sqrt(1 + hu * hu + hv * hv));
            }
            double drift = 0;
            for (size_t k = 0; k < checks.size() && err <= tol && drift <= kDrift * tol; ++k) {
                const double s = checks[k][0], t = checks[k][1];
                const Vec3 r = o + e1 * (u0 + s * (u1 - u0)) + e2 * (v0 + t * (v1 - v0)) + axis * height(s, t);
                drift = std::max(drift, distance(int(checks[k][2]), r));
            }
            if (err > tol || drift > kDrift * tol) {
                if (K == kMax) return false;
                continue;
            }
            maxErr = std::max(maxErr, err);

            // Control points: x and y at the Greville abscissae reproduce the
            // linear map exactly, so the surface is o + u e1 + v e2 + h axis.
            out.ctrl.assign(size_t(Ku), std::vector<Vec3>(size_t(Kv)));
            for (int a = 0; a < Ku; ++a)
                for (int b = 0; b < Kv; ++b) {
                    const double ga = (U[size_t(a + 1)] + U[size_t(a + 2)] + U[size_t(a + 3)]) / 3;
                    const double gb = (V[size_t(b + 1)] + V[size_t(b + 2)] + V[size_t(b + 3)]) / 3;
                    out.ctrl[size_t(a)][size_t(b)] = o + e1 * (u0 + ga * (u1 - u0)) + e2 * (v0 + gb * (v1 - v0)) +
                                                     axis * h[size_t(a * Kv + b)];
                }
            out.knots[0] = U, out.knots[1] = V;
            break;
        }

        // The outer loop runs counter-clockwise in (u, v), holes clockwise.
        int outer = -1;
        for (size_t k = 0; k < loops.size(); ++k) {
            double area = 0;
            std::vector<std::array<double, 2>> poly;
            for (const auto& [e, same] : loops[k]) {
                const auto& pts = edges_[size_t(e)].mm;
                for (size_t q = 0; q + 1 < pts.size(); ++q) {
                    const Vec3 d = pts[same ? q : pts.size() - 1 - q] - o;
                    poly.push_back({dot(d, e1), dot(d, e2)});
                }
            }
            for (size_t q = 0; q < poly.size(); ++q) {
                const auto& A = poly[q];
                const auto& B = poly[(q + 1) % poly.size()];
                area += A[0] * B[1] - B[0] * A[1];
            }
            if (area > 0) {
                if (outer >= 0) return false;
                outer = int(k);
            }
        }
        if (outer < 0) return false;
        std::swap(loops[0], loops[size_t(outer)]);
        out.loops = loops;
        out.volume = 0;
        for (int p : members) out.volume += faces_[size_t(p)].volume;
        return true;
    }

    int64_t node(const int i[3]) const { return i[0] + sx_ * i[1] + sxy_ * i[2]; }
    double nodePhi(int w, const int i[3]) const { return phi_[size_t(w)][size_t(node(i))]; }
    Vec3 nodePos(const int i[3]) const { return {i[0] * H_[0], i[1] * H_[1], i[2] * H_[2]}; }

    double phi(int w, Vec3 p) const {
        const Wall& W = walls_[size_t(w)];
        return W.sign * (lat_.value(p.x, p.y, p.z) - W.levelSign * lat_.levelAt(p.x, p.y, p.z));
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
        verts_.push_back({p, faces, key});
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
        verts_.push_back({p, faces, key});
        vertOf_.emplace(key, int(verts_.size()) - 1);
        return int(verts_.size()) - 1;
    }

    int edge(int a, int b, int wall) {
        const uint64_t key = (uint64_t(std::min(a, b)) << 32) | uint64_t(std::max(a, b));
        if (auto it = edgeOf_.find(key); it != edgeOf_.end()) return it->second;
        edges_.push_back({std::min(a, b), std::max(a, b), wall, {}, {}, {}, {}, {}});
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
        e.mm = mapped(e.pts);
        e.ctrl = sp.fit(e.mm);
        e.pctrl = sp.fit(e.pts);
    }

    // A segment of a box edge: straight in the parameter box. In mm a line,
    // or on a curved shape a curve through 32 samples per cell.
    void sampleLine(Edge& e) const {
        const Vec3 A = verts_[size_t(e.a)].p, B = verts_[size_t(e.b)].p;
        const int k = shape_.bent() ? std::max(8, int(32 * norm(B - A))) : 1;
        e.pts.clear();
        for (int i = 0; i <= k; ++i) e.pts.push_back(A + (B - A) * (double(i) / k));
        e.mm = mapped(e.pts);
        fitLine(e);
    }

    // A line if the mm samples are straight, else a curve within tol / 4.
    void fitLine(Edge& e) const {
        const Vec3 a = e.mm.front(), d = e.mm.back() - a;
        bool straight = norm(d) > 0;  // closed: a curve
        for (const Vec3& q : e.mm) straight = straight && norm(cross(q - a, d)) <= 1e-9 * dot(d, d);
        e.knots.clear();
        if (!straight) return fitCurve(e.mm, 0.25 * tol_, e.ctrl, e.knots);
        e.mm = {a, e.mm.back()};
        e.ctrl.clear();
    }

    std::vector<Vec3> side(const Patch& p, int k, bool inMm) const {  // samples from v[k] to v[k + 1]
        const Edge& e = edges_[size_t(p.e[size_t(k)])];
        std::vector<Vec3> s = inMm ? e.mm : e.pts;
        if (e.a != p.v[size_t(k)]) std::reverse(s.begin(), s.end());
        return s;
    }

    // Control points of the bicubic spline through a grid of samples.
    std::vector<std::vector<Vec3>> interpolate(const std::vector<std::vector<Vec3>>& pts, const Spline& sp) const {
        const size_t m = pts.size();
        std::vector<std::vector<Vec3>> tmp(m, std::vector<Vec3>(m)), ctrl(m, std::vector<Vec3>(m));
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < m; ++j)
                for (size_t k = 0; k < m; ++k) tmp[i][j] = tmp[i][j] + pts[k][j] * sp.inv[i][k];
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < m; ++j)
                for (size_t k = 0; k < m; ++k) ctrl[i][j] = ctrl[i][j] + tmp[i][k] * sp.inv[j][k];
        return ctrl;
    }

    // Coons grid from the four boundary curves, projected onto the surface,
    // then interpolated. Returns the max distance at the cell mid points.
    double fitPatch(Patch& p, const Spline& sp) {
        const int n = sp.n;
        const size_t m = size_t(n) + 1;
        // Boundary samples from the edges; the mm boundary comes from the
        // edges' own mm samples, so the patch boundary equals the edge curve.
        for (bool inMm : {false, true}) {
            const auto B = side(p, 0, inMm), R = side(p, 1, inMm);
            auto T = side(p, 2, inMm), Lf = side(p, 3, inMm);
            std::reverse(T.begin(), T.end());    // v3 -> v2
            std::reverse(Lf.begin(), Lf.end());  // v0 -> v3
            auto& grid = inMm ? p.mm : p.pts;
            grid.assign(m, std::vector<Vec3>(m));
            for (size_t i = 0; i < m; ++i)
                for (size_t j = 0; j < m; ++j) {
                    Vec3& P = grid[i][j];
                    if (j == 0) P = B[i];
                    else if (j == m - 1) P = T[i];
                    else if (i == 0) P = Lf[j];
                    else if (i == m - 1) P = R[j];
                    else if (inMm) P = mm(p.pts[i][j]);
                    else {
                        const double u = double(i) / n, v = double(j) / n;
                        P = B[i] * (1 - v) + T[i] * v + Lf[j] * (1 - u) + R[j] * u -
                            (B[0] * ((1 - u) * (1 - v)) + B[m - 1] * (u * (1 - v)) + T[m - 1] * (u * v) +
                             T[0] * ((1 - u) * v));
                        P = project(p.wall, P, 0);
                    }
                }
        }
        p.ctrl = interpolate(p.mm, sp);
        p.pctrl = interpolate(p.pts, sp);

        // Distance at the cell mid points. Their projections onto the exact
        // surface are kept as data for merging.
        double err = 0;
        p.mids.clear();
        for (size_t i = 0; i < m - 1; ++i)
            for (size_t j = 0; j < m - 1; ++j) {
                Vec3 q;
                for (size_t k = 0; k < m; ++k)
                    for (size_t l = 0; l < m; ++l) q = q + p.ctrl[k][l] * (sp.mid[i][k] * sp.mid[j][l]);
                err = std::max(err, distance(p.wall, q));
                p.mids.push_back(mm(project(p.wall, mapped_ ? shape_.unmap(q) : q, 0)));
            }
        return err;
    }

    // Points along an edge: in the parameter box at the fine parameters (a
    // box edge segment: its samples), or in mm.
    std::vector<Vec3> finePoints(const Edge& e, bool inMm) const {
        if (e.wall < 0) return inMm ? e.mm : e.pts;
        const auto& ctrl = inMm ? e.ctrl : e.pctrl;
        std::vector<Vec3> out(spline_.fine.size());
        for (size_t i = 0; i < out.size(); ++i)
            for (size_t k = 0; k < ctrl.size(); ++k) out[i] = out[i] + ctrl[k] * spline_.fine[i][k];
        return out;
    }

    // A patch spline at half steps: the samples and the points between them.
    std::vector<std::vector<Vec3>> halfSteps(const std::vector<std::vector<Vec3>>& ctrl) const {
        const size_t m = ctrl.size(), f = 2 * m - 1;
        std::vector<std::vector<double>> B(f);
        for (size_t i = 0; i < f; ++i) B[i] = spline_.basis(double(i) / double(f - 1));
        std::vector<std::vector<Vec3>> tmp(f, std::vector<Vec3>(m)), s(f, std::vector<Vec3>(f));
        for (size_t i = 0; i < f; ++i)
            for (size_t l = 0; l < m; ++l)
                for (size_t k = 0; k < m; ++k) tmp[i][l] = tmp[i][l] + ctrl[k][l] * B[i][k];
        for (size_t i = 0; i < f; ++i)
            for (size_t j = 0; j < f; ++j)
                for (size_t l = 0; l < m; ++l) s[i][j] = s[i][j] + tmp[i][l] * B[j][l];
        return s;
    }

    // A patch spline at the fine parameters.
    std::vector<std::vector<Vec3>> finePatch(const std::vector<std::vector<Vec3>>& ctrl) const {
        const auto& F = spline_.fine;
        const size_t f = F.size(), m = ctrl.size();
        std::vector<std::vector<Vec3>> tmp(f, std::vector<Vec3>(m)), s(f, std::vector<Vec3>(f));
        for (size_t i = 0; i < f; ++i)
            for (size_t l = 0; l < m; ++l)
                for (size_t k = 0; k < m; ++k) tmp[i][l] = tmp[i][l] + ctrl[k][l] * F[i][k];
        for (size_t i = 0; i < f; ++i)
            for (size_t j = 0; j < f; ++j)
                for (size_t l = 0; l < m; ++l) s[i][j] = s[i][j] + tmp[i][l] * F[j][l];
        return s;
    }

    // Volume term of a patch: the flux (Phi, 0, 0) through its spline in the
    // parameter box, on triangles at the fine parameters. Its boundary points
    // equal finePoints() of its edges, so with the caps the triangles close a
    // polyhedron.
    double patchVolume(const Patch& p) const {
        const auto s = finePatch(p.pctrl);
        double v = 0;
        auto triangle = [&](Vec3 a, Vec3 b, Vec3 c) { v += shape_.flux((a + b + c) * (1.0 / 3)) * cross(b - a, c - a).x; };
        for (size_t i = 0; i + 1 < s.size(); ++i)
            for (size_t j = 0; j + 1 < s.size(); ++j) {
                triangle(s[i][j], s[i + 1][j], s[i + 1][j + 1]);
                triangle(s[i][j], s[i + 1][j + 1], s[i][j + 1]);
            }
        return v / 2;
    }

    // Volume term of a cap polygon (s, t) on box face F: the flux of
    // (Phi, 0, 0) through it. Only the face p.x = L_x has any: Phi = 0 at
    // p.x = 0, and the other faces are parallel to the flux. By Green, the
    // integral of Phi over the polygon is the loop integral of Psi dt
    // (Shape::fluxArea), with (s, t) = (p.y, p.z).
    double capVolume(int F, const std::vector<std::array<double, 2>>& poly) const {
        if (F != 1) return 0;
        double v = 0;
        for (size_t q = 0; q < poly.size(); ++q) {
            const auto& A = poly[q];
            const auto& B = poly[(q + 1) % poly.size()];
            const double s = 0.5 * (A[0] + B[0]), t = 0.5 * (A[1] + B[1]);
            v += shape_.fluxArea({L_[0], s, t}) * (B[1] - A[1]);
        }
        return v;
    }

    // Caps on box face F: the region inside the lattice, bounded by the patch
    // boundary curves (reversed) and inside parts of the perimeter. A face on
    // a periodic axis is a seam, not a cap.
    bool buildCaps(int F, std::string& error) {
        const int a = F / 2, side = F % 2;
        if (periodic_ >> a & 1) return true;
        const int s = side ? (a + 1) % 3 : (a + 2) % 3, t = side ? (a + 2) % 3 : (a + 1) % 3;  // e_s x e_t = outward
        std::vector<Directed> de;
        for (const Patch& p : patches_)
            for (size_t k = 0; k < 4; ++k) {
                const int x = p.v[k], y = p.v[(k + 1) % 4];
                if (verts_[size_t(x)].faces & verts_[size_t(y)].faces & (1 << F))
                    de.push_back({y, x, p.e[k], edges_[size_t(p.e[k])].a == y});
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
            sampleLine(edges_[size_t(e)]);
            de.push_back({x, y, e, edges_[size_t(e)].a == x});
        }
        return closeCaps(F, de, false, error);
    }

    // Walks directed edges on box face F into loops and adds the cap faces:
    // counter-clockwise loops are outer bounds, each clockwise loop is a hole
    // in the smallest outer loop around it. Loops are measured in the
    // parameter box, or in mm for a flat cap welded across the seam (its
    // parameter rectangle is cut there).
    bool closeCaps(int F, const std::vector<Directed>& de, bool inMm, std::string& error) {
        const int a = F / 2, side = F % 2;
        const int s = side ? (a + 1) % 3 : (a + 2) % 3, t = side ? (a + 2) % 3 : (a + 1) % 3;  // e_s x e_t = outward
        Vec3 nrm, e1{1, 0, 0};
        nrm[a] = side ? 1.0 : -1.0;  // in mm: the tube's ends, normal along z
        const Vec3 e2 = cross(nrm, e1);
        auto uv = [&](const Vec3& P) -> std::array<double, 2> {
            if (inMm) return {dot(P, e1), dot(P, e2)};
            return {P[s], P[t]};
        };

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
                const bool same = de[k].same;
                loop.push_back({de[k].edge, same});
                const std::vector<Vec3> pts = finePoints(e, inMm);
                for (size_t q = 0; q + 1 < pts.size(); ++q) poly.push_back(uv(pts[same ? q : pts.size() - 1 - q]));
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

        auto inside = [&](const std::array<double, 2>& p, const std::vector<std::array<double, 2>>& poly) {
            bool in = false;
            for (size_t q = 0, r = poly.size() - 1; q < poly.size(); r = q++)
                if ((poly[q][1] > p[1]) != (poly[r][1] > p[1]) &&
                    p[0] < (poly[r][0] - poly[q][0]) * (p[1] - poly[q][1]) / (poly[r][1] - poly[q][1]) + poly[q][0])
                    in = !in;
            return in;
        };
        // Volumes need parameter coordinates; the flat caps welded in mm
        // are never on face 1, so they have none.
        std::unordered_map<size_t, size_t> faceOf;
        for (size_t k = 0; k < loops.size(); ++k) {
            if (areas[k] <= 0) continue;
            faceOf[k] = faces_.size();
            Face f;
            f.box = F;
            f.loops.push_back(loops[k]);
            f.volume = inMm ? 0 : capVolume(F, polys[k]);
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
            f.volume += inMm ? 0 : capVolume(F, polys[k]);
        }
        return true;
    }
};

}  // namespace

SmoothReport writeSmoothStep(const std::string& path, const std::string& name, const Lattice& lat,
                             const int cells[3], int patches, double tol, double minFragment, bool exact,
                             double fieldVolume, double shapeVolume) {
    // A coarse patch grid can miss small features: a tunnel inside two cells,
    // or a corner that does not converge. The lattice repeats per cell, so a
    // patch count that is valid for one cell is valid for all. Step it up
    // until the topology is closed and the fit is within tolerance.
    // Folded patches can still close up, so also compare the enclosed volume
    // with the field's.
    SmoothReport r;
    for (int p = patches; p <= std::max(patches + 7, 16); ++p) {
        r = SmoothReport{};
        r.patches = p;
        Builder b(lat, cells, p, tol);
        b.buildPatches();
        for (int n = 4; r.error.empty(); n *= 2) {
            r.maxError = b.run(n, r.error);
            r.samples = n;
            // A large error is folded geometry, which more samples do not fix.
            if (r.maxError <= tol || r.maxError > 20 * tol || n >= 32) break;
        }
        if (r.error.empty() && r.maxError > tol) r.error = "fit error above tolerance";
        std::vector<std::pair<std::vector<int>, double>> shells;
        if (!r.error.empty() || !b.weld(r.error) || !b.shells(shells, r.error)) continue;
        double total = 0, maxVol = 0;
        for (const auto& s : shells) total += s.second, maxVol = std::max(maxVol, s.second);
        if (std::fabs(total - fieldVolume) > 0.02 * shapeVolume) {
            r.error = "enclosed volume does not match the field (folded patches)";
            continue;
        }
        r.patchFaces = 0;
        for (const auto& s : shells) r.patchFaces += s.first.size();
        double mergeErr = 0;
        r.merged = !exact && b.merge(mergeErr, r.warning);
        if (r.merged) {
            r.maxError = std::max(r.maxError, mergeErr);
            shells.clear();
            b.shells(shells, r.error);  // merge() already checked it
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
        // Tolerant edges: merged faces, and edges on curved caps (fitted
        // curves on a cylinder or torus, within tol / 4).
        const double uncertainty = r.merged || lat.shape.bent() ? tol : 1e-7;
        if (kept.empty()) r.error = "nothing to export";
        else if (!writeStepFile(path, name, [&](StepFile& st) { return b.write(st, kept); }, uncertainty))
            r.error = "cannot write '" + path + "'";
        return r;
    }
    return r;
}

}  // namespace tpms
