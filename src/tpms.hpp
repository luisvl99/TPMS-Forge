#pragma once
// Level-set (nodal) approximations of triply periodic minimal surfaces, and the
// signed field of a lattice built from one of them.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace tpms {

constexpr double kTwoPi = 6.283185307179586;

// Surface function f(x, y, z) with period 2*pi along each axis. The surface is f = 0.
using SurfaceFn = double (*)(double x, double y, double z);

inline double gyroid(double x, double y, double z) {
    return std::sin(x) * std::cos(y) + std::sin(y) * std::cos(z) + std::sin(z) * std::cos(x);
}

inline double schwarzP(double x, double y, double z) {
    return std::cos(x) + std::cos(y) + std::cos(z);
}

inline double diamond(double x, double y, double z) {
    const double sx = std::sin(x), sy = std::sin(y), sz = std::sin(z);
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz;
}

inline double neovius(double x, double y, double z) {
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return 3.0 * (cx + cy + cz) + 4.0 * cx * cy * cz;
}

inline double iwp(double x, double y, double z) {
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return 2.0 * (cx * cy + cy * cz + cz * cx) -
           (std::cos(2 * x) + std::cos(2 * y) + std::cos(2 * z));
}

inline double frd(double x, double y, double z) {
    const double c2x = std::cos(2 * x), c2y = std::cos(2 * y), c2z = std::cos(2 * z);
    return 4.0 * std::cos(x) * std::cos(y) * std::cos(z) - (c2x * c2y + c2y * c2z + c2z * c2x);
}

inline double fischerKochS(double x, double y, double z) {
    return std::cos(2 * x) * std::sin(y) * std::cos(z) +
           std::cos(x) * std::cos(2 * y) * std::sin(z) +
           std::sin(x) * std::cos(y) * std::cos(2 * z);
}

struct Surface {
    std::string name;
    std::vector<std::string> aliases;
    std::string formula;
    SurfaceFn f;
};

inline const std::vector<Surface>& surfaces() {
    static const std::vector<Surface> list = {
        {"gyroid", {"g"}, "sin x cos y + sin y cos z + sin z cos x", gyroid},
        {"schwarz-p", {"p", "primitive"}, "cos x + cos y + cos z", schwarzP},
        {"diamond", {"d", "schwarz-d"},
         "sin x sin y sin z + sin x cos y cos z + cos x sin y cos z + cos x cos y sin z", diamond},
        {"neovius", {}, "3 (cos x + cos y + cos z) + 4 cos x cos y cos z", neovius},
        {"iwp", {"schoen-iwp"},
         "2 (cos x cos y + cos y cos z + cos z cos x) - (cos 2x + cos 2y + cos 2z)", iwp},
        {"frd", {"schoen-frd"},
         "4 cos x cos y cos z - (cos 2x cos 2y + cos 2y cos 2z + cos 2z cos 2x)", frd},
        {"fischer-koch-s", {"fks"},
         "cos 2x sin y cos z + cos x cos 2y sin z + sin x cos y cos 2z", fischerKochS},
    };
    return list;
}

inline const Surface* findSurface(std::string name) {
    for (char& ch : name) ch = char(std::tolower(static_cast<unsigned char>(ch)));
    for (const Surface& s : surfaces()) {
        if (s.name == name) return &s;
        for (const std::string& a : s.aliases)
            if (a == name) return &s;
    }
    return nullptr;
}

enum class Form { Network, Sheet };  // network: f <= c, sheet: |f| <= c
enum class Part { Solid, Void };     // void = shape minus solid (fluid domain)

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

// The outer shape of the part.
// - Box: the lattice fills [0, cells x cell].
// - Tube, torus shell: the lattice is built in a parameter box, one unit per
//   cell, and bent into the shape, so cells follow the curvature. Parameter
//   axes: tube = radial, around, along the axis; torus shell = radial,
//   around the ring, around the tube.
// - Box with a graded cell size: also built in a parameter box, and
//   stretched along one axis, where the cell length changes linearly.
// - Cylinder, sphere, torus (trimmed): a straight lattice cut by the shape.
// "Mapped" shapes are the ones built in a parameter box.
// Round shapes are centered at the origin, axis along z.
enum class ShapeKind { Box, Tube, TorusShell, Cylinder, Sphere, Torus };

struct Shape {
    ShapeKind kind = ShapeKind::Box;
    double r0 = 0, r1 = 0;        // tube, torus shell: inner and outer radius; trimmed: r1 = (tube) radius
    double height = 0;            // tube, cylinder
    double major = 0;             // torus shell, torus: from the axis to the tube center
    bool logRadial = true;        // tube, torus shell: r = r0 (r1 / r0)^t keeps cells square; false: linear in r
    double cells[3] = {1, 1, 1};  // cells along the (parameter) axes
    double cellSize[3] = {10, 10, 10};  // box: cell size per axis, mm
    int cellGrade = -1;                 // box: the axis whose cell length changes, -1: none
    double cellEnd = 0;                 // box: the cell length at the far end of that axis

    bool bent() const { return kind == ShapeKind::Tube || kind == ShapeKind::TorusShell; }
    bool mapped() const { return bent() || (kind == ShapeKind::Box && cellGrade >= 0); }
    bool trimmed() const {
        return kind == ShapeKind::Cylinder || kind == ShapeKind::Sphere || kind == ShapeKind::Torus;
    }
    int periodic() const {  // bit per parameter axis that wraps around
        return kind == ShapeKind::Tube ? 2 : kind == ShapeKind::TorusShell ? 6 : 0;
    }

    // Box: position in mm along axis a at parameter u (cell units), and the
    // local cell length (its derivative), which grows linearly on cellGrade.
    double boxPos(int a, double u) const {
        const double k = a == cellGrade ? (cellEnd - cellSize[a]) / cells[a] : 0;
        return (cellSize[a] + 0.5 * k * u) * u;
    }
    double boxCell(int a, double u) const {
        return a == cellGrade ? cellSize[a] + (cellEnd - cellSize[a]) / cells[a] * u : cellSize[a];
    }

    // Tube, torus shell: radius in mm at radial parameter u (cell units).
    double radius(double u) const {
        const double t = u / cells[0];
        return logRadial ? r0 * std::pow(r1 / r0, t) : r0 + (r1 - r0) * t;
    }

    // Mapped: parameter point (cell units) to mm. Orientation preserving.
    Vec3 map(Vec3 p) const {
        if (kind == ShapeKind::Box) return {boxPos(0, p.x), boxPos(1, p.y), boxPos(2, p.z)};
        const double r = radius(p.x), a = kTwoPi * p.y / cells[1];
        if (kind == ShapeKind::Tube) return {r * std::cos(a), r * std::sin(a), height * (p.z / cells[2] - 0.5)};
        const double b = kTwoPi * p.z / cells[2], d = major + r * std::cos(b);
        return {d * std::cos(a), d * std::sin(a), r * std::sin(b)};
    }

    // Mapped: inverse of map(). Angles come back in [0, 2 pi).
    Vec3 unmap(Vec3 q) const {
        if (kind == ShapeKind::Box) {  // the root of boxPos(a, u) = q[a], in a form stable for k = 0
            Vec3 p;
            for (int a = 0; a < 3; ++a) {
                const double c = cellSize[a], k = a == cellGrade ? (cellEnd - c) / cells[a] : 0;
                p[a] = 2 * q[a] / (c + std::sqrt(std::max(0.0, c * c + 2 * k * q[a])));
            }
            return p;
        }
        auto turn = [](double angle) { return (angle < 0 ? angle + kTwoPi : angle) / kTwoPi; };
        double r;
        Vec3 p;
        p.y = turn(std::atan2(q.y, q.x)) * cells[1];
        if (kind == ShapeKind::Tube) {
            r = std::hypot(q.x, q.y);
            p.z = (q.z / height + 0.5) * cells[2];
        } else {
            const double d = std::hypot(q.x, q.y) - major;
            r = std::hypot(d, q.z);
            p.z = turn(std::atan2(q.z, d)) * cells[2];
        }
        p.x = (logRadial ? std::log(r / r0) / std::log(r1 / r0) : (r - r0) / (r1 - r0)) * cells[0];
        return p;
    }

    // Phi(p) with dPhi/dp.x = det J and Phi = 0 at p.x = 0. The flux
    // (Phi, 0, 0) has divergence det J, so a closed surface in the parameter
    // box encloses the mm volume of the flux through it. Box (not mapped):
    // the parameter box is the part in mm, Phi = p.x.
    double flux(Vec3 p) const {
        if (kind == ShapeKind::Box) return mapped() ? boxPos(0, p.x) * boxCell(1, p.y) * boxCell(2, p.z) : p.x;
        const double pi = 0.5 * kTwoPi, r = radius(p.x);
        if (kind == ShapeKind::Tube) return pi * height / (cells[1] * cells[2]) * (r * r - r0 * r0);
        const double b = kTwoPi * p.z / cells[2];
        return 4 * pi * pi / (cells[1] * cells[2]) *
               (major * (r * r - r0 * r0) / 2 + std::cos(b) * (r * r * r - r0 * r0 * r0) / 3);
    }

    // Psi(p): Phi integrated over p.y from 0, for the flux through caps
    // (Green's theorem). Phi does not depend on p.y except on a graded box.
    double fluxArea(Vec3 p) const {
        if (kind == ShapeKind::Box && mapped()) return boxPos(0, p.x) * boxPos(1, p.y) * boxCell(2, p.z);
        return p.y * flux(p);
    }

    // Mapped: the smallest cell edge in mm (tube and torus shell: radial at
    // r0; around at r0 for the tube and the torus shell's tube, at
    // major - r1 for its ring).
    double minCell() const {
        if (kind == ShapeKind::Box) {
            double m = 1e300;
            for (int a = 0; a < 3; ++a) m = std::min({m, cellSize[a], a == cellGrade ? cellEnd : 1e300});
            return m;
        }
        const double radial = radius(1) - r0;
        if (kind == ShapeKind::Tube) return std::min({radial, kTwoPi * r0 / cells[1], height / cells[2]});
        return std::min({radial, kTwoPi * (major - r1) / cells[1], kTwoPi * r0 / cells[2]});
    }

    // Mapped: columns of the Jacobian of map() at p.
    void jacobian(Vec3 p, Vec3 J[3]) const {
        if (kind == ShapeKind::Box) {
            J[0] = {boxCell(0, p.x), 0, 0}, J[1] = {0, boxCell(1, p.y), 0}, J[2] = {0, 0, boxCell(2, p.z)};
            return;
        }
        const double r = radius(p.x), a = kTwoPi * p.y / cells[1];
        const double dr = (logRadial ? r * std::log(r1 / r0) : r1 - r0) / cells[0];
        const double ca = std::cos(a), sa = std::sin(a), da = kTwoPi / cells[1];
        if (kind == ShapeKind::Tube) {
            J[0] = Vec3{ca, sa, 0} * dr;
            J[1] = Vec3{-sa, ca, 0} * (r * da);
            J[2] = {0, 0, height / cells[2]};
            return;
        }
        const double b = kTwoPi * p.z / cells[2], cb = std::cos(b), sb = std::sin(b);
        J[0] = Vec3{cb * ca, cb * sa, sb} * dr;
        J[1] = Vec3{-sa, ca, 0} * ((major + r * cb) * da);
        J[2] = Vec3{-sb * ca, -sb * sa, cb} * (r * kTwoPi / cells[2]);
    }

    // Mapped: a field gradient g in parameter space as the gradient in mm,
    // J^-T g.
    Vec3 gradient(Vec3 p, Vec3 g) const {
        Vec3 J[3];
        jacobian(p, J);
        const Vec3 c = cross(J[1], J[2]) * g.x + cross(J[2], J[0]) * g.y + cross(J[0], J[1]) * g.z;
        return c * (1.0 / dot(J[0], cross(J[1], J[2])));
    }

    // Trimmed: half extents, and a function that is negative inside.
    Vec3 half() const {
        if (kind == ShapeKind::Cylinder) return {r1, r1, 0.5 * height};
        if (kind == ShapeKind::Sphere) return {r1, r1, r1};
        return {major + r1, major + r1, r1};
    }
    double inside(double x, double y, double z) const {
        const double rho = std::hypot(x, y);
        if (kind == ShapeKind::Cylinder) return std::max(rho - r1, std::fabs(z) - 0.5 * height);
        if (kind == ShapeKind::Sphere) return std::sqrt(rho * rho + z * z) - r1;
        return std::hypot(rho - major, z) - r1;
    }

    // Volume in mm^3.
    double volume() const {
        const double pi = 0.5 * kTwoPi;
        switch (kind) {
        case ShapeKind::Box: return boxPos(0, cells[0]) * boxPos(1, cells[1]) * boxPos(2, cells[2]);
        case ShapeKind::Tube: return pi * (r1 * r1 - r0 * r0) * height;
        case ShapeKind::TorusShell: return 2 * pi * pi * major * (r1 * r1 - r0 * r0);
        case ShapeKind::Cylinder: return pi * r1 * r1 * height;
        case ShapeKind::Sphere: return 4.0 / 3.0 * pi * r1 * r1 * r1;
        case ShapeKind::Torus: return 2 * pi * pi * major * r1 * r1;
        default: return 0;
        }
    }
};

struct Lattice {
    SurfaceFn f = gyroid;
    Form form = Form::Sheet;
    Part part = Part::Solid;
    double level = 0.0;               // iso level c (graded: at the start)
    double cell[3] = {10, 10, 10};    // unit cell size per axis: mm, or 1 in a parameter box
    double shift[3] = {0, 0, 0};      // offset of the lattice, in cells
    Shape shape;                      // trimmed shapes cut the field
    // Graded: c along lattice axis `grade`, from coordinate gradeFrom to
    // gradeTo, as evenly spaced values (linear in between).
    int grade = -1;
    double gradeFrom = 0, gradeTo = 1;
    std::vector<double> levels;

    // Surface function f at a point. The phase (in cells) is reduced to one
    // cell, and points within 1e-12 of a cell boundary are put on it, so a
    // point and its copy one cell over give the same value bit for bit:
    // opposite faces of a cell match exactly.
    double value(double x, double y, double z) const {
        auto phase = [](double u, double offset) {
            const double k = std::round(u);
            if (std::fabs(u - k) < 1e-12) u = k;
            u = u - std::floor(u) + offset;  // offset after the reduction, so copies stay identical
            return kTwoPi * (u - std::floor(u));
        };
        return f(phase(x / cell[0], shift[0]), phase(y / cell[1], shift[1]), phase(z / cell[2], shift[2]));
    }

    // Iso level c at a point.
    double levelAt(double x, double y, double z) const {
        if (grade < 0) return level;
        const double p = grade == 0 ? x : grade == 1 ? y : z;
        const double t = std::clamp((p - gradeFrom) / (gradeTo - gradeFrom), 0.0, 1.0) * double(levels.size() - 1);
        const size_t i = std::min(size_t(t), levels.size() - 2);
        return levels[i] + (levels[i + 1] - levels[i]) * (t - double(i));
    }

    // Signed field at a point. Negative inside the exported region.
    double operator()(double x, double y, double z) const {
        const double v = value(x, y, z), c = levelAt(x, y, z);
        double g = form == Form::Sheet ? std::fabs(v) - c : v - c;
        if (part == Part::Void) g = -g;
        return shape.trimmed() ? std::max(g, shape.inside(x, y, z)) : g;
    }
};

// Sorted values of f (network) or |f| (sheet) at n^3 points of one unit
// cell. The iso level for solid fraction R is their R-quantile.
inline std::vector<double> cellValues(SurfaceFn f, Form form, int n = 100) {
    std::vector<double> vals;
    vals.reserve(size_t(n) * n * n);
    const double h = kTwoPi / n;
    for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                const double v = f((i + 0.5) * h, (j + 0.5) * h, (k + 0.5) * h);
                vals.push_back(form == Form::Sheet ? std::fabs(v) : v);
            }
    std::sort(vals.begin(), vals.end());
    return vals;
}

// Iso level c that gives solid volume fraction R (0..1) in one unit cell.
inline double levelForDensity(const std::vector<double>& cellValues, double density) {
    return cellValues[std::min(cellValues.size() - 1, size_t(density * double(cellValues.size())))];
}

}  // namespace tpms
