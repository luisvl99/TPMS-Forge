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
// - Tube, torus shell (conformal): the lattice is built in a parameter box,
//   one unit per cell, and bent into the shape, so cells follow the
//   curvature. Parameter axes: tube = radial, around, along the axis;
//   torus shell = radial, around the ring, around the tube.
// - Cylinder, sphere, torus (trimmed): a straight lattice cut by the shape.
// Round shapes are centered at the origin, axis along z.
enum class ShapeKind { Box, Tube, TorusShell, Cylinder, Sphere, Torus };

struct Shape {
    ShapeKind kind = ShapeKind::Box;
    double r0 = 0, r1 = 0;        // tube, torus shell: inner and outer radius; trimmed: r1 = (tube) radius
    double height = 0;            // tube, cylinder
    double major = 0;             // torus shell, torus: from the axis to the tube center
    bool logRadial = true;        // conformal: r = r0 (r1 / r0)^t keeps cells square; false: linear in r
    double cells[3] = {1, 1, 1};  // conformal: cells along the parameter axes

    bool conformal() const { return kind == ShapeKind::Tube || kind == ShapeKind::TorusShell; }
    bool trimmed() const {
        return kind == ShapeKind::Cylinder || kind == ShapeKind::Sphere || kind == ShapeKind::Torus;
    }
    int periodic() const {  // bit per parameter axis that wraps around
        return kind == ShapeKind::Tube ? 2 : kind == ShapeKind::TorusShell ? 6 : 0;
    }

    // Conformal: parameter point (cell units) to mm. Orientation preserving.
    Vec3 map(Vec3 p) const {
        const double t = p.x / cells[0];
        const double r = logRadial ? r0 * std::pow(r1 / r0, t) : r0 + (r1 - r0) * t;
        const double a = kTwoPi * p.y / cells[1];
        if (kind == ShapeKind::Tube) return {r * std::cos(a), r * std::sin(a), height * (p.z / cells[2] - 0.5)};
        const double b = kTwoPi * p.z / cells[2], d = major + r * std::cos(b);
        return {d * std::cos(a), d * std::sin(a), r * std::sin(b)};
    }

    // Conformal: columns of the Jacobian of map() at p.
    void jacobian(Vec3 p, Vec3 J[3]) const {
        const double h = 1e-6;
        for (int k = 0; k < 3; ++k) {
            Vec3 q = p, r = p;
            q[k] += h, r[k] -= h;
            J[k] = (map(q) - map(r)) * (0.5 / h);
        }
    }

    // Conformal: a field gradient in parameter space as a direction in mm
    // (J^-T g, up to the positive factor 1 / det J).
    Vec3 normal(Vec3 p, Vec3 g) const {
        Vec3 J[3];
        jacobian(p, J);
        return cross(J[1], J[2]) * g.x + cross(J[2], J[0]) * g.y + cross(J[0], J[1]) * g.z;
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

    // Volume in mm^3 (not for the box).
    double volume() const {
        const double pi = 0.5 * kTwoPi;
        switch (kind) {
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
    double level = 0.0;               // iso level c
    double cell[3] = {10, 10, 10};    // unit cell size per axis: mm, or 1 in a conformal parameter box
    Shape shape;                      // trimmed shapes cut the field

    // Surface function f at a point.
    double value(double x, double y, double z) const {
        return f(kTwoPi * x / cell[0], kTwoPi * y / cell[1], kTwoPi * z / cell[2]);
    }

    // Signed field at a point. Negative inside the exported region.
    double operator()(double x, double y, double z) const {
        const double v = value(x, y, z);
        double g = form == Form::Sheet ? std::fabs(v) - level : v - level;
        if (part == Part::Void) g = -g;
        return shape.trimmed() ? std::max(g, shape.inside(x, y, z)) : g;
    }
};

// Iso level c that gives the requested solid volume fraction (0..1) in one unit
// cell. Takes the quantile of f (network) or |f| (sheet) over n^3 cell samples.
inline double levelForDensity(SurfaceFn f, Form form, double density, int n = 100) {
    std::vector<double> vals;
    vals.reserve(size_t(n) * n * n);
    const double h = kTwoPi / n;
    for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                const double v = f((i + 0.5) * h, (j + 0.5) * h, (k + 0.5) * h);
                vals.push_back(form == Form::Sheet ? std::fabs(v) : v);
            }
    const size_t idx = std::min(vals.size() - 1, size_t(density * double(vals.size())));
    std::nth_element(vals.begin(), vals.begin() + std::ptrdiff_t(idx), vals.end());
    return vals[idx];
}

}  // namespace tpms
