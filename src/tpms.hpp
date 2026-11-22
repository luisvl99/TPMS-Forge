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

// Surface function f(x, y, z) with period 2*pi along each axis, and its
// gradient. The surface is f = 0.
using SurfaceFn = double (*)(double x, double y, double z);
using GradientFn = Vec3 (*)(double x, double y, double z);

inline double gyroid(double x, double y, double z) {
    return std::sin(x) * std::cos(y) + std::sin(y) * std::cos(z) + std::sin(z) * std::cos(x);
}
inline Vec3 gyroidGrad(double x, double y, double z) {
    const double sx = std::sin(x), sy = std::sin(y), sz = std::sin(z);
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return {cx * cy - sz * sx, cy * cz - sx * sy, cz * cx - sy * sz};
}

inline double schwarzP(double x, double y, double z) {
    return std::cos(x) + std::cos(y) + std::cos(z);
}
inline Vec3 schwarzPGrad(double x, double y, double z) { return {-std::sin(x), -std::sin(y), -std::sin(z)}; }

inline double diamond(double x, double y, double z) {
    const double sx = std::sin(x), sy = std::sin(y), sz = std::sin(z);
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz;
}
inline Vec3 diamondGrad(double x, double y, double z) {
    const double sx = std::sin(x), sy = std::sin(y), sz = std::sin(z);
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return {cx * sy * sz + cx * cy * cz - sx * sy * cz - sx * cy * sz,
            sx * cy * sz - sx * sy * cz + cx * cy * cz - cx * sy * sz,
            sx * sy * cz - sx * cy * sz - cx * sy * sz + cx * cy * cz};
}

inline double neovius(double x, double y, double z) {
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return 3.0 * (cx + cy + cz) + 4.0 * cx * cy * cz;
}
inline Vec3 neoviusGrad(double x, double y, double z) {
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return {-std::sin(x) * (3 + 4 * cy * cz), -std::sin(y) * (3 + 4 * cz * cx), -std::sin(z) * (3 + 4 * cx * cy)};
}

inline double iwp(double x, double y, double z) {
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return 2.0 * (cx * cy + cy * cz + cz * cx) -
           (std::cos(2 * x) + std::cos(2 * y) + std::cos(2 * z));
}
inline Vec3 iwpGrad(double x, double y, double z) {
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    return {-2 * std::sin(x) * (cy + cz) + 2 * std::sin(2 * x), -2 * std::sin(y) * (cz + cx) + 2 * std::sin(2 * y),
            -2 * std::sin(z) * (cx + cy) + 2 * std::sin(2 * z)};
}

inline double frd(double x, double y, double z) {
    const double c2x = std::cos(2 * x), c2y = std::cos(2 * y), c2z = std::cos(2 * z);
    return 4.0 * std::cos(x) * std::cos(y) * std::cos(z) - (c2x * c2y + c2y * c2z + c2z * c2x);
}
inline Vec3 frdGrad(double x, double y, double z) {
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    const double c2x = std::cos(2 * x), c2y = std::cos(2 * y), c2z = std::cos(2 * z);
    return {-4 * std::sin(x) * cy * cz + 2 * std::sin(2 * x) * (c2y + c2z),
            -4 * std::sin(y) * cz * cx + 2 * std::sin(2 * y) * (c2z + c2x),
            -4 * std::sin(z) * cx * cy + 2 * std::sin(2 * z) * (c2x + c2y)};
}

inline double fischerKochS(double x, double y, double z) {
    return std::cos(2 * x) * std::sin(y) * std::cos(z) +
           std::cos(x) * std::cos(2 * y) * std::sin(z) +
           std::sin(x) * std::cos(y) * std::cos(2 * z);
}
inline Vec3 fischerKochSGrad(double x, double y, double z) {
    const double sx = std::sin(x), sy = std::sin(y), sz = std::sin(z);
    const double cx = std::cos(x), cy = std::cos(y), cz = std::cos(z);
    const double s2x = std::sin(2 * x), s2y = std::sin(2 * y), s2z = std::sin(2 * z);
    const double c2x = std::cos(2 * x), c2y = std::cos(2 * y), c2z = std::cos(2 * z);
    return {-2 * s2x * sy * cz - sx * c2y * sz + cx * cy * c2z, c2x * cy * cz - 2 * cx * s2y * sz - sx * sy * c2z,
            -c2x * sy * sz + cx * c2y * cz - 2 * sx * cy * s2z};
}

struct Surface {
    std::string name;
    std::vector<std::string> aliases;
    std::string formula;
    SurfaceFn f;
    GradientFn grad;
};

inline const std::vector<Surface>& surfaces() {
    static const std::vector<Surface> list = {
        {"gyroid", {"g"}, "sin x cos y + sin y cos z + sin z cos x", gyroid, gyroidGrad},
        {"schwarz-p", {"p", "primitive"}, "cos x + cos y + cos z", schwarzP, schwarzPGrad},
        {"diamond", {"d", "schwarz-d"},
         "sin x sin y sin z + sin x cos y cos z + cos x sin y cos z + cos x cos y sin z", diamond, diamondGrad},
        {"neovius", {}, "3 (cos x + cos y + cos z) + 4 cos x cos y cos z", neovius, neoviusGrad},
        {"iwp", {"schoen-iwp"},
         "2 (cos x cos y + cos y cos z + cos z cos x) - (cos 2x + cos 2y + cos 2z)", iwp, iwpGrad},
        {"frd", {"schoen-frd"},
         "4 cos x cos y cos z - (cos 2x cos 2y + cos 2y cos 2z + cos 2z cos 2x)", frd, frdGrad},
        {"fischer-koch-s", {"fks"},
         "cos 2x sin y cos z + cos x cos 2y sin z + sin x cos y cos 2z", fischerKochS, fischerKochSGrad},
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

// The outer shape of the part.
// - Box: the lattice fills [0, cells x cell].
// - Tube, torus shell, sphere shell: the lattice is built in a parameter
//   box, one unit per cell, and bent into the shape, so cells follow the
//   curvature. Parameter axes: tube = radial, around, along the axis; torus
//   shell = radial, around the ring, around the tube; sphere shell = radial,
//   around (longitude), from south to north. A sphere shell stops short of
//   the poles: cells that follow it cannot reach them (the lattice would
//   differ at one point for each longitude).
// - Box with a graded cell size: also built in a parameter box, and
//   stretched along one axis, where the cell length changes linearly.
// - Cylinder, sphere, torus (trimmed): a straight lattice cut by the shape.
// "Mapped" shapes are the ones built in a parameter box.
// Round shapes are centered at the origin, axis along z.
enum class ShapeKind { Box, Tube, TorusShell, Cylinder, Sphere, Torus, SphereShell };

struct Shape {
    ShapeKind kind = ShapeKind::Box;
    double r0 = 0, r1 = 0;        // tube, torus/sphere shell: inner and outer radius; trimmed: r1 = (tube) radius
    double height = 0;            // tube, cylinder
    double major = 0;             // torus shell, torus: from the axis to the tube center
    bool logRadial = true;        // tube, torus shell: r = r0 (r1 / r0)^t keeps cells square; false: linear in r
    double cells[3] = {1, 1, 1};  // cells along the (parameter) axes
    double cellSize[3] = {10, 10, 10};  // box: cell size per axis, mm
    int cellGrade = -1;                 // box: the axis whose cell length changes, -1: none
    double cellEnd = 0;                 // box: the cell length at the far end of that axis
    int slabAxis = -1;                  // box: plates or buffers on the faces across this axis
    double slab[2] = {0, 0};            // their thickness at the start and the end of it, mm
    // Box, sheet: channel buffers (cross flow). Each fluid channel, A (f > c)
    // and B (f < -c), has an inlet and an outlet buffer on the two faces
    // across its own axis, behind a header plate with holes where that
    // channel opens: the plate closes the other channel there.
    int channelAxis[2] = {-1, -1};  // axis of channel A and of channel B, -1: none
    double header = 0, buffer = 0;  // plate thickness and buffer length, mm
    Vec3 frame[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};  // sphere: axes of map(), poles along frame[2]
    double polar[2] = {0, 0};  // sphere shell: angles from the south pole (-z) where it starts and ends, radians

    bool bent() const {
        return kind == ShapeKind::Tube || kind == ShapeKind::TorusShell || kind == ShapeKind::SphereShell;
    }
    // Sphere shell: the angle from the south pole at parameter w (cells).
    double polarAngle(double w) const { return polar[0] + (polar[1] - polar[0]) * w / cells[2]; }

    // The fluid channel whose buffer lies on box face F (0: A, 1: B), or -1.
    int headerChannel(int F) const { return F / 2 == channelAxis[0] ? 0 : F / 2 == channelAxis[1] ? 1 : -1; }

    // Name of face F = 2 * axis + side of the (parameter) box, where a cut
    // face lies: x-, x+, ... for the box; inner, outer, bottom, top for the
    // tube and torus and sphere shell ("" for the faces where they wrap
    // around).
    const char* faceName(int F) const {
        static const char* const box[6] = {"x-", "x+", "y-", "y+", "z-", "z+"};
        static const char* const bentNames[6] = {"inner", "outer", "", "", "bottom", "top"};
        static const char* const cylinder[6] = {"", "side", "", "", "bottom", "top"};
        if (kind == ShapeKind::Cylinder) return cylinder[F];
        if (kind == ShapeKind::Sphere && r0 > 0) return F < 2 ? bentNames[F] : "";  // hollow: inner, outer
        if (trimmed()) return F == 1 ? "surface" : "";
        if (kind == ShapeKind::TorusShell && F >= 4) return "";
        return bent() ? bentNames[F] : box[F];
    }
    bool mapped() const { return bent() || (kind == ShapeKind::Box && cellGrade >= 0); }
    bool trimmed() const {
        return kind == ShapeKind::Cylinder || kind == ShapeKind::Sphere || kind == ShapeKind::Torus;
    }
    // Bit per parameter axis that wraps around. Trimmed shapes: of the
    // parameter box of their smooth B-rep (see map()).
    int periodic() const {
        return kind == ShapeKind::Tube || kind == ShapeKind::Cylinder || kind == ShapeKind::Sphere ||
                       kind == ShapeKind::SphereShell
                   ? 2
               : kind == ShapeKind::TorusShell || kind == ShapeKind::Torus                        ? 6
                                                                                                  : 0;
    }
    // Trimmed: face F of the parameter box maps to a line or a point (the
    // axis, the center line of a torus, the center or a pole of a sphere).
    bool pole(int F) const {
        if (kind == ShapeKind::Sphere) return (F == 0 && r0 == 0) || F >= 4;  // hollow: no center
        return (kind == ShapeKind::Cylinder || kind == ShapeKind::Torus) && F == 0;
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
    // Trimmed: the parameter box of the smooth B-rep, from the axis, center
    // line or center: cylindrical (radius, angle, z), toroidal (radius from
    // the center line, angle around the ring, angle around the tube) or
    // spherical (radius, longitude, angle from the south pole) coordinates.
    // `cells` steps along each; the straight lattice is evaluated at map(p).
    Vec3 map(Vec3 p) const {
        if (kind == ShapeKind::Box) return {boxPos(0, p.x), boxPos(1, p.y), boxPos(2, p.z)};
        if (trimmed()) {
            const double r = r0 + (r1 - r0) * p.x / cells[0], a = kTwoPi * p.y / cells[1], ca = std::cos(a),
                         sa = std::sin(a);
            if (kind == ShapeKind::Cylinder) return {r * ca, r * sa, height * (p.z / cells[2] - 0.5)};
            if (kind == ShapeKind::Sphere) {
                const double t = 0.5 * kTwoPi * p.z / cells[2], s = r * std::sin(t);
                return frame[0] * (s * ca) + frame[1] * (s * sa) + frame[2] * (-r * std::cos(t));
            }
            const double b = kTwoPi * p.z / cells[2], d = major + r * std::cos(b);
            return {d * ca, d * sa, r * std::sin(b)};
        }
        const double r = radius(p.x), a = kTwoPi * p.y / cells[1];
        if (kind == ShapeKind::Tube) return {r * std::cos(a), r * std::sin(a), height * (p.z / cells[2] - 0.5)};
        if (kind == ShapeKind::SphereShell) {
            const double t = polarAngle(p.z), s = r * std::sin(t);
            return {s * std::cos(a), s * std::sin(a), -r * std::cos(t)};
        }
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
        } else if (kind == ShapeKind::SphereShell) {
            r = norm(q);
            p.z = (std::acos(std::clamp(-q.z / r, -1.0, 1.0)) - polar[0]) / (polar[1] - polar[0]) * cells[2];
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
        const double pi = 0.5 * kTwoPi;
        if (trimmed()) {  // as the tube and torus shell with r0 = 0, and the sphere's r^3 / 3 sin t
            const double r = r0 + (r1 - r0) * p.x / cells[0], k = pi / (cells[1] * cells[2]);
            if (kind == ShapeKind::Cylinder) return k * height * r * r;
            if (kind == ShapeKind::Sphere)
                return k * pi * std::sin(pi * p.z / cells[2]) * (r * r * r - r0 * r0 * r0) * (2.0 / 3);
            return 4 * k * pi * (major * r * r / 2 + std::cos(kTwoPi * p.z / cells[2]) * r * r * r / 3);
        }
        const double r = radius(p.x);
        if (kind == ShapeKind::Tube) return pi * height / (cells[1] * cells[2]) * (r * r - r0 * r0);
        if (kind == ShapeKind::SphereShell)
            return kTwoPi / cells[1] * (polar[1] - polar[0]) / cells[2] * std::sin(polarAngle(p.z)) *
                   (r * r * r - r0 * r0 * r0) / 3;
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
        if (kind == ShapeKind::SphereShell)  // around: narrowest at the end nearest a pole
            return std::min({radial, kTwoPi * r0 * std::min(std::sin(polar[0]), std::sin(polar[1])) / cells[1],
                             r0 * (polar[1] - polar[0]) / cells[2]});
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
        if (kind == ShapeKind::SphereShell) {
            const double t = polarAngle(p.z), st = std::sin(t), ct = std::cos(t);
            J[0] = Vec3{st * ca, st * sa, -ct} * dr;
            J[1] = Vec3{-sa, ca, 0} * (r * st * da);
            J[2] = Vec3{ct * ca, ct * sa, st} * (r * (polar[1] - polar[0]) / cells[2]);
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
        if (kind == ShapeKind::Sphere) return std::max(std::hypot(rho, z) - r1, r0 - std::hypot(rho, z));
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
        case ShapeKind::Sphere: return 4.0 / 3.0 * pi * (r1 * r1 * r1 - r0 * r0 * r0);
        case ShapeKind::Torus: return 2 * pi * pi * major * r1 * r1;
        case ShapeKind::SphereShell:
            return kTwoPi / 3 * (r1 * r1 * r1 - r0 * r0 * r0) * (std::cos(polar[0]) - std::cos(polar[1]));
        default: return 0;
        }
    }
};

struct Lattice {
    SurfaceFn f = gyroid;
    GradientFn grad = gyroidGrad;  // of f
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
    // Sheet of a set wall thickness: value() is f / |grad f| in mm, about the
    // distance to the surface f = 0, and the level is half the wall.
    bool wall = false;
    // Rotated lattice (box and trimmed shapes): value() samples it at
    // turn * p; the rows of turn are those of the inverse rotation.
    bool rotated = false;
    Vec3 turn[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    // Surface function f at a point (or f / |grad f|, see `wall`). The phase
    // (in cells) is reduced to one cell, and points within 1e-12 of a cell
    // boundary are put on it, so a point and its copy one cell over give the
    // same value bit for bit: opposite faces of a cell match exactly.
    double value(double x, double y, double z) const {
        if (rotated) {
            const Vec3 p{x, y, z};
            x = dot(turn[0], p), y = dot(turn[1], p), z = dot(turn[2], p);
        }
        auto phase = [](double u, double offset) {
            const double k = std::round(u);
            if (std::fabs(u - k) < 1e-12) u = k;
            u = u - std::floor(u) + offset;  // offset after the reduction, so copies stay identical
            return kTwoPi * (u - std::floor(u));
        };
        const double u = phase(x / cell[0], shift[0]), v = phase(y / cell[1], shift[1]), w = phase(z / cell[2], shift[2]);
        const double fv = f(u, v, w);
        if (!wall) return fv;
        // The gradient at the reduced phase, so copies stay identical; in a
        // parameter box, mapped to mm.
        Vec3 g = grad(u, v, w);
        for (int a = 0; a < 3; ++a) g[a] *= kTwoPi / cell[a];
        if (shape.mapped()) g = shape.gradient({x, y, z}, g);
        return fv / std::max(norm(g), 1e-12);
    }

    // The field without a trimming shape. A mapped shape stays: value()
    // needs its mapping for a wall in mm.
    Lattice untrimmed() const {
        Lattice l = *this;
        if (l.shape.trimmed()) l.shape = Shape{};
        return l;
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
