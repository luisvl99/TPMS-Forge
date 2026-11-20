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
enum class Part { Solid, Void };     // void = box minus solid (fluid domain)

struct Lattice {
    SurfaceFn f = gyroid;
    Form form = Form::Sheet;
    Part part = Part::Solid;
    double level = 0.0;               // iso level c
    double cell[3] = {10, 10, 10};    // unit cell size per axis, mm

    // Surface function f at a point in mm.
    double value(double x, double y, double z) const {
        return f(kTwoPi * x / cell[0], kTwoPi * y / cell[1], kTwoPi * z / cell[2]);
    }

    // Signed field at a point in mm. Negative inside the exported region.
    double operator()(double x, double y, double z) const {
        const double v = value(x, y, z);
        const double g = form == Form::Sheet ? std::fabs(v) - level : v - level;
        return part == Part::Solid ? g : -g;
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
