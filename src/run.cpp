#include "run.hpp"

#include "export.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <map>

namespace tpms {
namespace {

const char* kUsage = R"(tpms - TPMS lattice to STEP / STL

Usage: tpms [options] -o FILE [-o FILE ...]

Lattice
  --type NAME           Surface type (default gyroid). See --list.
  --form sheet|network  sheet: |f| <= c, a wall around the surface.
                        network: f <= c, one of the two channel networks.
                        Default sheet.
  --density R [R1]      Target solid volume fraction of one unit cell,
                        0 < R < 1 (default 0.2). With R1: graded, from R at
                        the start of the --grade axis to R1 at its end.
  --level C [C1]        Iso level c (or graded from C to C1). Use instead of
                        --density.
  --grade x|y|z         Axis of a graded --density or --level (default x).
                        Tube, torus-shell: x radial, y around, z along the
                        axis or around the tube. Linear in the cell count.
  --part solid|void     solid, or void = shape minus solid (fluid domain for
                        CFD). --density always means the solid fraction.

Shape (default box)
  --shape NAME          box | tube | torus-shell   (cells follow the shape)
                        cylinder | sphere | torus  (straight cells, trimmed)
  --cell L [Ly Lz]      Unit cell size in mm (box and trimmed; default 10).
  --unit-cell           box: a periodic cell for homogenization (RVE).
                        --cells defaults to 1; opposite faces match exactly;
                        edges are exact and every piece is kept (small ones
                        connect to the next cell). No grading.
  --shift S [Sy Sz]     Offset of the lattice in cells (default 0): moves
                        where the faces of the shape cut the cells.
  --plates x|y|z T [T1] box, solid: plates T mm thick on the faces at the
                        start and end of the axis (T1: at the end; 0: none),
                        fused with the lattice. Sandwich panels, walls.
  --buffers x|y|z L [L1] box, void: empty buffer zones L mm long before and
                        after the lattice along the axis, fused with the
                        fluid: inlet and outlet regions for CFD.
  --channel-buffers A B L [T]
                        box, sheet: separate inlets and outlets for the two
                        fluid channels (heat exchanger, cross flow). Channel
                        A (f > c) flows along axis A, channel B (f < -c)
                        along axis B (x, y or z). Each gets buffers L mm long
                        on its two faces, behind header plates T mm thick
                        (default 1) with holes where it opens: the plates
                        close the other channel there. --part void writes
                        each channel's fluid to its own file; --part solid
                        writes the lattice with the plates.
  --cell-grade x|y|z L1 box: the cell length along that axis changes
                        linearly from --cell (start) to L1 mm (end).
  --cells N [Ny Nz]     Cell counts (default 2).
                          box: along x, y, z
                          tube: radial, around, along the axis
                          torus-shell: radial, around the ring, around the tube
  --radius R | R0 R1    tube, torus-shell: inner and outer radius.
                        cylinder, sphere: radius. torus: tube radius.
  --height H            tube, cylinder: length along z.
  --major R             torus-shell, torus: from the axis to the tube center.
  --radial log|linear   tube, torus-shell: radial spacing (default log, which
                        keeps cells square; linear keeps them equally thick).
  Round shapes are centered at the origin with their axis along z.

STEP (smooth B-spline faces; every shape)
  --patches N           Patch grid steps per unit cell edge to start from
                        (default 3). The tool steps it up until the geometry
                        is valid and within --tol.
  --tol T               Max fit error in mm (default 0.001 x smallest cell
                        edge).
  --exact               Keep one face per patch, with exactly shared edges.
                        Default: merge patches into large trimmed faces whose
                        edges lie within --tol of the faces.
  --faceted             Write one planar face per triangle instead (v0).

STL and faceted STEP
  --res N               Grid steps per unit cell edge (default 24).
                        Thin walls need more. The tool warns when too low.

Both
  --min-fragment F      Drop pieces smaller than F x largest piece volume
                        (default 0.01). Use 0 to keep all pieces.

Output
  -o, --out FILE        .step / .stp, .stl (binary), or .json / .csv (the
                        property report). Repeat for several files.
  --report              Print the property report: solid and fluid volumes,
                        wetted area, hydraulic diameters, wall thickness,
                        permeability, stiffness bounds, and which faces of
                        the shape each piece touches.
  --material RHO [E]    Report: solid density in kg/m^3 (mass), and Young's
                        modulus in GPa (stiffness estimates in GPa).
  --flow U [MU RHO]     Report: superficial (inlet) velocity in m/s, fluid
                        viscosity in Pa s and density in kg/m^3 (default
                        water: 0.001, 998). Adds Reynolds numbers and
                        pressure drop estimates per fluid channel.
  --list                List surface types.
  -h, --help            Show this help.

Examples
  tpms --type gyroid --form sheet --density 0.25 --cell 8 --cells 4 -o gyroid.step -o gyroid.stl
  tpms --shape tube --radius 15 25 --height 30 --cells 1 12 3 -o tube.stl)";

struct ShapeName {
    const char* name;
    ShapeKind kind;
};
const ShapeName kShapes[] = {{"box", ShapeKind::Box},           {"tube", ShapeKind::Tube},
                             {"torus-shell", ShapeKind::TorusShell}, {"cylinder", ShapeKind::Cylinder},
                             {"sphere", ShapeKind::Sphere},     {"torus", ShapeKind::Torus}};

// x, y or z as 0, 1, 2; -1 otherwise.
int axisOf(const std::string& s) {
    const std::string a = s;
    return a == "x" ? 0 : a == "y" ? 1 : a == "z" ? 2 : -1;
}

bool number(const std::string& s, double& out) {
    char* end = nullptr;
    out = std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0' && std::isfinite(out);
}

// Reads up to `max` numbers after the flag at args[i]. Returns how many;
// the rest of `out` keeps its values (defaults).
int numbers(const std::vector<std::string>& args, size_t& i, double* out, int max) {
    int n = 0;
    double v;
    while (n < max && i + 1 < args.size() && number(args[i + 1], v)) out[n++] = v, ++i;
    return n;
}

std::string lower(std::string s) {
    for (char& ch : s) ch = char(std::tolower(static_cast<unsigned char>(ch)));
    return s;
}

std::string ext(const std::string& path) { return lower(std::filesystem::path(path).extension().string()); }

void say(const Print& print, const char* fmt, ...) {
    char buf[2048];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    print(buf);
}

// The first axis whose opposite faces of the box [0, L] differ in a mesh,
// or -1: the mesh points on face 0, and on face L moved by -L, must be the
// same set.
int unmatchedFaces(const Mesh& m, const double L[3]) {
    for (int a = 0; a < 3; ++a) {
        const int u = (a + 1) % 3, v = (a + 2) % 3;
        std::vector<std::array<double, 2>> near, far;
        for (const Vec3& p : m.vertices) {
            if (std::fabs(p[a]) < 1e-9 * L[a]) near.push_back({p[u], p[v]});
            else if (std::fabs(p[a] - L[a]) < 1e-9 * L[a]) far.push_back({p[u], p[v]});
        }
        std::sort(near.begin(), near.end()), std::sort(far.begin(), far.end());
        if (near.size() != far.size()) return a;
        for (size_t k = 0; k < near.size(); ++k)
            if (std::fabs(near[k][0] - far[k][0]) + std::fabs(near[k][1] - far[k][1]) > 1e-9 * (L[u] + L[v])) return a;
    }
    return -1;
}

// Name of a mesh triangle (mm) for the faceted STEP: the face of the shape
// it lies on, or a wall (wall A / wall B: the sides f = c and f = -c of a
// sheet).
std::string triangleName(const Options& o, const Lattice& lat, const Grid& g, Vec3 a, Vec3 b, Vec3 c) {
    const Shape& s = o.shape;
    if (s.mapped()) a = s.unmap(a), b = s.unmap(b), c = s.unmap(c);
    double h = 1e300;
    for (int k = 0; k < 3; ++k) h = std::min(h, g.size[k] / g.steps[k]);
    const Vec3 m = (a + b + c) * (1.0 / 3);
    if (s.trimmed()) {  // on the shape where its field, not the lattice's, is zero
        const double v = lat.value(m.x, m.y, m.z), lc = lat.levelAt(m.x, m.y, m.z);
        if (std::fabs(s.inside(m.x, m.y, m.z)) < std::fabs(lat.form == Form::Sheet ? std::fabs(v) - lc : v - lc)) {
            if (s.kind != ShapeKind::Cylinder) return "surface";
            if (std::fabs(std::hypot(m.x, m.y) - s.r1) < std::fabs(std::fabs(m.z) - 0.5 * s.height)) return "side";
            return m.z < 0 ? "bottom" : "top";
        }
    } else {
        // Plates or buffers: their outer faces are the shape's faces there;
        // on the lattice's face behind them lies a wall. These planes can
        // fall on grid nodes, where the mesher moves crossings by 2 % of a
        // step: hence the looser test.
        const int sa = s.slabAxis;
        auto onPlane = [&](int k, double at) {
            return std::fabs(a[k] - at) < 0.05 * h && std::fabs(b[k] - at) < 0.05 * h && std::fabs(c[k] - at) < 0.05 * h;
        };
        // Channel buffers of channel K on face F: on F the plate covers the
        // other channel (its wall; slivers of walls that touch F go by the
        // field, as walls do); beyond F, the plate's outer face and the
        // holes through it are K's walls, the buffer's end is the port (F's
        // name), and its sides are named after the faces they lie on.
        const Grid box = gridOf(o, 1);
        for (int F = 0; F < 6; ++F) {
            const int K = s.headerChannel(F), ax = F / 2;
            if (K < 0) continue;
            bool over = true;  // the centroid lies over the face's rectangle
            for (int k = 0; k < 3; ++k) over &= k == ax || (m[k] > -0.05 * h && m[k] < box.size[k] + 0.05 * h);
            const double face = F % 2 ? box.size[ax] : 0, out = F % 2 ? 1 : -1;
            if (!over) continue;
            if (onPlane(ax, face)) return lat.value(m.x, m.y, m.z) > 0 ? "wall A" : "wall B";
            if ((m[ax] - face) * out <= 0) continue;  // inside the box
            if (onPlane(ax, face + out * s.header)) return K ? "wall B" : "wall A";
            if (onPlane(ax, face + out * (s.header + s.buffer))) return s.faceName(F);
            for (int k = 0; k < 3; ++k)
                for (int side = 0; side < 2 && k != ax; ++side)
                    if (onPlane(k, side * box.size[k])) return s.faceName(2 * k + side);
            return K ? "wall B" : "wall A";
        }
        if (sa >= 0) {
            const double hi = gridOf(o, 1).size[sa];
            if (onPlane(sa, -s.slab[0])) return s.faceName(2 * sa);
            if (onPlane(sa, hi + s.slab[1])) return s.faceName(2 * sa + 1);
            if (onPlane(sa, 0) || onPlane(sa, hi)) {
                if (o.part == Part::Void || lat.form != Form::Sheet) return "wall";
                return lat.value(m.x, m.y, m.z) > 0 ? "wall A" : "wall B";  // the fluid channel it faces
            }
        }
        for (int k = 0; k < 3; ++k)
            for (int side = 0; side < 2 && !(g.periodic >> k & 1); ++side) {
                const double at = g.origin[k] + side * g.size[k];
                if (std::fabs(a[k] - at) < 1e-6 * h && std::fabs(b[k] - at) < 1e-6 * h && std::fabs(c[k] - at) < 1e-6 * h)
                    return s.faceName(2 * k + side);
            }
    }
    if (lat.form != Form::Sheet) return "wall";
    return lat.value(a.x, a.y, a.z) > 0 ? "wall A" : "wall B";  // a vertex lies on the wall: f = +-c
}

}  // namespace

const char* shapeName(ShapeKind k) {
    for (const auto& s : kShapes)
        if (s.kind == k) return s.name;
    return "?";
}

bool parseArgs(const std::vector<std::string>& args, Options& o, std::string& error) {
    double cells[3] = {double(o.cells[0]), double(o.cells[1]), double(o.cells[2])};
    bool cellGiven = false, cellsGiven = false, densityGiven = false;
    int radii = 0;
    for (size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        const bool hasValue = i + 1 < args.size();
        double v, r[2];
        int n;
        if (a == "--type" && hasValue) {
            o.type = args[++i];
        } else if (a == "--form" && hasValue) {
            const std::string f = lower(args[++i]);
            if (f == "sheet") o.form = Form::Sheet;
            else if (f == "network" || f == "skeletal") o.form = Form::Network;
            else return error = "--form must be sheet or network", false;
        } else if (a == "--part" && hasValue) {
            const std::string p = lower(args[++i]);
            if (p == "solid") o.part = Part::Solid;
            else if (p == "void" || p == "fluid") o.part = Part::Void;
            else return error = "--part must be solid or void", false;
        } else if (a == "--density" || a == "--level") {
            n = numbers(args, i, r, 2);
            if (n == 0) return error = a + " needs 1 or 2 numbers", false;
            const bool isLevel = a == "--level";
            (isLevel ? o.level : o.density) = r[0];
            (isLevel ? o.level1 : o.density1) = r[n - 1];
            o.hasLevel |= isLevel, densityGiven |= !isLevel;
            if (o.hasLevel && densityGiven) return error = "use --density or --level, not both", false;
            o.graded = n == 2;
        } else if (a == "--grade" && hasValue) {
            if ((o.grade = axisOf(lower(args[++i]))) < 0) return error = "--grade must be x, y or z", false;
        } else if (a == "--cell-grade" && i + 2 < args.size() && number(args[i + 2], v)) {
            if ((o.shape.cellGrade = axisOf(lower(args[i + 1]))) < 0) return error = "--cell-grade needs x, y or z", false;
            o.shape.cellEnd = v;
            i += 2;
        } else if (a == "--report") {
            o.report = true;
        } else if ((a == "--plates" || a == "--buffers") && hasValue) {
            if ((o.shape.slabAxis = axisOf(lower(args[++i]))) < 0) return error = a + " needs x, y or z", false;
            n = numbers(args, i, o.shape.slab, 2);
            if (n == 0) return error = a + " needs 1 or 2 thicknesses", false;
            if (n == 1) o.shape.slab[1] = o.shape.slab[0];
            (a == "--plates" ? o.plates : o.buffers) = true;
        } else if (a == "--channel-buffers" && i + 2 < args.size()) {
            int* axis = o.shape.channelAxis;
            axis[0] = axisOf(lower(args[i + 1])), axis[1] = axisOf(lower(args[i + 2]));
            if (axis[0] < 0 || axis[1] < 0 || axis[0] == axis[1])
                return error = "--channel-buffers needs two different axes (x, y or z): channel A's, then B's", false;
            i += 2;
            n = numbers(args, i, r, 2);
            if (n == 0) return error = "--channel-buffers needs a buffer length (and a plate thickness)", false;
            o.shape.buffer = r[0], o.shape.header = n == 2 ? r[1] : 1;
        } else if (a == "--material") {
            if (numbers(args, i, o.material, 2) == 0 || !(o.material[0] > 0) || o.material[1] < 0)
                return error = "--material needs a density > 0 in kg/m^3 (and a Young's modulus in GPa)", false;
        } else if (a == "--flow") {
            if (numbers(args, i, o.flow, 3) == 0 || !(o.flow[0] > 0 && o.flow[1] > 0 && o.flow[2] > 0))
                return error = "--flow needs a velocity in m/s (and a viscosity in Pa s and a density in kg/m^3), all > 0",
                       false;
        } else if (a == "--unit-cell") {
            o.unitCell = true;
        } else if (a == "--shift") {
            n = numbers(args, i, o.shift, 3);
            if (n != 1 && n != 3) return error = "--shift needs 1 or 3 numbers", false;
            if (n == 1) o.shift[1] = o.shift[2] = o.shift[0];
        } else if (a == "--shape" && hasValue) {
            const std::string k = lower(args[++i]);
            const auto it = std::find_if(std::begin(kShapes), std::end(kShapes), [&](const ShapeName& s) { return k == s.name; });
            if (it == std::end(kShapes)) return error = "unknown shape '" + k + "'", false;
            o.shape.kind = it->kind;
        } else if (a == "--cell") {
            n = numbers(args, i, o.cell, 3);
            if (n != 1 && n != 3) return error = "--cell needs 1 or 3 numbers", false;
            if (n == 1) o.cell[1] = o.cell[2] = o.cell[0];
            cellGiven = true;
        } else if (a == "--cells") {
            n = numbers(args, i, cells, 3);
            if (n != 1 && n != 3) return error = "--cells needs 1 or 3 numbers", false;
            if (n == 1) cells[1] = cells[2] = cells[0];
            cellsGiven = true;
        } else if (a == "--radius") {
            radii = numbers(args, i, r, 2);
            if (radii == 0) return error = "--radius needs 1 or 2 numbers", false;
            o.shape.r0 = radii == 2 ? r[0] : 0, o.shape.r1 = r[radii - 1];
        } else if (a == "--height" && hasValue && number(args[i + 1], v)) {
            o.shape.height = v;
            ++i;
        } else if (a == "--major" && hasValue && number(args[i + 1], v)) {
            o.shape.major = v;
            ++i;
        } else if (a == "--radial" && hasValue) {
            const std::string s = lower(args[++i]);
            if (s == "log" || s == "linear") o.shape.logRadial = s == "log";
            else return error = "--radial must be log or linear", false;
        } else if (a == "--res" && hasValue && number(args[i + 1], v)) {
            o.res = int(std::clamp(v, 0.0, 1e4));  // int() of a huge value is undefined; range checked below
            ++i;
        } else if (a == "--patches" && hasValue && number(args[i + 1], v)) {
            o.patches = int(std::clamp(v, 0.0, 1e4));
            ++i;
        } else if (a == "--tol" && hasValue && number(args[i + 1], v)) {
            o.tol = v;
            ++i;
        } else if (a == "--exact") {
            o.exact = true;
        } else if (a == "--faceted") {
            o.faceted = true;
        } else if (a == "--min-fragment" && hasValue && number(args[i + 1], v)) {
            o.minFragment = v;
            ++i;
        } else if ((a == "-o" || a == "--out") && hasValue) {
            o.outputs.push_back(args[++i]);
        } else {
            return error = "unknown or incomplete option '" + a + "'", false;
        }
    }

    if (!findSurface(o.type)) return error = "unknown surface type '" + o.type + "'. Run 'tpms --list'.", false;
    if (!o.graded) o.density1 = o.density, o.level1 = o.level;
    if (!o.hasLevel && !(o.density > 0 && o.density < 1 && o.density1 > 0 && o.density1 < 1))
        return error = "--density must be in (0, 1)", false;
    if (o.hasLevel && o.form == Form::Sheet && !(o.level > 0 && o.level1 > 0))
        return error = "sheet needs --level > 0", false;

    Shape& s = o.shape;
    const std::string name = shapeName(s.kind);
    if (s.bent() && cellGiven) return error = "--cell does not apply to --shape " + name + ": set --cells", false;
    if (s.cellGrade >= 0 && s.kind != ShapeKind::Box) return error = "--cell-grade works for the box only", false;
    if (o.plates || o.buffers) {
        const std::string opt = o.plates ? "--plates" : "--buffers";
        if (o.plates && o.buffers) return error = "use --plates or --buffers, not both", false;
        if (s.kind != ShapeKind::Box || s.cellGrade >= 0 || o.unitCell)
            return error = opt + " works for the box (without --cell-grade or --unit-cell)", false;
        if (o.plates && o.part != Part::Solid) return error = "--plates need --part solid (for the fluid: --buffers)", false;
        if (o.buffers && o.part != Part::Void) return error = "--buffers need --part void (for the solid: --plates)", false;
        if (s.slab[0] < 0 || s.slab[1] < 0 || s.slab[0] + s.slab[1] <= 0)
            return error = opt + " needs thicknesses >= 0, not both 0", false;
    }
    if (s.channelAxis[0] >= 0) {
        if (s.kind != ShapeKind::Box || s.cellGrade >= 0 || o.unitCell || o.plates || o.buffers)
            return error = "--channel-buffers works for the box (without --cell-grade, --unit-cell, --plates or --buffers)",
                   false;
        if (o.form != Form::Sheet) return error = "--channel-buffers needs --form sheet (two fluid channels)", false;
        if (!(s.buffer > 0 && s.header > 0)) return error = "--channel-buffers needs a length and a thickness > 0", false;
    }
    if (o.unitCell) {
        if (s.kind != ShapeKind::Box) return error = "--unit-cell works for the box only", false;
        if (o.graded || s.cellGrade >= 0) return error = "--unit-cell needs a lattice without grading", false;
        if (!cellsGiven) cells[0] = cells[1] = cells[2] = 1;
        o.exact = true;      // merged faces on opposite sides would differ within --tol
        o.minFragment = 0;   // a piece cut off here continues in the next cell
    }
    if (s.cellGrade >= 0 && !(s.cellEnd > 0)) return error = "--cell-grade needs a cell length > 0", false;
    if (s.trimmed() && cellsGiven) return error = "--cells does not apply to --shape " + name + ": set --cell", false;
    if (s.kind == ShapeKind::Box && (radii || s.height || s.major))
        return error = "--radius, --height and --major need a --shape", false;
    const bool shell = s.bent(), needsHeight = s.kind == ShapeKind::Tube || s.kind == ShapeKind::Cylinder;
    const bool needsMajor = s.kind == ShapeKind::TorusShell || s.kind == ShapeKind::Torus;
    if (s.kind != ShapeKind::Box) {
        if (shell && !(radii == 2 && s.r0 > 0 && s.r1 > s.r0))
            return error = "--shape " + name + " needs --radius R0 R1 with 0 < R0 < R1", false;
        if (!shell && !(radii == 1 && s.r1 > 0)) return error = "--shape " + name + " needs --radius R > 0", false;
        if (needsHeight && !(s.height > 0)) return error = "--shape " + name + " needs --height > 0", false;
        if (!needsHeight && s.height) return error = "--height does not apply to --shape " + name, false;
        if (needsMajor && !(s.major > s.r1)) return error = "--shape " + name + " needs --major > the radius", false;
        if (!needsMajor && s.major) return error = "--major does not apply to --shape " + name, false;
    }
    for (int a = 0; a < 3; ++a) {
        if (!(o.cell[a] > 0)) return error = "--cell must be > 0", false;
        if (cells[a] < 1 || cells[a] > 1e5 || cells[a] != std::floor(cells[a]))
            return error = "--cells must be whole numbers from 1 to 100000", false;
        o.cells[a] = int(cells[a]);
        s.cells[a] = cells[a];
        s.cellSize[a] = o.cell[a];
    }
    if (o.res < 4 || o.res > 1000) return error = "--res must be from 4 to 1000", false;
    if (o.patches < 2 || o.patches > 1000) return error = "--patches must be from 2 to 1000", false;
    if (o.tol < 0) return error = "--tol must be > 0", false;
    if (o.tol == 0) o.tol = autoTol(o);
    for (const auto& out : o.outputs) {
        const std::string e = ext(out);
        if (e != ".step" && e != ".stp" && e != ".stl" && e != ".json" && e != ".csv")
            return error = "output '" + out + "' must end in .step, .stp, .stl, .json or .csv", false;
    }
    return true;
}

double autoTol(const Options& o) {
    return 0.001 * (o.shape.mapped() ? o.shape.minCell() : std::min({o.cell[0], o.cell[1], o.cell[2]}));
}

Lattice makeLattice(const Options& o) {
    const Surface* surface = findSurface(o.type);
    Lattice lat;
    lat.f = surface->f;
    lat.form = o.form;
    lat.part = o.part;
    lat.shape = o.shape;
    for (int a = 0; a < 3; ++a) lat.cell[a] = o.shape.mapped() ? 1 : o.cell[a], lat.shift[a] = o.shift[a];
    // Iso level at t = 0 (start of a grading) to 1 (end): linear in the
    // density (through the cell's quantiles), or in the level itself.
    const std::vector<double> values = o.hasLevel ? std::vector<double>{} : cellValues(surface->f, o.form);
    auto levelAt = [&](double t) {
        if (o.hasLevel) return o.level + (o.level1 - o.level) * t;
        return levelForDensity(values, o.density + (o.density1 - o.density) * t);
    };
    lat.level = levelAt(0);
    if (o.graded) {
        const int a = lat.grade = o.grade;
        if (o.shape.trimmed()) lat.gradeFrom = -o.shape.half()[a], lat.gradeTo = o.shape.half()[a];
        else lat.gradeFrom = 0, lat.gradeTo = gridOf(o, 1).size[a];  // mm, or cells in a parameter box
        for (int k = 0; k <= 64; ++k) lat.levels.push_back(levelAt(k / 64.0));
    }
    return lat;
}

Grid gridOf(const Options& o, int res) {
    Grid g;
    const Shape& s = o.shape;
    for (int a = 0; a < 3; ++a) {
        if (s.trimmed()) {
            const double h = o.cell[a] / res;
            g.steps[a] = int(std::ceil(2 * s.half()[a] / h)) + 2;  // one step of margin on each side
            g.size[a] = g.steps[a] * h;
            g.origin[a] = -0.5 * g.size[a];
        } else {
            g.steps[a] = o.cells[a] * res;
            g.size[a] = s.mapped() ? o.cells[a] : o.cells[a] * o.cell[a];
        }
    }
    g.periodic = s.trimmed() ? 0 : s.periodic();  // a trimmed shape's mesh grid is a plain box in mm
    return g;
}

Grid meshGrid(const Options& o, int res) {
    Grid g = gridOf(o, res);
    const Shape& s = o.shape;
    // Whole steps added before (F even) or after face F along its axis.
    auto extend = [&](int F, double length) {
        const int a = F / 2;
        const double h = g.size[a] / g.steps[a];
        const int n = int(std::ceil(length / h - 1e-9));
        if (F % 2 == 0) g.origin[a] -= n * h;
        g.steps[a] += n, g.size[a] = g.steps[a] * h;
    };
    if (s.slabAxis >= 0) extend(2 * s.slabAxis, s.slab[0]), extend(2 * s.slabAxis + 1, s.slab[1]);
    for (int F = 0; F < 6; ++F)
        if (s.headerChannel(F) >= 0) extend(F, s.header + (o.part == Part::Void ? s.buffer : 0));
    return g;
}

double slabVolume(const Options& o, const Lattice& lat) {
    const Shape& s = o.shape;
    double v = 0;
    if (s.slabAxis >= 0) v += (s.slab[0] + s.slab[1]) * s.volume() / s.boxPos(s.slabAxis, o.cells[s.slabAxis]);
    const Grid box = gridOf(o, 1);
    for (int F = 0; F < 6; ++F) {
        const int K = s.headerChannel(F), a = F / 2, u = (a + 1) % 3, w = (a + 2) % 3;
        if (K < 0) continue;
        // Where channel K opens on the face (the holes through its plate):
        // quasi-random R2 points, as sampledFraction does in 3D.
        const double g = 1.32471795724474602596, step[2] = {1 / g, 1 / (g * g)};
        const int samples = 1 << 14;
        int open = 0;
        for (int k = 0; k < samples; ++k) {
            Vec3 p;
            p[a] = F % 2 ? box.size[a] : 0;
            p[u] = std::fmod(0.5 + k * step[0], 1.0) * box.size[u], p[w] = std::fmod(0.5 + k * step[1], 1.0) * box.size[w];
            const double f = lat.value(p.x, p.y, p.z), c = lat.levelAt(p.x, p.y, p.z);
            open += K ? f < -c : f > c;
        }
        const double area = box.size[u] * box.size[w], holes = area * open / samples;
        v += o.part == Part::Solid ? s.header * (area - holes) : s.header * holes + s.buffer * area;
    }
    return v;
}

namespace {

// Field of a box with channel buffers (negative inside): `channel` -1 for
// the whole part, or one fluid channel (0: A, 1: B) of a void part. Beyond
// each face across a channel's axis lies a header plate, solid except where
// that channel opens on the face (the face's pattern, extruded), and for
// the void part the buffer beyond it. On the face itself the lattice takes
// the plate's cover, so the plate's face over the other channel lies
// exactly there. kSteep turns mm into field units: large, so that change
// is limited to the face. Lateral bounds are needed only along the other
// channel's axis: elsewhere the grid ends at the box.
double headerField(const Options& o, const Lattice& lat, const double L[3], Vec3 p, int channel) {
    const Shape& s = o.shape;
    const bool solid = o.part == Part::Solid;
    constexpr double kSteep = 1e4;
    auto plate = [&](int K, Vec3 q) {  // < 0 where channel K's plate is solid at a point q of its face
        const double f = lat.value(q.x, q.y, q.z), c = lat.levelAt(q.x, q.y, q.z);
        return K ? -f - c : f - c;
    };
    int out = -1, outside = 0;
    double beyond = 0;
    for (int a = 0; a < 3; ++a) {
        const double d = std::max(-p[a], p[a] - L[a]);
        if (d > 0) out = a, ++outside;
        beyond = std::max(beyond, d);
    }
    if (outside == 0) {
        double g;
        if (channel < 0) {
            g = lat(p.x, p.y, p.z);
        } else {
            const double f = lat.value(p.x, p.y, p.z), c = lat.levelAt(p.x, p.y, p.z);
            g = channel ? f + c : c - f;
        }
        for (int F = 0; F < 6; ++F) {
            const int K = s.headerChannel(F), a = F / 2;
            if (K < 0) continue;
            Vec3 q = p;
            q[a] = F % 2 ? L[a] : 0;
            const double depth = std::fabs(p[a] - q[a]);
            g = solid ? std::min(g, std::max(plate(K, q), kSteep * depth))
                      : std::max(g, std::min(-plate(K, q), -kSteep * depth));
        }
        return g;
    }
    const int a = out, F = 2 * a + int(p[a] > L[a]), K = outside == 1 ? s.headerChannel(F) : -1;
    if (K < 0 || (channel >= 0 && K != channel)) return beyond;
    Vec3 q = p;
    q[a] = F % 2 ? L[a] : 0;
    const double t = std::fabs(p[a] - q[a]);  // beyond the face, mm
    double g = solid ? std::max(plate(K, q), t - s.header)
                     : std::max(std::min(-plate(K, q), kSteep * (s.header - t)), t - s.header - s.buffer);
    for (int b = 0; b < 3; ++b)
        if (b != a && s.headerChannel(2 * b) >= 0) g = std::max(g, std::max(-p[b], p[b] - L[b]));
    return g;
}

}  // namespace

Mesh meshSolid(const Options& o, const Lattice& lat, int res, size_t* droppedDegenerate) {
    // Plates or buffers: the lattice in its box, joined with a slab on each
    // side of the axis that has one. As fields: inside the outer ends, either
    // the lattice or the slabs (cell coordinate p, box from 0 to hi).
    // Written so that the lattice's own region stays inside on the faces
    // where it meets a slab, which a min of two box fields would not do.
    const Shape& s = o.shape;
    const int a = s.slabAxis;
    const bool headers = s.channelAxis[0] >= 0;
    const Grid box = gridOf(o, 1);
    const double hi = a < 0 ? 0 : box.size[a], none = 1e300;
    const Field field = [&](double x, double y, double z) {
        if (headers) return headerField(o, lat, box.size, {x, y, z}, -1);
        const double g = lat(x, y, z);
        if (a < 0) return g;
        const double p = a == 0 ? x : a == 1 ? y : z;
        const double ends = std::max(s.slab[0] > 0 ? -s.slab[0] - p : -none, s.slab[1] > 0 ? p - hi - s.slab[1] : -none);
        const double inBox = std::min(s.slab[0] > 0 ? p : none, s.slab[1] > 0 ? hi - p : none);
        return std::max(ends, std::min(g, inBox));
    };
    Mesh m;
    if (o.part == Part::Void && lat.form == Form::Sheet && a < 0) {
        // A sheet's fluid: each channel (f > c, f < -c) meshed alone. In one
        // mesh of |f| > c, a sheet thinner than the grid leaks between two
        // nodes and joins them. (Buffers join them anyway.)
        size_t dropped[2] = {0, 0};
        for (int c = 0; c < 2; ++c) {
            const double sign = c ? 1 : -1;
            const Mesh part = meshField(
                [&](double x, double y, double z) {
                    if (headers) return headerField(o, lat, box.size, {x, y, z}, c);
                    const double g = sign * lat.value(x, y, z) + lat.levelAt(x, y, z);
                    return s.trimmed() ? std::max(g, s.inside(x, y, z)) : g;
                },
                meshGrid(o, res), &dropped[c]);
            const auto base = uint32_t(m.vertices.size());
            m.vertices.insert(m.vertices.end(), part.vertices.begin(), part.vertices.end());
            for (auto t : part.triangles) m.triangles.push_back({t[0] + base, t[1] + base, t[2] + base});
        }
        if (droppedDegenerate) *droppedDegenerate = dropped[0] + dropped[1];
    } else {
        m = meshField(field, meshGrid(o, res), droppedDegenerate);
    }
    if (o.shape.mapped())
        for (Vec3& v : m.vertices) v = o.shape.map(v);
    return m;
}

double shapeVolume(const Options& o) { return o.shape.volume(); }

double sampledFraction(const Options& o, const Lattice& lat) {
    // Quasi-random R3 points: a regular grid can line up with the lattice.
    // Points in a parameter box are weighted by the volume the mapping gives them.
    const Grid g = gridOf(o, 1);
    const double phi = 1.2207440846057595, step[3] = {1 / phi, 1 / (phi * phi), 1 / (phi * phi * phi)};
    double in = 0, all = 0;
    for (int k = 0; k < (1 << 17); ++k) {
        Vec3 p;
        for (int a = 0; a < 3; ++a) p[a] = g.origin[a] + std::fmod(0.5 + k * step[a], 1.0) * g.size[a];
        double w = 1;
        if (o.shape.mapped()) {
            Vec3 J[3];
            o.shape.jacobian(p, J);
            w = dot(J[0], cross(J[1], J[2]));
        } else if (o.shape.trimmed() && o.shape.inside(p.x, p.y, p.z) >= 0) {
            continue;
        }
        all += w;
        if (lat(p.x, p.y, p.z) < 0) in += w;
    }
    return all > 0 ? in / all : 0;
}

Pieces keepPieces(const Mesh& mesh, double minFragment) {
    Pieces out;
    std::vector<Mesh> pieces = splitComponents(mesh);
    std::vector<double> vols;
    double maxVol = 0;
    for (const Mesh& p : pieces) {
        vols.push_back(volume(p));
        maxVol = std::max(maxVol, vols.back());
    }
    std::vector<int> keptAs(pieces.size(), -1);
    for (size_t i = 0; i < pieces.size(); ++i) {
        if (vols[i] <= 0) continue;
        if (vols[i] < minFragment * maxVol) ++out.fragments;
        else keptAs[i] = int(out.kept.size()), out.kept.push_back(pieces[i]);
    }
    // An inner cavity joins the smallest piece around it (ray parity from a
    // point on it). The cavity of a dropped fragment goes with it; a tiny
    // cavity is a fragment too, and is filled.
    for (size_t i = 0; i < pieces.size(); ++i) {
        if (vols[i] > 0) continue;
        if (-vols[i] < minFragment * maxVol) {
            ++out.fragments;
            continue;
        }
        int host = -1, outer = 0;
        for (size_t j = 0; j < pieces.size(); ++j)
            if (vols[j] > 0) host = int(j), ++outer;
        if (outer > 1) host = -1;
        // Several pieces: a ray from a few of its triangles (one along a grid
        // line can graze edges and miscount).
        for (size_t k = 0; k < pieces[i].triangles.size() && k < 5 && host < 0; ++k) {
            const auto& t = pieces[i].triangles[k * pieces[i].triangles.size() / 5];
            const Vec3 q = (pieces[i].vertices[t[0]] + pieces[i].vertices[t[1]] + pieces[i].vertices[t[2]]) * (1.0 / 3);
            for (size_t j = 0; j < pieces.size(); ++j) {
                if (vols[j] <= 0 || (host >= 0 && vols[j] >= vols[size_t(host)])) continue;
                int hits = 0;
                for (const auto& u : pieces[j].triangles)
                    hits += rayHits(pieces[j].vertices[u[0]], pieces[j].vertices[u[1]], pieces[j].vertices[u[2]], q);
                if (hits % 2) host = int(j);
            }
        }
        if (host < 0) {
            ++out.cavities;
            continue;
        }
        if (keptAs[size_t(host)] < 0) continue;
        Mesh& m = out.kept[size_t(keptAs[size_t(host)])];
        const auto base = uint32_t(m.vertices.size());
        m.vertices.insert(m.vertices.end(), pieces[i].vertices.begin(), pieces[i].vertices.end());
        for (auto tri : pieces[i].triangles) m.triangles.push_back({tri[0] + base, tri[1] + base, tri[2] + base});
    }
    for (size_t i = 0; i < out.kept.size(); ++i) out.volume += volume(out.kept[i]), out.triangles += out.kept[i].triangles.size();
    return out;
}

int run(const std::vector<std::string>& args, const Print& print) {
    for (const auto& a : args) {
        if (a == "-h" || a == "--help") {
            print(kUsage);
            return 0;
        }
        if (a == "--list") {
            for (const Surface& s : surfaces()) {
                std::string names = s.name;
                for (const auto& al : s.aliases) names += ", " + al;
                say(print, "  %-28s f = %s", names.c_str(), s.formula.c_str());
            }
            return 0;
        }
    }
    Options o;
    std::string error;
    if (!parseArgs(args, o, error)) return print("error: " + error), 2;
    if (o.outputs.empty() && !o.report) return print("error: no output file. Add -o lattice.step (or --report)"), 2;
    bool needMesh = false, needReport = o.report;
    auto isReport = [](const std::string& out) { return ext(out) == ".json" || ext(out) == ".csv"; };
    for (const auto& out : o.outputs) {
        needMesh |= ext(out) == ".stl" || (o.faceted && !isReport(out));
        needReport |= isReport(out);
    }

    const Lattice lat = makeLattice(o);
    const int* cells = o.cells;
    const Shape& s = o.shape;
    const double shapeVol = shapeVolume(o), slabVol = slabVolume(o, lat);
    auto solidFraction = [&](double exported) { return o.part == Part::Solid ? exported : 1.0 - exported; };
    // Solid fraction of a volume of the exported part: in the lattice box,
    // without plates or buffers.
    auto boxFraction = [&](double vol) { return solidFraction((vol - slabVol) / shapeVol); };
    // What the field gives, from samples. Bending or grading a lattice moves
    // it away from --density, so compare meshes with this, not the target.
    const double expected = sampledFraction(o, lat);
    auto checkFraction = [&](double vol, const char* fix) {
        vol -= slabVol;
        if (std::fabs(vol / shapeVol - expected) > 0.05 * expected)
            say(print, "warning: solid fraction %.4f is far from the field's %.4f. %s",
                solidFraction(vol / shapeVol), solidFraction(expected), fix);
    };
    auto dropped = [&](size_t fragments, size_t cavities) {
        if (fragments)
            say(print, "note: dropped %zu small fragment(s) or filled tiny cavities (< %g x largest).", fragments,
                o.minFragment);
        if (cavities) say(print, "warning: dropped %zu inner cavity shell(s) found inside no solid.", cavities);
    };
    auto wrote = [&](const std::string& path) {
        std::error_code ec;
        const auto bytes = std::filesystem::file_size(path, ec);
        say(print, "Wrote     : %s (%.1f MB)", path.c_str(), ec ? 0.0 : double(bytes) / 1e6);
    };

    say(print, "Surface   : %s, %s, %s part", findSurface(o.type)->name.c_str(),
        o.form == Form::Sheet ? "sheet" : "network", o.part == Part::Solid ? "solid" : "void");
    const char* const kAxes = "xyz";
    if (o.graded && o.hasLevel)
        say(print, "Level c   : %.5g to %.5g along %c", lat.levels.front(), lat.levels.back(), kAxes[o.grade]);
    else if (o.graded)
        say(print, "Level c   : %.5g to %.5g along %c (from target solid fraction %.3f to %.3f)", lat.levels.front(),
            lat.levels.back(), kAxes[o.grade], o.density, o.density1);
    else if (o.hasLevel)
        say(print, "Level c   : %.5g", lat.level);
    else
        say(print, "Level c   : %.5g (from target solid fraction %.3f)", lat.level, o.density);
    switch (s.kind) {
    case ShapeKind::Box:
        say(print, "Box       : %g x %g x %g mm (%d x %d x %d cells of %g x %g x %g mm)", s.boxPos(0, cells[0]),
            s.boxPos(1, cells[1]), s.boxPos(2, cells[2]), cells[0], cells[1], cells[2], o.cell[0], o.cell[1], o.cell[2]);
        if (o.unitCell) say(print, "Unit cell : periodic; exact edges, every piece kept");
        if (s.slabAxis >= 0)
            say(print, "%s   : %g mm at %s, %g mm at %s (%.6g mm^3, not in the fractions)",
                o.plates ? "Plates " : "Buffers", s.slab[0], s.faceName(2 * s.slabAxis), s.slab[1],
                s.faceName(2 * s.slabAxis + 1), slabVol);
        if (s.channelAxis[0] >= 0)
            say(print, "Channels  : A along %c, B along %c; buffers %g mm behind header plates %g mm thick "
                       "(%.6g mm^3, not in the fractions)",
                kAxes[s.channelAxis[0]], kAxes[s.channelAxis[1]], s.buffer, s.header, slabVol);
        if (s.cellGrade >= 0)
            say(print, "Cell size : along %c from %g to %g mm (graded)", kAxes[s.cellGrade], o.cell[s.cellGrade],
                s.cellEnd);
        break;
    case ShapeKind::Tube:
    case ShapeKind::TorusShell:
        say(print, "Shape     : %s, radius %g to %g mm, %s %g mm, %s radial spacing", shapeName(s.kind), s.r0, s.r1,
            s.kind == ShapeKind::Tube ? "height" : "major radius", s.kind == ShapeKind::Tube ? s.height : s.major,
            s.logRadial ? "log" : "linear");
        say(print, "Cells     : %d radial x %d around x %d %s", cells[0], cells[1], cells[2],
            s.kind == ShapeKind::Tube ? "along the axis" : "around the tube");
        break;
    default: {
        char extra[64] = "";
        if (s.kind == ShapeKind::Cylinder) std::snprintf(extra, sizeof extra, ", height %g mm", s.height);
        if (s.kind == ShapeKind::Torus) std::snprintf(extra, sizeof extra, ", major radius %g mm", s.major);
        say(print, "Shape     : %s, radius %g mm%s, trimmed from cells of %g x %g x %g mm", shapeName(s.kind), s.r1,
            extra, o.cell[0], o.cell[1], o.cell[2]);
    }
    }
    if (o.graded)
        say(print, "note: the field gives an overall solid fraction of %.4f.", solidFraction(expected));
    else if (std::fabs(solidFraction(expected) - o.density) > 0.01 && !o.hasLevel)
        say(print, "note: the field gives solid fraction %.4f here, not %.3f: the shape bends or cuts the cells.",
            solidFraction(expected), o.density);

    if (o.shift[0] || o.shift[1] || o.shift[2])
        say(print, "Shift     : %g x %g x %g cells", o.shift[0], o.shift[1], o.shift[2]);
    const auto t0 = std::chrono::steady_clock::now();
    bool ok = true;
    if (o.unitCell) {  // self-check: opposite faces must match
        const double L[3] = {cells[0] * o.cell[0], cells[1] * o.cell[1], cells[2] * o.cell[2]};
        const int bad = unmatchedFaces(meshSolid(o, lat, o.res), L);
        if (bad < 0) say(print, "Periodic  : opposite faces match");
        else say(print, "error: the faces across %c do not match", kAxes[bad]), ok = false;
    }
    Pieces mesh;
    if (needMesh) {
        size_t degenerate = 0;
        const Mesh full = meshSolid(o, lat, o.res, &degenerate);
        const ClosedCheck check = checkClosed(full);
        mesh = keepPieces(full, o.minFragment);
        if (mesh.kept.empty()) return print("error: nothing to export. Check --density / --level."), 2;
        say(print, "Mesh      : %d steps per cell, %zu triangles, %zu solid(s), closed: %s, solid fraction %.4f", o.res,
            mesh.triangles, mesh.kept.size(), check.closed() ? "yes" : "NO", boxFraction(mesh.volume));
        if (!check.closed()) {
            say(print, "warning: mesh not closed (%zu open, %zu non-manifold edges). Try another --res.",
                check.openEdges, check.nonManifoldEdges);
            ok = false;
        }
        if (degenerate) say(print, "warning: dropped %zu zero-area triangles.", degenerate);
        dropped(mesh.fragments, mesh.cavities);
        checkFraction(mesh.volume, "Raise --res.");
        if (mesh.triangles > 300000) say(print, "note: %zu triangles is heavy for CAD. Lower --res.", mesh.triangles);
    }

    // A sheet's fluid is two channels: each goes to a STEP file of its own.
    const bool splitChannels = o.part == Part::Void && o.form == Form::Sheet;
    auto channelNote = [](bool both) {
        return std::string(both ? "a fluid solid lies in both channels (they meet in a buffer, or the sheet has holes "
                                  "at this --res), so they stay in one file"
                                : "only one fluid channel is left, so it is one file");
    };
    for (const auto& out : o.outputs) {
        const std::filesystem::path p(out);
        if (isReport(out)) continue;  // the report, below
        if (ext(out) == ".stl" || o.faceted) {
            // Faceted STEP: solids named by part and fluid channel, faces as
            // walls or by the face of the shape they lie on.
            const Grid grid = meshGrid(o, o.res);
            auto name = [&](size_t k, size_t t) {
                const Mesh& m = mesh.kept[k];
                const auto& tri = m.triangles[t];
                return triangleName(o, lat, grid, m.vertices[tri[0]], m.vertices[tri[1]], m.vertices[tri[2]]);
            };
            std::vector<std::string> solidNames;
            std::map<std::string, int> seen;
            for (size_t k = 0; k < mesh.kept.size() && ext(out) != ".stl"; ++k) {
                bool wallA = false, wallB = false;
                for (size_t t = 0; t < mesh.kept[k].triangles.size(); ++t) {
                    const std::string n = name(k, t);
                    wallA |= n == "wall A", wallB |= n == "wall B";
                }
                const std::string base = o.part == Part::Solid ? "solid" : wallA && !wallB ? "fluid A"
                                                                        : wallB && !wallA ? "fluid B"
                                                                                          : "fluid";
                const int count = ++seen[base];
                solidNames.push_back(count == 1 ? base : base + " " + std::to_string(count));
            }
            if (ext(out) == ".stl") {
                if (writeStl(out, mesh.kept)) wrote(out);
                else say(print, "error: cannot write '%s'", out.c_str()), ok = false;
                continue;
            }
            // A sheet's fluid: channel A and channel B in a file each.
            std::vector<std::vector<size_t>> groups(1);
            std::vector<std::string> suffixes{""};
            for (size_t k = 0; k < mesh.kept.size(); ++k) groups[0].push_back(k);
            if (splitChannels) {
                std::vector<size_t> channel[2];
                bool both = false;
                for (size_t k = 0; k < solidNames.size(); ++k) {
                    if (solidNames[k].rfind("fluid A", 0) == 0) channel[0].push_back(k);
                    else if (solidNames[k].rfind("fluid B", 0) == 0) channel[1].push_back(k);
                    else both = true;
                }
                if (both || channel[0].empty() || channel[1].empty()) say(print, "note: %s", channelNote(both).c_str());
                else groups = {channel[0], channel[1]}, suffixes = {"_A", "_B"};
            }
            for (size_t g = 0; g < groups.size(); ++g) {
                const std::string file = withSuffix(out, suffixes[g]);
                std::vector<Mesh> parts;
                std::vector<std::string> names;
                for (size_t k : groups[g]) parts.push_back(mesh.kept[k]), names.push_back(solidNames[k]);
                const auto& ks = groups[g];
                if (writeStep(file, parts, p.stem().string() + suffixes[g], names,
                              [&](size_t i, size_t t) { return name(ks[i], t); }))
                    wrote(file);
                else say(print, "error: cannot write '%s'", file.c_str()), ok = false;
            }
            continue;
        }
        const SmoothReport r =
            writeSmoothStep(out, p.stem().string(), lat, cells, o.patches, o.tol, o.minFragment, o.exact,
                            expected * shapeVol + slabVol, shapeVol + slabVol, splitChannels);
        if (!r.error.empty()) {
            if (r.error == "nothing to export") say(print, "error: nothing to export. Check --density / --level.");
            else if (r.error.rfind("cannot write", 0) == 0) say(print, "error: %s", r.error.c_str());
            else say(print, "error: %s: %s at %d patches per cell. Try a higher --patches%s.", out.c_str(),
                     r.error.c_str(), r.patches,
                     s.trimmed() ? ", or --shift 0.1 (walls through the axis or center are hard to fit)" : "");
            ok = false;
            continue;
        }
        say(print, "Smooth    : %d patches per cell, %zu faces%s, %zu solid(s), solid fraction %.4f", r.patches,
            r.faces, r.merged ? (" (merged from " + std::to_string(r.patchFaces) + ")").c_str() : "", r.solids,
            boxFraction(r.volume));
        if (!r.warning.empty()) say(print, "warning: kept patch faces, merge failed: %s.", r.warning.c_str());
        say(print, "Fit error : %.3g mm max (tolerance %.3g mm, %d samples per patch edge)", r.maxError, o.tol,
            r.samples);
        dropped(r.fragments, r.cavities);
        checkFraction(r.volume, "Raise --patches.");
        if (!r.note.empty()) say(print, "note: %s.", r.note.c_str());
        for (const auto& file : r.files) wrote(file);
    }
    if (needReport) {
        const Report rep = makeReport(o, lat, o.res);
        if (o.report) printReport(rep, o, print);
        for (const auto& out : o.outputs) {
            if (!isReport(out)) continue;
            if (ext(out) == ".json" ? writeReportJson(out, rep, o, lat) : writeReportCsv(out, rep, o, lat)) wrote(out);
            else say(print, "error: cannot write '%s'", out.c_str()), ok = false;
        }
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    say(print, "Time      : %.2f s", secs);
    return ok ? 0 : 1;
}

}  // namespace tpms
