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

STEP (smooth B-spline faces; box, tube and torus-shell)
  --patches N           Patch grid steps per unit cell edge to start from
                        (default 3). The tool steps it up until the geometry
                        is valid and within --tol.
  --tol T               Max fit error in mm (default 0.001 x smallest cell
                        edge).
  --exact               Keep one face per patch, with exactly shared edges.
                        Default: merge patches into large trimmed faces whose
                        edges lie within --tol of the faces.
  --faceted             Write one planar face per triangle instead (v0).
                        Works for every shape.

STL and faceted STEP
  --res N               Grid steps per unit cell edge (default 24).
                        Thin walls need more. The tool warns when too low.

Both
  --min-fragment F      Drop pieces smaller than F x largest piece volume
                        (default 0.01). Use 0 to keep all pieces.

Output
  -o, --out FILE        .step / .stp, .stl (binary), or .json (the property
                        report). Repeat for several files.
  --report              Print the property report: solid and fluid volumes,
                        wetted area, hydraulic diameters, wall thickness,
                        and which faces of the shape each piece touches.
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

// Reads up to `max` numbers after the flag at args[i]. Returns how many.
int numbers(const std::vector<std::string>& args, size_t& i, double* out, int max) {
    int n = 0;
    while (n < max && i + 1 < args.size() && number(args[i + 1], out[n])) ++i, ++n;
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
    bool cellGiven = false, cellsGiven = false;
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
            o.hasLevel |= isLevel;
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
            o.res = int(v);
            ++i;
        } else if (a == "--patches" && hasValue && number(args[i + 1], v)) {
            o.patches = int(v);
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
        if (cells[a] < 1 || cells[a] != std::floor(cells[a])) return error = "--cells must be whole numbers >= 1", false;
        o.cells[a] = int(cells[a]);
        s.cells[a] = cells[a];
        s.cellSize[a] = o.cell[a];
    }
    if (o.res < 4) return error = "--res must be >= 4", false;
    if (o.patches < 2) return error = "--patches must be >= 2", false;
    if (o.tol < 0) return error = "--tol must be > 0", false;
    if (o.tol == 0) o.tol = autoTol(o);
    for (const auto& out : o.outputs) {
        const std::string e = ext(out);
        if (e != ".step" && e != ".stp" && e != ".stl" && e != ".json")
            return error = "output '" + out + "' must end in .step, .stp, .stl or .json", false;
        if ((e == ".step" || e == ".stp") && !o.faceted && s.trimmed())
            return error = std::string("smooth STEP does not support --shape ") + shapeName(s.kind) +
                           " yet: write .stl or add --faceted", false;
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
    g.periodic = s.periodic();
    return g;
}

Grid meshGrid(const Options& o, int res) {
    Grid g = gridOf(o, res);
    const int a = o.shape.slabAxis;
    if (a < 0) return g;
    const double h = g.size[a] / g.steps[a];
    const int n0 = int(std::ceil(o.shape.slab[0] / h - 1e-9)), n1 = int(std::ceil(o.shape.slab[1] / h - 1e-9));
    g.origin[a] -= n0 * h, g.steps[a] += n0 + n1, g.size[a] = g.steps[a] * h;
    return g;
}

double slabVolume(const Options& o) {
    const int a = o.shape.slabAxis;
    if (a < 0) return 0;
    return (o.shape.slab[0] + o.shape.slab[1]) * o.shape.volume() / o.shape.boxPos(a, o.cells[a]);
}

Mesh meshSolid(const Options& o, const Lattice& lat, int res, size_t* droppedDegenerate) {
    // Plates or buffers: the lattice in its box, joined with a slab on each
    // side of the axis that has one. As fields: inside the outer ends, either
    // the lattice or the slabs (cell coordinate p, box from 0 to hi).
    // Written so that the lattice's own region stays inside on the faces
    // where it meets a slab, which a min of two box fields would not do.
    const Shape& s = o.shape;
    const int a = s.slabAxis;
    const double hi = a < 0 ? 0 : gridOf(o, 1).size[a], none = 1e300;
    const Field field = [&](double x, double y, double z) {
        const double g = lat(x, y, z);
        if (a < 0) return g;
        const double p = a == 0 ? x : a == 1 ? y : z;
        const double ends = std::max(s.slab[0] > 0 ? -s.slab[0] - p : -none, s.slab[1] > 0 ? p - hi - s.slab[1] : -none);
        const double inBox = std::min(s.slab[0] > 0 ? p : none, s.slab[1] > 0 ? hi - p : none);
        return std::max(ends, std::min(g, inBox));
    };
    Mesh m = meshField(field, meshGrid(o, res), droppedDegenerate);
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
    for (size_t i = 0; i < pieces.size(); ++i) {
        if (vols[i] <= 0) ++out.cavities;
        else if (vols[i] < minFragment * maxVol) ++out.fragments;
        else {
            out.volume += vols[i];
            out.triangles += pieces[i].triangles.size();
            out.kept.push_back(std::move(pieces[i]));
        }
    }
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
    for (const auto& out : o.outputs) {
        needMesh |= ext(out) == ".stl" || (o.faceted && ext(out) != ".json");
        needReport |= ext(out) == ".json";
    }

    const Lattice lat = makeLattice(o);
    const int* cells = o.cells;
    const Shape& s = o.shape;
    const double shapeVol = shapeVolume(o), slabVol = slabVolume(o);
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
            say(print, "note: dropped %zu small fragment(s) cut off by the shape (< %g x largest).", fragments,
                o.minFragment);
        if (cavities) say(print, "warning: dropped %zu inner cavity shell(s). Not supported yet.", cavities);
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

    for (const auto& out : o.outputs) {
        const std::filesystem::path p(out);
        if (ext(out) == ".json") continue;  // the report, below
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
            const bool written = ext(out) == ".stl"
                                     ? writeStl(out, mesh.kept)
                                     : writeStep(out, mesh.kept, p.stem().string(), solidNames, name);
            if (!written) {
                say(print, "error: cannot write '%s'", out.c_str());
                ok = false;
                continue;
            }
            wrote(out);
            continue;
        }
        const SmoothReport r =
            writeSmoothStep(out, p.stem().string(), lat, cells, o.patches, o.tol, o.minFragment, o.exact,
                            expected * shapeVol + slabVol, shapeVol + slabVol);
        if (!r.error.empty()) {
            say(print, "error: %s: %s at %d patches per cell. Try a higher --patches.", out.c_str(), r.error.c_str(),
                r.patches);
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
        wrote(out);
    }
    if (needReport) {
        const Report rep = makeReport(o, lat, o.res);
        if (o.report) printReport(rep, o, print);
        for (const auto& out : o.outputs) {
            if (ext(out) != ".json") continue;
            if (writeReportJson(out, rep, o, lat)) wrote(out);
            else say(print, "error: cannot write '%s'", out.c_str()), ok = false;
        }
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    say(print, "Time      : %.2f s", secs);
    return ok ? 0 : 1;
}

}  // namespace tpms
