#include "run.hpp"

#include "export.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <filesystem>

namespace tpms {
namespace {

const char* kUsage = R"(tpms - TPMS lattice to STEP / STL

Usage: tpms [options] -o FILE [-o FILE ...]

Lattice
  --type NAME           Surface type (default gyroid). See --list.
  --form sheet|network  sheet: |f| <= c, a wall around the surface.
                        network: f <= c, one of the two channel networks.
                        Default sheet.
  --density R           Target solid volume fraction of one unit cell,
                        0 < R < 1 (default 0.2).
  --level C             Iso level c. Use instead of --density.
  --part solid|void     solid, or void = shape minus solid (fluid domain for
                        CFD). --density always means the solid fraction.

Shape (default box)
  --shape NAME          box | tube | torus-shell   (cells follow the shape)
                        cylinder | sphere | torus  (straight cells, trimmed)
  --cell L [Ly Lz]      Unit cell size in mm (box and trimmed; default 10).
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

STEP (smooth B-spline faces; box only so far)
  --patches N           Patch grid steps per unit cell edge to start from
                        (default 3). The tool steps it up until the geometry
                        is valid and within --tol.
  --tol T               Max fit error in mm (default 0.001 x cell size).
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
  -o, --out FILE        .step / .stp or .stl (binary). Repeat for several files.
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

const char* shapeName(ShapeKind k) {
    for (const auto& s : kShapes)
        if (s.kind == k) return s.name;
    return "?";
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

}  // namespace

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
        } else if (a == "--density" && hasValue && number(args[i + 1], v)) {
            o.density = v;
            ++i;
        } else if (a == "--level" && hasValue && number(args[i + 1], v)) {
            o.level = v;
            o.hasLevel = true;
            ++i;
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
    if (!o.hasLevel && !(o.density > 0 && o.density < 1)) return error = "--density must be in (0, 1)", false;
    if (o.hasLevel && o.form == Form::Sheet && o.level <= 0) return error = "sheet needs --level > 0", false;

    Shape& s = o.shape;
    const std::string name = shapeName(s.kind);
    if (s.conformal() && cellGiven) return error = "--cell does not apply to --shape " + name + ": set --cells", false;
    if (s.trimmed() && cellsGiven) return error = "--cells does not apply to --shape " + name + ": set --cell", false;
    if (s.kind == ShapeKind::Box && (radii || s.height || s.major))
        return error = "--radius, --height and --major need a --shape", false;
    const bool shell = s.conformal(), needsHeight = s.kind == ShapeKind::Tube || s.kind == ShapeKind::Cylinder;
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
    }
    if (o.res < 4) return error = "--res must be >= 4", false;
    if (o.patches < 2) return error = "--patches must be >= 2", false;
    if (o.tol < 0) return error = "--tol must be > 0", false;
    if (o.tol == 0) o.tol = 0.001 * std::min({o.cell[0], o.cell[1], o.cell[2]});
    for (const auto& out : o.outputs) {
        const std::string e = ext(out);
        if (e != ".step" && e != ".stp" && e != ".stl")
            return error = "output '" + out + "' must end in .step, .stp or .stl", false;
        if (e != ".stl" && !o.faceted && s.kind != ShapeKind::Box)
            return error = "smooth STEP supports --shape box only so far: write .stl or add --faceted", false;
    }
    return true;
}

Lattice makeLattice(const Options& o) {
    const Surface* surface = findSurface(o.type);
    Lattice lat;
    lat.f = surface->f;
    lat.form = o.form;
    lat.part = o.part;
    lat.shape = o.shape;
    for (int a = 0; a < 3; ++a) lat.cell[a] = o.shape.conformal() ? 1 : o.cell[a];
    lat.level = o.hasLevel ? o.level : levelForDensity(surface->f, o.form, o.density);
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
            g.size[a] = s.conformal() ? o.cells[a] : o.cells[a] * o.cell[a];
        }
    }
    g.periodic = s.periodic();
    return g;
}

Mesh meshSolid(const Options& o, const Lattice& lat, int res, size_t* droppedDegenerate) {
    Mesh m = meshField([&](double x, double y, double z) { return lat(x, y, z); }, gridOf(o, res), droppedDegenerate);
    if (o.shape.conformal())
        for (Vec3& v : m.vertices) v = o.shape.map(v);
    return m;
}

double shapeVolume(const Options& o) {
    if (o.shape.kind != ShapeKind::Box) return o.shape.volume();
    return o.cells[0] * o.cell[0] * o.cells[1] * o.cell[1] * o.cells[2] * o.cell[2];
}

double sampledFraction(const Options& o, const Lattice& lat) {
    // Quasi-random R3 points: a regular grid can line up with the lattice.
    // Conformal points are weighted by the volume the mapping gives them.
    const Grid g = gridOf(o, 1);
    const double phi = 1.2207440846057595, step[3] = {1 / phi, 1 / (phi * phi), 1 / (phi * phi * phi)};
    double in = 0, all = 0;
    for (int k = 0; k < (1 << 17); ++k) {
        Vec3 p;
        for (int a = 0; a < 3; ++a) p[a] = g.origin[a] + std::fmod(0.5 + k * step[a], 1.0) * g.size[a];
        double w = 1;
        if (o.shape.conformal()) {
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
    if (o.outputs.empty()) return print("error: no output file. Add -o lattice.step"), 2;
    bool needMesh = false;
    for (const auto& out : o.outputs) needMesh |= ext(out) == ".stl" || o.faceted;

    const Lattice lat = makeLattice(o);
    const int* cells = o.cells;
    const Shape& s = o.shape;
    const double shapeVol = shapeVolume(o);
    auto solidFraction = [&](double exported) { return o.part == Part::Solid ? exported : 1.0 - exported; };
    // What the field gives, from samples. Bending a conformal lattice moves
    // it away from --density, so compare meshes with this, not the target.
    const double expected = sampledFraction(o, lat);
    auto checkFraction = [&](double vol, const char* fix) {
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
    if (o.hasLevel) say(print, "Level c   : %.5g", lat.level);
    else say(print, "Level c   : %.5g (from target solid fraction %.3f)", lat.level, o.density);
    switch (s.kind) {
    case ShapeKind::Box:
        say(print, "Box       : %g x %g x %g mm (%d x %d x %d cells of %g x %g x %g mm)", cells[0] * o.cell[0],
            cells[1] * o.cell[1], cells[2] * o.cell[2], cells[0], cells[1], cells[2], o.cell[0], o.cell[1], o.cell[2]);
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
    if (std::fabs(solidFraction(expected) - o.density) > 0.01 && !o.hasLevel)
        say(print, "note: the field gives solid fraction %.4f here, not %.3f: the shape bends or cuts the cells.",
            solidFraction(expected), o.density);

    const auto t0 = std::chrono::steady_clock::now();
    bool ok = true;
    Pieces mesh;
    if (needMesh) {
        size_t degenerate = 0;
        const Mesh full = meshSolid(o, lat, o.res, &degenerate);
        const ClosedCheck check = checkClosed(full);
        mesh = keepPieces(full, o.minFragment);
        if (mesh.kept.empty()) return print("error: nothing to export. Check --density / --level."), 2;
        say(print, "Mesh      : %d steps per cell, %zu triangles, %zu solid(s), closed: %s, solid fraction %.4f", o.res,
            mesh.triangles, mesh.kept.size(), check.closed() ? "yes" : "NO", solidFraction(mesh.volume / shapeVol));
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
        if (ext(out) == ".stl" || o.faceted) {
            const bool written =
                ext(out) == ".stl" ? writeStl(out, mesh.kept) : writeStep(out, mesh.kept, p.stem().string());
            if (!written) {
                say(print, "error: cannot write '%s'", out.c_str());
                ok = false;
                continue;
            }
            wrote(out);
            continue;
        }
        const SmoothReport r =
            writeSmoothStep(out, p.stem().string(), lat, cells, o.patches, o.tol, o.minFragment, o.exact);
        if (!r.error.empty()) {
            say(print, "error: %s: %s at %d patches per cell. Try a higher --patches.", out.c_str(), r.error.c_str(),
                r.patches);
            ok = false;
            continue;
        }
        say(print, "Smooth    : %d patches per cell, %zu faces%s, %zu solid(s), solid fraction %.4f", r.patches,
            r.faces, r.merged ? (" (merged from " + std::to_string(r.patchFaces) + ")").c_str() : "", r.solids,
            solidFraction(r.volume / shapeVol));
        if (!r.warning.empty()) say(print, "warning: kept patch faces, merge failed: %s.", r.warning.c_str());
        say(print, "Fit error : %.3g mm max (tolerance %.3g mm, %d samples per patch edge)", r.maxError, o.tol,
            r.samples);
        dropped(r.fragments, r.cavities);
        checkFraction(r.volume, "Raise --patches.");
        wrote(out);
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    say(print, "Time      : %.2f s", secs);
    return ok ? 0 : 1;
}

}  // namespace tpms
