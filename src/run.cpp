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

Geometry
  --type NAME           Surface type (default gyroid). See --list.
  --form sheet|network  sheet: |f| <= c, a wall around the surface.
                        network: f <= c, one of the two channel networks.
                        Default sheet.
  --density R           Target solid volume fraction, 0 < R < 1 (default 0.2).
  --level C             Iso level c. Use instead of --density.
  --cell L [Ly Lz]      Unit cell size in mm (default 10).
  --cells N [Ny Nz]     Number of cells per axis (default 2).
  --part solid|void     solid, or void = box minus solid (fluid domain for CFD).
                        --density always means the solid fraction.

STEP (smooth B-spline faces)
  --patches N           Patch grid steps per unit cell edge to start from
                        (default 3). The tool steps it up until the geometry
                        is valid and within --tol.
  --tol T               Max fit error in mm (default 0.001 x cell size).
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
  -o, --out FILE        .step / .stp or .stl (binary). Repeat for several files.
  --list                List surface types.
  -h, --help            Show this help.

Example
  tpms --type gyroid --form sheet --density 0.25 --cell 8 --cells 4 -o gyroid.step -o gyroid.stl)";

bool number(const std::string& s, double& out) {
    char* end = nullptr;
    out = std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0' && std::isfinite(out);
}

// Reads 1 or 3 numbers after the flag at args[i].
bool numbers(const std::vector<std::string>& args, size_t& i, double out[3]) {
    double v, v1, v2;
    if (i + 1 >= args.size() || !number(args[i + 1], v)) return false;
    out[0] = out[1] = out[2] = v;
    ++i;
    if (i + 2 < args.size() && number(args[i + 1], v1) && number(args[i + 2], v2)) {
        out[1] = v1;
        out[2] = v2;
        i += 2;
    }
    return true;
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
    for (size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        const bool hasValue = i + 1 < args.size();
        double v;
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
        } else if (a == "--cell") {
            if (!numbers(args, i, o.cell)) return error = "--cell needs 1 or 3 numbers", false;
        } else if (a == "--cells") {
            if (!numbers(args, i, cells)) return error = "--cells needs 1 or 3 numbers", false;
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
    for (int a = 0; a < 3; ++a) {
        if (!(o.cell[a] > 0)) return error = "--cell must be > 0", false;
        if (cells[a] < 1 || cells[a] != std::floor(cells[a])) return error = "--cells must be whole numbers >= 1", false;
        o.cells[a] = int(cells[a]);
    }
    if (o.res < 4) return error = "--res must be >= 4", false;
    if (o.patches < 2) return error = "--patches must be >= 2", false;
    if (o.tol < 0) return error = "--tol must be > 0", false;
    if (o.tol == 0) o.tol = 0.001 * std::min({o.cell[0], o.cell[1], o.cell[2]});
    for (const auto& out : o.outputs) {
        const std::string e = ext(out);
        if (e != ".step" && e != ".stp" && e != ".stl")
            return error = "output '" + out + "' must end in .step, .stp or .stl", false;
    }
    return true;
}

Lattice makeLattice(const Options& o) {
    const Surface* surface = findSurface(o.type);
    Lattice lat;
    lat.f = surface->f;
    lat.form = o.form;
    lat.part = o.part;
    for (int a = 0; a < 3; ++a) lat.cell[a] = o.cell[a];
    lat.level = o.hasLevel ? o.level : levelForDensity(surface->f, o.form, o.density);
    return lat;
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
    const double boxVol = cells[0] * o.cell[0] * cells[1] * o.cell[1] * cells[2] * o.cell[2];
    auto solidFraction = [&](double vol) { return o.part == Part::Solid ? vol / boxVol : 1.0 - vol / boxVol; };
    // The --density sample covers whole cells, so a large gap means a bad
    // geometry: walls thinner than the grid, or many dropped pieces.
    auto checkFraction = [&](double frac, const char* fix) {
        if (!o.hasLevel && std::fabs(frac - o.density) > 0.05 * o.density)
            say(print, "warning: solid fraction %.4f is far from target %.4f. %s", frac, o.density, fix);
    };
    auto dropped = [&](size_t fragments, size_t cavities) {
        if (fragments)
            say(print, "note: dropped %zu small fragment(s) cut off by the box (< %g x largest).", fragments,
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
    say(print, "Box       : %g x %g x %g mm (%d x %d x %d cells of %g x %g x %g mm)", cells[0] * o.cell[0],
        cells[1] * o.cell[1], cells[2] * o.cell[2], cells[0], cells[1], cells[2], o.cell[0], o.cell[1], o.cell[2]);

    const auto t0 = std::chrono::steady_clock::now();
    bool ok = true;
    Pieces mesh;
    if (needMesh) {
        size_t degenerate = 0;
        const Mesh full = meshLattice(lat, cells, o.res, &degenerate);
        const ClosedCheck check = checkClosed(full);
        mesh = keepPieces(full, o.minFragment);
        if (mesh.kept.empty()) return print("error: nothing to export. Check --density / --level."), 2;
        say(print, "Mesh      : %d steps per cell, %zu triangles, %zu solid(s), closed: %s, solid fraction %.4f", o.res,
            mesh.triangles, mesh.kept.size(), check.closed() ? "yes" : "NO", solidFraction(mesh.volume));
        if (!check.closed()) {
            say(print, "warning: mesh not closed (%zu open, %zu non-manifold edges). Try another --res.",
                check.openEdges, check.nonManifoldEdges);
            ok = false;
        }
        if (degenerate) say(print, "warning: dropped %zu zero-area triangles.", degenerate);
        dropped(mesh.fragments, mesh.cavities);
        checkFraction(solidFraction(mesh.volume), "Raise --res.");
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
            solidFraction(r.volume));
        if (!r.warning.empty()) say(print, "warning: kept patch faces, merge failed: %s.", r.warning.c_str());
        say(print, "Fit error : %.3g mm max (tolerance %.3g mm, %d samples per patch edge)", r.maxError, o.tol,
            r.samples);
        dropped(r.fragments, r.cavities);
        checkFraction(solidFraction(r.volume), "Raise --patches.");
        wrote(out);
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    say(print, "Time      : %.2f s", secs);
    return ok ? 0 : 1;
}

}  // namespace tpms
