#include "export.hpp"
#include "mesher.hpp"
#include "tpms.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

using namespace tpms;

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
  tpms --type gyroid --form sheet --density 0.25 --cell 8 --cells 4 -o gyroid.step -o gyroid.stl
)";

struct Options {
    std::string type = "gyroid";
    Form form = Form::Sheet;
    Part part = Part::Solid;
    double density = 0.2;
    double level = 0;
    bool hasLevel = false;
    double cell[3] = {10, 10, 10};
    double cells[3] = {2, 2, 2};
    int res = 24;
    int patches = 3;
    double tol = 0;  // 0 = 0.001 x smallest cell size
    bool faceted = false;
    double minFragment = 0.01;
    std::vector<std::string> outputs;
};

bool number(const char* s, double& out) {
    char* end = nullptr;
    out = std::strtod(s, &end);
    return end != s && *end == '\0' && std::isfinite(out);
}

// Reads 1 or 3 numbers after the flag at argv[i].
bool numbers(int argc, char** argv, int& i, double out[3]) {
    double v, v1, v2;
    if (i + 1 >= argc || !number(argv[i + 1], v)) return false;
    out[0] = out[1] = out[2] = v;
    ++i;
    if (i + 2 < argc && number(argv[i + 1], v1) && number(argv[i + 2], v2)) {
        out[1] = v1;
        out[2] = v2;
        i += 2;
    }
    return true;
}

int fail(const std::string& msg) {
    std::fprintf(stderr, "error: %s\nRun 'tpms --help' for usage.\n", msg.c_str());
    return 2;
}

std::string lower(std::string s) {
    for (char& ch : s) ch = char(std::tolower(static_cast<unsigned char>(ch)));
    return s;
}

std::string ext(const std::string& path) { return lower(std::filesystem::path(path).extension().string()); }

void wrote(const std::string& path) {
    std::error_code ec;
    const auto bytes = std::filesystem::file_size(path, ec);
    std::printf("Wrote     : %s (%.1f MB)\n", path.c_str(), ec ? 0.0 : double(bytes) / 1e6);
}

}  // namespace

int main(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const bool hasValue = i + 1 < argc;
        double v;
        if (a == "-h" || a == "--help") {
            std::fputs(kUsage, stdout);
            return 0;
        } else if (a == "--list") {
            for (const Surface& s : surfaces()) {
                std::string names = s.name;
                for (const auto& al : s.aliases) names += ", " + al;
                std::printf("  %-28s f = %s\n", names.c_str(), s.formula.c_str());
            }
            return 0;
        } else if (a == "--type" && hasValue) {
            o.type = argv[++i];
        } else if (a == "--form" && hasValue) {
            const std::string f = lower(argv[++i]);
            if (f == "sheet") o.form = Form::Sheet;
            else if (f == "network" || f == "skeletal") o.form = Form::Network;
            else return fail("--form must be sheet or network");
        } else if (a == "--part" && hasValue) {
            const std::string p = lower(argv[++i]);
            if (p == "solid") o.part = Part::Solid;
            else if (p == "void" || p == "fluid") o.part = Part::Void;
            else return fail("--part must be solid or void");
        } else if (a == "--density" && hasValue && number(argv[i + 1], v)) {
            o.density = v;
            ++i;
        } else if (a == "--level" && hasValue && number(argv[i + 1], v)) {
            o.level = v;
            o.hasLevel = true;
            ++i;
        } else if (a == "--cell") {
            if (!numbers(argc, argv, i, o.cell)) return fail("--cell needs 1 or 3 numbers");
        } else if (a == "--cells") {
            if (!numbers(argc, argv, i, o.cells)) return fail("--cells needs 1 or 3 numbers");
        } else if (a == "--res" && hasValue && number(argv[i + 1], v)) {
            o.res = int(v);
            ++i;
        } else if (a == "--patches" && hasValue && number(argv[i + 1], v)) {
            o.patches = int(v);
            ++i;
        } else if (a == "--tol" && hasValue && number(argv[i + 1], v)) {
            o.tol = v;
            ++i;
        } else if (a == "--faceted") {
            o.faceted = true;
        } else if (a == "--min-fragment" && hasValue && number(argv[i + 1], v)) {
            o.minFragment = v;
            ++i;
        } else if ((a == "-o" || a == "--out") && hasValue) {
            o.outputs.push_back(argv[++i]);
        } else {
            return fail("unknown or incomplete option '" + a + "'");
        }
    }

    const Surface* surface = findSurface(o.type);
    if (!surface) return fail("unknown surface type '" + o.type + "'. Run 'tpms --list'.");
    if (!o.hasLevel && !(o.density > 0 && o.density < 1)) return fail("--density must be in (0, 1)");
    int cells[3];
    for (int a = 0; a < 3; ++a) {
        if (!(o.cell[a] > 0)) return fail("--cell must be > 0");
        if (o.cells[a] < 1 || o.cells[a] != std::floor(o.cells[a])) return fail("--cells must be whole numbers >= 1");
        cells[a] = int(o.cells[a]);
    }
    if (o.res < 4) return fail("--res must be >= 4");
    if (o.patches < 2) return fail("--patches must be >= 2");
    if (o.tol < 0) return fail("--tol must be > 0");
    if (o.tol == 0) o.tol = 0.001 * std::min({o.cell[0], o.cell[1], o.cell[2]});
    if (o.outputs.empty()) return fail("no output file. Add -o lattice.step");
    bool needMesh = false;
    for (const auto& out : o.outputs) {
        const std::string e = ext(out);
        if (e != ".step" && e != ".stp" && e != ".stl") return fail("output '" + out + "' must end in .step, .stp or .stl");
        needMesh |= e == ".stl" || o.faceted;
    }

    Lattice lat;
    lat.f = surface->f;
    lat.form = o.form;
    lat.part = o.part;
    for (int a = 0; a < 3; ++a) lat.cell[a] = o.cell[a];
    lat.level = o.hasLevel ? o.level : levelForDensity(surface->f, o.form, o.density);
    if (o.form == Form::Sheet && lat.level <= 0) return fail("sheet needs --level > 0");

    const double boxVol = cells[0] * o.cell[0] * cells[1] * o.cell[1] * cells[2] * o.cell[2];
    auto solidFraction = [&](double vol) { return o.part == Part::Solid ? vol / boxVol : 1.0 - vol / boxVol; };
    // The --density sample covers whole cells, so a large gap means a bad
    // geometry: walls thinner than the grid, or many dropped pieces.
    auto checkFraction = [&](double frac, const char* fix) {
        if (!o.hasLevel && std::fabs(frac - o.density) > 0.05 * o.density)
            std::printf("warning: solid fraction %.4f is far from target %.4f. %s\n", frac, o.density, fix);
    };
    auto dropped = [&](size_t fragments, size_t cavities) {
        if (fragments)
            std::printf("note: dropped %zu small fragment(s) cut off by the box (< %g x largest).\n", fragments,
                        o.minFragment);
        if (cavities) std::printf("warning: dropped %zu inner cavity shell(s). Not supported yet.\n", cavities);
    };

    std::printf("Surface   : %s, %s, %s part\n", surface->name.c_str(),
                o.form == Form::Sheet ? "sheet" : "network", o.part == Part::Solid ? "solid" : "void");
    if (o.hasLevel) std::printf("Level c   : %.5g\n", lat.level);
    else std::printf("Level c   : %.5g (from target solid fraction %.3f)\n", lat.level, o.density);
    std::printf("Box       : %g x %g x %g mm (%d x %d x %d cells of %g x %g x %g mm)\n",
                cells[0] * o.cell[0], cells[1] * o.cell[1], cells[2] * o.cell[2], cells[0], cells[1], cells[2],
                o.cell[0], o.cell[1], o.cell[2]);

    const auto t0 = std::chrono::steady_clock::now();
    bool ok = true;
    std::vector<Mesh> kept;
    if (needMesh) {
        size_t degenerate = 0;
        const Mesh mesh = meshLattice(lat, cells, o.res, &degenerate);
        const ClosedCheck check = checkClosed(mesh);
        std::vector<Mesh> pieces = splitComponents(mesh);
        std::vector<double> vols;
        double maxVol = 0;
        for (const Mesh& p : pieces) {
            vols.push_back(volume(p));
            maxVol = std::max(maxVol, vols.back());
        }
        size_t fragments = 0, cavities = 0, triangles = 0;
        double keptVol = 0;
        for (size_t i = 0; i < pieces.size(); ++i) {
            if (vols[i] <= 0) ++cavities;
            else if (vols[i] < o.minFragment * maxVol) ++fragments;
            else {
                keptVol += vols[i];
                triangles += pieces[i].triangles.size();
                kept.push_back(std::move(pieces[i]));
            }
        }
        if (kept.empty()) return fail("nothing to export. Check --density / --level.");
        std::printf("Mesh      : %d steps per cell, %zu triangles, %zu solid(s), closed: %s, solid fraction %.4f\n",
                    o.res, triangles, kept.size(), check.closed() ? "yes" : "NO", solidFraction(keptVol));
        if (!check.closed()) {
            std::printf("warning: mesh not closed (%zu open, %zu non-manifold edges). Try another --res.\n",
                        check.openEdges, check.nonManifoldEdges);
            ok = false;
        }
        if (degenerate) std::printf("warning: dropped %zu zero-area triangles.\n", degenerate);
        dropped(fragments, cavities);
        checkFraction(solidFraction(keptVol), "Raise --res.");
        if (triangles > 300000)
            std::printf("note: %zu triangles is heavy for CAD. Lower --res.\n", triangles);
    }

    for (const auto& out : o.outputs) {
        const std::filesystem::path p(out);
        if (ext(out) == ".stl" || o.faceted) {
            const bool written = ext(out) == ".stl" ? writeStl(out, kept) : writeStep(out, kept, p.stem().string());
            if (!written) {
                std::fprintf(stderr, "error: cannot write '%s'\n", out.c_str());
                ok = false;
                continue;
            }
            wrote(out);
            continue;
        }
        const SmoothReport r = writeSmoothStep(out, p.stem().string(), lat, cells, o.patches, o.tol, o.minFragment);
        if (!r.error.empty()) {
            std::fprintf(stderr, "error: %s: %s at %d patches per cell. Try a higher --patches.\n", out.c_str(),
                         r.error.c_str(), r.patches);
            ok = false;
            continue;
        }
        std::printf("Smooth    : %d patches per cell, %zu faces, %zu solid(s), solid fraction %.4f\n", r.patches,
                    r.faces, r.solids, solidFraction(r.volume));
        std::printf("Fit error : %.3g mm max (tolerance %.3g mm, %d samples per patch edge)\n", r.maxError, o.tol,
                    r.samples);
        dropped(r.fragments, r.cavities);
        checkFraction(solidFraction(r.volume), "Raise --patches.");
        wrote(out);
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("Time      : %.2f s\n", secs);
    return ok ? 0 : 1;
}
