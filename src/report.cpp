// Property report: volumes, areas and connectivity of the solid and the
// fluid (void), from closed meshes of both.

#include "run.hpp"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>

namespace tpms {
namespace {

void say(const Print& print, const char* fmt, ...) {
    char buf[1024];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    print(buf);
}

const char* channelName(int c) { return c ? "f < -c" : "f > c"; }

std::string joined(const std::vector<std::string>& v, const char* sep) {
    std::string s;
    for (const auto& x : v) s += (s.empty() ? "" : sep) + x;
    return s;
}

// A number for JSON (null when not finite, e.g. no wetted area).
std::string num(double v) {
    if (!std::isfinite(v)) return "null";
    char b[32];
    std::snprintf(b, sizeof b, "%.6g", v);
    return b;
}

std::string quoted(const std::string& s) { return "\"" + s + "\""; }

// Estimates per fluid channel (or the fluid of a network), from its
// measured volume and wetted area, by packed-bed correlations: rough
// values, for comparing designs. Permeability by Kozeny-Carman with
// Kozeny's constant 5: K = porosity Dh^2 / 80. With --flow (superficial
// velocity u): the
// Reynolds number on Dh and the speed in the pores, and the pressure
// gradient mu u / K (Darcy, viscous) + 7/6 rho u^2 / (porosity^2 Dh)
// (Ergun's inertial term 1.75 rho (1 - e) u^2 / (e^3 dp), with the packed
// bed's dp = 1.5 Dh (1 - e) / e). The drop is over the box along the
// channel's flow axis: its channel buffers' or buffers' axis, else x.
struct ChannelFlow {
    std::string name;                 // "f > c", "f < -c", or "fluid"
    double volume = 0, area = 0;      // mm^3, mm^2
    double porosity = 0, dh = 0;      // of the shape volume; hydraulic diameter, mm
    double permeability = 0;          // m^2
    int axis = -1;                    // box: flow axis
    double length = 0;                // box: along it, mm
    double re = 0, viscous = 0, gradient = 0, drop = 0;  // with --flow: Reynolds number; Pa/m (viscous part, total); Pa
};

std::vector<ChannelFlow> channelFlows(const Report& r, const Options& o) {
    std::vector<ChannelFlow> out;
    const Shape& s = o.shape;
    for (int c = 0; c < r.channels; ++c) {
        ChannelFlow f;
        f.name = r.channels == 2 ? channelName(c) : "fluid";
        f.volume = r.channels == 2 ? r.channelVolume[c] : r.fluidVolume;
        f.area = r.channels == 2 ? r.channelArea[c] : r.wettedArea;
        f.porosity = f.volume / r.shapeVolume, f.dh = 4 * f.volume / f.area;
        f.permeability = f.porosity * f.dh * f.dh / 80 * 1e-6;
        if (s.kind == ShapeKind::Box) {
            f.axis = s.channelAxis[c] >= 0 ? s.channelAxis[c] : std::max(s.slabAxis, 0);
            f.length = gridOf(o, 1).size[f.axis];
        }
        if (o.flow[0] > 0) {
            const double u = o.flow[0], mu = o.flow[1], rho = o.flow[2];
            f.re = rho * (u / f.porosity) * f.dh * 1e-3 / mu;
            f.viscous = mu * u / f.permeability;
            f.gradient = f.viscous + 7.0 / 6 * rho * u * u / (f.porosity * f.porosity * f.dh * 1e-3);
            f.drop = f.gradient * f.length * 1e-3;
        }
        out.push_back(f);
    }
    return out;
}

// Gibson-Ashby bounds of the relative Young's modulus E / Es at solid
// fraction phi: bending-dominated phi^2, stretching-dominated phi / 3.
std::pair<double, double> stiffness(const Report& r) {
    const double phi = r.solidVolume / r.shapeVolume;
    return {phi * phi, phi / 3};
}

const char* const kAxisNames[3] = {"x", "y", "z"};

// Pairs of opposite faces that each piece connects, as "solid" and the
// fluid channels.
struct Through {
    bool solid = false;
    std::vector<std::string> fluid;
};
std::vector<Through> through(const Report& r) {
    std::vector<Through> out;
    auto touches = [](const ReportPiece& p, const std::string& f) {
        return std::find(p.faces.begin(), p.faces.end(), f) != p.faces.end();
    };
    for (const auto& [a, b] : r.opposite) {
        Through t;
        for (const auto& p : r.solids) t.solid |= touches(p, a) && touches(p, b);
        for (int c = 0; c < r.channels; ++c)
            for (const auto& p : r.fluids)
                if (p.channel == c && touches(p, a) && touches(p, b)) {
                    t.fluid.push_back(r.channels == 2 ? channelName(c) : "fluid");
                    break;
                }
        out.push_back(t);
    }
    return out;
}

}  // namespace

Report makeReport(const Options& o, const Lattice& lat, int res) {
    Report r;
    r.res = res;
    r.shapeVolume = shapeVolume(o);
    const Shape& s = o.shape;
    const Grid grid = gridOf(o, res);
    double h = 1e300;
    for (int a = 0; a < 3; ++a) h = std::min(h, grid.size[a] / grid.steps[a]);
    r.step = std::cbrt(r.shapeVolume / (double(grid.steps[0]) * grid.steps[1] * grid.steps[2]));
    if (s.trimmed()) r.step = h;  // the grid covers more than the shape
    const bool sheet = lat.form == Form::Sheet;
    r.channels = sheet ? 2 : 1;
    Lattice open = lat;  // the solid's lattice field without the trimming shape
    open.shape = Shape{};
    open.part = Part::Solid;
    auto toMm = [&](Vec3 p) { return s.mapped() ? s.map(p) : p; };

    // Faces of the shape, and the ones a grid point lies on (bitmask).
    std::vector<std::string> names;
    std::function<int(Vec3)> facesAt;
    if (!s.trimmed()) {
        for (int f = 0; f < 6; ++f) names.push_back(s.faceName(f));
        for (int a = 0; a < 3; ++a)
            if (!names[size_t(2 * a)].empty()) r.opposite.push_back({names[size_t(2 * a)], names[size_t(2 * a + 1)]});
        facesAt = [&](Vec3 p) {
            int m = 0;
            for (int a = 0; a < 3; ++a)
                for (int side = 0; side < 2 && !(grid.periodic >> a & 1); ++side)
                    if (std::fabs(p[a] - grid.origin[a] - side * grid.size[a]) < 1e-6 * h) m |= 1 << (2 * a + side);
            return m;
        };
    } else if (s.kind == ShapeKind::Cylinder) {
        names = {"side", "bottom", "top"};
        r.opposite.push_back({"bottom", "top"});
        facesAt = [&](Vec3 p) {
            const double tol = 0.05 * h;
            return int(std::fabs(std::hypot(p.x, p.y) - s.r1) < tol) | int(std::fabs(p.z + 0.5 * s.height) < tol) << 1 |
                   int(std::fabs(p.z - 0.5 * s.height) < tol) << 2;
        };
    } else {
        names = {"surface"};
        facesAt = [&](Vec3 p) { return int(std::fabs(s.inside(p.x, p.y, p.z)) < 0.05 * h); };
    }

    // Pieces of a closed mesh: volume in mm, faces touched. Volumes include
    // inner cavities (negative pieces), which are counted apart.
    auto pieces = [&](const Mesh& mesh, int channel, std::vector<ReportPiece>& out) {
        double total = 0;
        for (Mesh m : splitComponents(mesh)) {
            int mask = 0;
            for (const Vec3& v : m.vertices) mask |= facesAt(v);
            for (Vec3& v : m.vertices) v = toMm(v);
            ReportPiece p;
            p.volume = volume(m), p.channel = channel;
            total += p.volume;
            if (p.volume <= 0) {
                ++r.cavities;
                continue;
            }
            for (size_t f = 0; f < names.size(); ++f)
                if ((mask >> f & 1) && !names[f].empty()) p.faces.push_back(names[f]);
            out.push_back(p);
        }
        return total;
    };
    // Area of a closed mesh in mm, without the cut faces on the shape;
    // `wall` marks the vertices of the rest.
    auto wallArea = [&](const Mesh& mesh, const Field& field, double& total, std::vector<char>& wall) {
        std::vector<Vec3> mm(mesh.vertices.size());
        for (size_t i = 0; i < mm.size(); ++i) mm[i] = toMm(mesh.vertices[i]);
        wall.assign(mm.size(), 0);
        double area = 0;
        for (const auto& t : mesh.triangles) {
            const Vec3 &a = mesh.vertices[t[0]], &b = mesh.vertices[t[1]], &c = mesh.vertices[t[2]];
            const double ta = 0.5 * norm(cross(mm[t[1]] - mm[t[0]], mm[t[2]] - mm[t[0]]));
            total += ta;
            const Vec3 m = (a + b + c) * (1.0 / 3);
            const bool cut = s.trimmed() ? std::fabs(s.inside(m.x, m.y, m.z)) < std::fabs(field(m.x, m.y, m.z))
                                         : (facesAt(a) & facesAt(b) & facesAt(c)) != 0;
            if (cut) continue;
            area += ta;
            for (uint32_t v : t) wall[v] = 1;
        }
        return area;
    };

    // The solid as exported.
    Lattice solidLat = lat;
    solidLat.part = Part::Solid;
    const Mesh solid = meshField([&](double x, double y, double z) { return solidLat(x, y, z); }, grid);
    r.solidVolume = pieces(solid, 0, r.solids);
    // Its cut faces are flat (or on the trimming shape) and meshed well;
    // its walls come from the channel meshes below, which follow thin
    // sheets better.
    double solidArea = 0;
    std::vector<char> unused;
    const double solidWalls = wallArea(solid, open, solidArea, unused);
    r.cutArea = solidArea - solidWalls;

    // The fluid, one mesh per channel: f > c, and for a sheet f < -c. Each
    // is a single level set, so a wall thinner than the grid cannot join
    // the two as it can in one mesh of the whole fluid. Their walls are the
    // wetted area; across a sheet, f changes by 2 c, so the local wall
    // thickness is 2 c / |grad f|.
    r.thickness[0] = 1e300, r.thickness[1] = 0;
    for (int c = 0; c < r.channels; ++c) {
        const double sgn = c ? 1 : -1;  // channel 0: c - f < 0, channel 1: f + c < 0
        const Field untrimmed = [&, sgn](double x, double y, double z) {
            return sgn * lat.value(x, y, z) + lat.levelAt(x, y, z);
        };
        const Field channel = [&](double x, double y, double z) {
            return s.trimmed() ? std::max(untrimmed(x, y, z), s.inside(x, y, z)) : untrimmed(x, y, z);
        };
        const Mesh mesh = meshField(channel, grid);
        r.channelVolume[c] = pieces(mesh, c, r.fluids);
        r.fluidVolume += r.channelVolume[c];
        double all = 0;
        std::vector<char> wall;
        r.channelArea[c] = wallArea(mesh, untrimmed, all, wall);
        r.wettedArea += r.channelArea[c];
        const double d = 1e-4 * h;
        for (size_t i = 0; i < wall.size() && sheet; ++i) {
            const Vec3 p = mesh.vertices[i];
            if (!wall[i] || std::fabs(untrimmed(p.x, p.y, p.z)) > 0.01 * lat.levelAt(p.x, p.y, p.z)) continue;
            Vec3 g;
            for (int a = 0; a < 3; ++a) {
                Vec3 q = p, q2 = p;
                q[a] += d, q2[a] -= d;
                g[a] = (lat.value(q.x, q.y, q.z) - lat.value(q2.x, q2.y, q2.z)) / (2 * d);
            }
            const double t = 2 * lat.levelAt(p.x, p.y, p.z) / norm(s.mapped() ? s.gradient(p, g) : g);
            r.thickness[0] = std::min(r.thickness[0], t), r.thickness[1] = std::max(r.thickness[1], t);
        }
    }
    if (r.thickness[0] > r.thickness[1]) r.thickness[0] = r.thickness[1] = 0;
    auto bigFirst = [](const ReportPiece& a, const ReportPiece& b) { return a.volume > b.volume; };
    std::sort(r.solids.begin(), r.solids.end(), bigFirst);
    std::sort(r.fluids.begin(), r.fluids.end(), bigFirst);
    return r;
}

std::vector<std::pair<std::string, std::string>> reportRows(const Report& r, const Options& o) {
    const double V = r.shapeVolume, Vs = r.solidVolume, Vf = r.fluidVolume, A = r.wettedArea;
    std::vector<std::pair<std::string, std::string>> rows;
    auto row = [&](const std::string& label, const char* fmt, ...) {
        char buf[512];
        va_list ap;
        va_start(ap, fmt);
        std::vsnprintf(buf, sizeof buf, fmt, ap);
        va_end(ap);
        rows.push_back({label, buf});
    };
    row("Shape volume", "%.6g mm^3", V);
    row("Solid", "%.6g mm^3, solid fraction %.4f, %zu piece(s)", Vs, Vs / V, r.solids.size());
    row("Fluid (void)", "%.6g mm^3, porosity %.4f, %zu piece(s)", Vf, Vf / V, r.fluids.size());
    row("Wetted area", "%.6g mm^2 (walls between solid and fluid)", A);
    row("Specific surface", "%.4g 1/mm = %.4g m^2/m^3 (wetted area / shape volume)", A / V, 1000 * A / V);
    row("Cut faces", "%.6g mm^2 (solid on the shape's faces)", r.cutArea);
    row("Hydraulic diameter", "%.4g mm (4 x fluid volume / wetted area)", 4 * Vf / A);
    for (int c = 0; c < r.channels && r.channels == 2; ++c)
        row(std::string("  channel ") + channelName(c), "%.6g mm^3, wetted area %.6g mm^2, hydraulic diameter %.4g mm",
            r.channelVolume[c], r.channelArea[c], 4 * r.channelVolume[c] / r.channelArea[c]);
    if (o.form == Form::Sheet)
        row("Wall thickness", "%.4g mm mean (2 x solid volume / wetted area), %.4g to %.4g mm local", 2 * Vs / A,
            r.thickness[0], r.thickness[1]);
    else
        row("Solid size", "%.4g mm (4 x solid volume / wetted area, like a strut diameter)", 4 * Vs / A);
    if (o.material[0] > 0) row("Mass", "%.4g g (solid density %g kg/m^3)", Vs * o.material[0] * 1e-6, o.material[0]);
    const auto [bending, stretching] = stiffness(r);
    char gpa[96] = "";
    if (o.material[1] > 0)
        std::snprintf(gpa, sizeof gpa, " = %.3g to %.3g GPa", bending * o.material[1], stretching * o.material[1]);
    row("Stiffness E/Es", "%.3g (bending-dominated, phi^2) to %.3g (stretching-dominated, phi/3)%s, Gibson-Ashby; "
        "%s usually lie near the %s value. For a real value, simulate a --unit-cell.",
        bending, stretching, gpa, o.form == Form::Sheet ? "sheets" : "networks", o.form == Form::Sheet ? "upper" : "lower");
    const auto flows = channelFlows(r, o);
    for (size_t k = 0; k < flows.size(); ++k) {
        const ChannelFlow& f = flows[k];
        const std::string who = flows.size() == 2 ? "channel " + f.name + ": " : "";
        row(k ? "" : "Permeability", "%s%.3g m^2 (Kozeny-Carman: porosity %.3f, Dh %.3g mm)", who.c_str(),
            f.permeability, f.porosity, f.dh);
    }
    for (size_t k = 0; k < flows.size() && o.flow[0] > 0; ++k) {
        const ChannelFlow& f = flows[k];
        const std::string who = flows.size() == 2 ? "channel " + f.name + ": " : "";
        char over[96] = "";
        if (f.axis >= 0) std::snprintf(over, sizeof over, ", %.3g Pa over %g mm along %s", f.drop, f.length, kAxisNames[f.axis]);
        row(k ? "" : "Pressure drop", "%sRe %.3g, %.3g Pa/m (%.0f%% viscous)%s", who.c_str(), f.re, f.gradient,
            100 * f.viscous / f.gradient, over);
    }
    if (o.flow[0] > 0)
        row("", "at %g m/s superficial velocity, viscosity %g Pa s, density %g kg/m^3; Kozeny-Carman and Ergun: "
                "rough packed-bed estimates, use CFD for design",
            o.flow[0], o.flow[1], o.flow[2]);
    auto list = [&](const char* label, const std::vector<ReportPiece>& ps, bool fluid) {
        for (size_t k = 0; k < ps.size() && k < 6; ++k) {
            const std::string chan = fluid && r.channels == 2 ? std::string(" (") + channelName(ps[k].channel) + ")" : "";
            row(k ? "" : label, "%.6g mm^3%s, touches %s", ps[k].volume, chan.c_str(),
                ps[k].faces.empty() ? "no face" : joined(ps[k].faces, " ").c_str());
        }
        if (ps.size() > 6) row("", "... and %zu smaller piece(s)", ps.size() - 6);
    };
    list("Solid pieces", r.solids, false);
    list("Fluid pieces", r.fluids, true);
    const auto th = through(r);
    for (size_t k = 0; k < th.size(); ++k)
        row(k ? "" : "Connects", "%s to %s: solid %s, fluid %s", r.opposite[k].first.c_str(),
            r.opposite[k].second.c_str(), th[k].solid ? "yes" : "no",
            th[k].fluid.empty() ? "no" : r.channels == 1 ? "yes" : joined(th[k].fluid, " and ").c_str());
    if (r.cavities) row("Cavities", "%zu inner shell(s), exported as voids of their solid", r.cavities);
    row("Accuracy", "from meshes at %d steps per cell: volumes and areas within about 1%%", r.res);
    if (o.shape.slabAxis >= 0 || o.shape.channelAxis[0] >= 0)
        row("Note", "for the lattice box: %s are not included",
            o.plates ? "plates" : o.buffers ? "buffers" : "channel buffers and their plates");
    if (r.thickness[0] > 0 && r.thickness[0] < r.step)
        row("Note", "walls down to %.3g mm are thinner than the mesh step (%.3g mm): piece counts and the solid "
                    "volume can be off there. Raise --res.", r.thickness[0], r.step);
    return rows;
}

void printReport(const Report& r, const Options& o, const Print& print) {
    print("Properties:");
    for (const auto& [label, value] : reportRows(r, o)) say(print, "  %-18s%s %s", label.c_str(), label.empty() ? " " : ":", value.c_str());
}

bool writeReportJson(const std::string& path, const Report& r, const Options& o, const Lattice& lat) {
    const double V = r.shapeVolume, Vs = r.solidVolume, Vf = r.fluidVolume, A = r.wettedArea;
    auto pieces = [&](const std::vector<ReportPiece>& ps, bool fluid) {
        std::string s = "[";
        for (size_t k = 0; k < ps.size(); ++k) {
            std::vector<std::string> faces;
            for (const auto& f : ps[k].faces) faces.push_back(quoted(f));
            s += std::string(k ? ", " : "") + "{\"volume_mm3\": " + num(ps[k].volume);
            if (fluid && r.channels == 2) s += ", \"channel\": " + quoted(channelName(ps[k].channel));
            s += ", \"faces\": [" + joined(faces, ", ") + "]}";
        }
        return s + "]";
    };
    const char* const kAxes[3] = {"x", "y", "z"};
    std::string level = num(lat.level);
    if (lat.grade >= 0)
        level = "{\"start\": " + num(lat.levels.front()) + ", \"end\": " + num(lat.levels.back()) +
                ", \"axis\": " + quoted(kAxes[lat.grade]) + "}";
    std::string channels = "[";
    const auto flows = channelFlows(r, o);
    for (size_t c = 0; c < flows.size(); ++c) {
        const ChannelFlow& f = flows[c];
        channels += std::string(c ? ",\n              " : "") + "{\"name\": " + quoted(f.name) + ", \"volume_mm3\": " +
                    num(f.volume) + ", \"wetted_area_mm2\": " + num(f.area) + ", \"porosity\": " + num(f.porosity) +
                    ", \"hydraulic_diameter_mm\": " + num(f.dh) + ", \"permeability_m2\": " + num(f.permeability);
        if (o.flow[0] > 0) {
            channels += ", \"reynolds\": " + num(f.re) + ", \"pressure_gradient_pa_per_m\": " + num(f.gradient) +
                        ", \"viscous_gradient_pa_per_m\": " + num(f.viscous);
            if (f.axis >= 0)
                channels += ", \"pressure_drop_pa\": " + num(f.drop) + ", \"flow_axis\": " + quoted(kAxisNames[f.axis]) +
                            ", \"flow_length_mm\": " + num(f.length);
        }
        channels += "}";
    }
    channels += "]";
    const auto [bending, stretching] = stiffness(r);
    std::string extras = "\"stiffness_ratio\": {\"bending\": " + num(bending) + ", \"stretching\": " + num(stretching) + "}";
    if (o.material[0] > 0)
        extras += ",\n  \"material\": {\"density_kg_m3\": " + num(o.material[0]) + ", \"mass_g\": " +
                  num(Vs * o.material[0] * 1e-6) +
                  (o.material[1] > 0 ? ", \"youngs_modulus_gpa\": " + num(o.material[1]) +
                                           ", \"stiffness_gpa\": {\"bending\": " + num(bending * o.material[1]) +
                                           ", \"stretching\": " + num(stretching * o.material[1]) + "}"
                                     : std::string()) +
                  "}";
    if (o.flow[0] > 0)
        extras += ",\n  \"flow\": {\"velocity_m_s\": " + num(o.flow[0]) + ", \"viscosity_pa_s\": " + num(o.flow[1]) +
                  ", \"density_kg_m3\": " + num(o.flow[2]) + "}";
    std::string connects = "[";
    const auto th = through(r);
    for (size_t k = 0; k < th.size(); ++k) {
        std::vector<std::string> fluid;
        for (const auto& f : th[k].fluid) fluid.push_back(quoted(f));
        connects += std::string(k ? ", " : "") + "{\"faces\": [" + quoted(r.opposite[k].first) + ", " +
                    quoted(r.opposite[k].second) + "], \"solid\": " + (th[k].solid ? "true" : "false") +
                    ", \"fluid\": [" + joined(fluid, ", ") + "]}";
    }
    connects += "]";
    const std::string size = o.form == Form::Sheet
                                 ? "\"wall_thickness_mm\": {\"mean\": " + num(2 * Vs / A) + ", \"min\": " +
                                       num(r.thickness[0]) + ", \"max\": " + num(r.thickness[1]) + "}"
                                 : "\"solid_size_mm\": " + num(4 * Vs / A);

    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::fprintf(f,
                 "{\n"
                 "  \"surface\": %s,\n  \"form\": %s,\n  \"shape\": %s,\n  \"level\": %s,\n"
                 "  \"mesh_steps_per_cell\": %d,\n"
                 "  \"shape_volume_mm3\": %s,\n"
                 "  \"solid\": {\"volume_mm3\": %s, \"fraction\": %s, \"pieces\": %s},\n"
                 "  \"fluid\": {\"volume_mm3\": %s, \"porosity\": %s, \"hydraulic_diameter_mm\": %s,\n"
                 "            \"channels\": %s,\n            \"pieces\": %s},\n"
                 "  \"wetted_area_mm2\": %s,\n  \"specific_surface_per_mm\": %s,\n  \"cut_area_mm2\": %s,\n"
                 "  %s,\n  %s,\n  \"connects\": %s,\n  \"cavities\": %zu\n}\n",
                 quoted(findSurface(o.type)->name).c_str(), quoted(o.form == Form::Sheet ? "sheet" : "network").c_str(),
                 quoted(shapeName(o.shape.kind)).c_str(), level.c_str(), r.res, num(V).c_str(), num(Vs).c_str(),
                 num(Vs / V).c_str(), pieces(r.solids, false).c_str(), num(Vf).c_str(), num(Vf / V).c_str(),
                 num(4 * Vf / A).c_str(), channels.c_str(), pieces(r.fluids, true).c_str(), num(A).c_str(),
                 num(A / V).c_str(), num(r.cutArea).c_str(), size.c_str(), extras.c_str(), connects.c_str(), r.cavities);
    const bool ok = !std::ferror(f);
    return std::fclose(f) == 0 && ok;
}

bool writeReportCsv(const std::string& path, const Report& r, const Options& o, const Lattice& lat) {
    const double V = r.shapeVolume, Vs = r.solidVolume, Vf = r.fluidVolume, A = r.wettedArea;
    const bool sheet = o.form == Form::Sheet, flow = o.flow[0] > 0;
    // Fixed columns, so rows of different runs line up; an empty value
    // where one does not apply. A network's fluid fills channel A's columns.
    std::vector<std::pair<std::string, std::string>> cols;
    auto text = [&](const std::string& name, const std::string& value) { cols.push_back({name, value}); };
    auto val = [&](const std::string& name, double v, bool applies = true) {
        text(name, applies && std::isfinite(v) ? num(v) : "");
    };
    text("surface", findSurface(o.type)->name);
    text("form", sheet ? "sheet" : "network");
    text("part", o.part == Part::Solid ? "solid" : "void");
    text("shape", shapeName(o.shape.kind));
    for (int a = 0; a < 3; ++a) val(std::string("cells_") + kAxisNames[a], o.cells[a]);
    for (int a = 0; a < 3; ++a) val(std::string("cell_") + kAxisNames[a] + "_mm", o.cell[a], !o.shape.bent());
    val("level", lat.grade >= 0 ? lat.levels.front() : lat.level);
    val("level_end", lat.grade >= 0 ? lat.levels.back() : lat.level);
    val("mesh_steps_per_cell", r.res);
    val("shape_volume_mm3", V);
    val("solid_volume_mm3", Vs);
    val("solid_fraction", Vs / V);
    val("fluid_volume_mm3", Vf);
    val("porosity", Vf / V);
    val("wetted_area_mm2", A);
    val("specific_surface_per_mm", A / V);
    val("cut_area_mm2", r.cutArea);
    val("hydraulic_diameter_mm", 4 * Vf / A);
    val("wall_thickness_mean_mm", 2 * Vs / A, sheet);
    val("wall_thickness_min_mm", r.thickness[0], sheet);
    val("wall_thickness_max_mm", r.thickness[1], sheet);
    val("solid_size_mm", 4 * Vs / A, !sheet);
    const auto flows = channelFlows(r, o);
    for (size_t c = 0; c < 2; ++c) {
        const bool has = c < flows.size();
        const ChannelFlow f = has ? flows[c] : ChannelFlow{};
        const std::string p = std::string("channel_") + "AB"[c] + "_";
        val(p + "volume_mm3", f.volume, has);
        val(p + "wetted_area_mm2", f.area, has);
        val(p + "porosity", f.porosity, has);
        val(p + "hydraulic_diameter_mm", f.dh, has);
        val(p + "permeability_m2", f.permeability, has);
        val(p + "reynolds", f.re, has && flow);
        val(p + "pressure_gradient_pa_per_m", f.gradient, has && flow);
        val(p + "pressure_drop_pa", f.drop, has && flow && f.axis >= 0);
    }
    const auto [bending, stretching] = stiffness(r);
    val("mass_g", Vs * o.material[0] * 1e-6, o.material[0] > 0);
    val("stiffness_ratio_bending", bending);
    val("stiffness_ratio_stretching", stretching);
    val("stiffness_bending_gpa", bending * o.material[1], o.material[1] > 0);
    val("stiffness_stretching_gpa", stretching * o.material[1], o.material[1] > 0);
    val("velocity_m_s", o.flow[0], flow);
    val("viscosity_pa_s", o.flow[1], flow);
    val("fluid_density_kg_m3", o.flow[2], flow);
    val("solid_pieces", double(r.solids.size()));
    val("fluid_pieces", double(r.fluids.size()));
    val("cavities", double(r.cavities));

    std::string head, row;
    for (size_t k = 0; k < cols.size(); ++k) head += (k ? "," : "") + cols[k].first, row += (k ? "," : "") + cols[k].second;
    std::FILE* file = std::fopen(path.c_str(), "wb");
    if (!file) return false;
    std::fprintf(file, "%s\n%s\n", head.c_str(), row.c_str());
    const bool ok = !std::ferror(file);
    return std::fclose(file) == 0 && ok;
}

}  // namespace tpms
