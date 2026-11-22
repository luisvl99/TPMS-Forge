#include "export.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <unordered_map>

namespace tpms {

bool writeStl(const std::string& path, const std::vector<Mesh>& parts) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    char header[80] = "TPMS lattice";
    std::fwrite(header, 1, sizeof header, f);
    uint32_t count = 0;
    for (const Mesh& m : parts) count += uint32_t(m.triangles.size());
    std::fwrite(&count, 4, 1, f);
    for (const Mesh& m : parts)
        for (const auto& t : m.triangles) {
            const Vec3 &a = m.vertices[t[0]], &b = m.vertices[t[1]], &c = m.vertices[t[2]];
            Vec3 n = cross(b - a, c - a);
            n = n * (1.0 / std::max(norm(n), 1e-300));
            const float rec[12] = {float(n.x), float(n.y), float(n.z), float(a.x), float(a.y), float(a.z),
                                   float(b.x), float(b.y), float(b.z), float(c.x), float(c.y), float(c.z)};
            const uint16_t attr = 0;
            std::fwrite(rec, 4, 12, f);
            std::fwrite(&attr, 2, 1, f);
        }
    const bool ok = !std::ferror(f);
    return std::fclose(f) == 0 && ok;
}

// STEP REAL: always has a decimal point, e.g. "0.", "1.5", "1.E-07". Like
// "%.12G", but to_chars is several times faster (merged files hold millions).
std::string real(double v) {
    char buf[32];
    std::string s(buf, std::to_chars(buf, buf + sizeof buf, v, std::chars_format::general, 12).ptr);
    const auto e = s.find('e');
    if (e != std::string::npos) s[e] = 'E';
    if (s.find('.') == std::string::npos) s.insert(e == std::string::npos ? s.size() : e, ".");
    return s;
}

std::string ref(long id) { return "#" + std::to_string(id); }

std::string refList(const std::vector<long>& ids) {
    std::string s = "(";
    for (size_t i = 0; i < ids.size(); ++i) {
        if (i) s += i % 10 == 0 ? ",\n" : ",";
        s += ref(ids[i]);
    }
    return s + ")";
}

std::string triple(Vec3 p) { return "(" + real(p.x) + "," + real(p.y) + "," + real(p.z) + ")"; }

namespace {

std::string quoted(const std::string& s) {
    std::string out = "'";
    for (char ch : s) {
        if (ch == '\'') out += "''";
        else if (ch >= 32 && ch < 127 && ch != '\\') out += ch;
    }
    return out + "'";
}

long writeSolid(StepFile& st, const Mesh& m, const std::string& name,
                const std::function<std::string(size_t)>& triangleName) {
    std::vector<long> point(m.vertices.size()), vertex(m.vertices.size());
    for (size_t i = 0; i < m.vertices.size(); ++i) {
        point[i] = st.add("CARTESIAN_POINT(''," + triple(m.vertices[i]) + ")");
        vertex[i] = st.add("VERTEX_POINT(''," + ref(point[i]) + ")");
    }

    struct Edge {
        long curve, dir;
    };
    std::unordered_map<uint64_t, Edge> edges;
    edges.reserve(m.triangles.size() * 3 / 2 + 1);
    auto key = [](uint32_t a, uint32_t b) {
        return (uint64_t(std::min(a, b)) << 32) | std::max(a, b);
    };
    for (const auto& t : m.triangles)
        for (int e = 0; e < 3; ++e) {
            const uint32_t lo = std::min(t[e], t[(e + 1) % 3]), hi = std::max(t[e], t[(e + 1) % 3]);
            if (edges.count(key(lo, hi))) continue;
            const Vec3 d = m.vertices[hi] - m.vertices[lo];
            const double len = norm(d);
            const long dir = st.add("DIRECTION(''," + triple(d * (1.0 / len)) + ")");
            const long vec = st.add("VECTOR(''," + ref(dir) + "," + real(len) + ")");
            const long line = st.add("LINE(''," + ref(point[lo]) + "," + ref(vec) + ")");
            const long curve = st.add("EDGE_CURVE(''," + ref(vertex[lo]) + "," + ref(vertex[hi]) + "," +
                                      ref(line) + ",.T.)");
            edges.emplace(key(lo, hi), Edge{curve, dir});
        }

    // Shells: connected sets of triangles. The one with positive volume is
    // the outside; the others are inner cavities, written as voids: their
    // triangles turned around, so each shell bounds its cavity, and the solid
    // refers to it reversed.
    std::vector<uint32_t> parent(m.vertices.size());
    for (uint32_t v = 0; v < parent.size(); ++v) parent[v] = v;
    auto find = [&](uint32_t v) {
        while (parent[v] != v) v = parent[v] = parent[parent[v]];
        return v;
    };
    for (const auto& t : m.triangles) parent[find(t[1])] = find(t[0]), parent[find(t[2])] = find(t[0]);
    std::unordered_map<uint32_t, double> shellVolume;
    for (const auto& t : m.triangles)
        shellVolume[find(t[0])] += dot(m.vertices[t[0]], cross(m.vertices[t[1]], m.vertices[t[2]]));

    std::unordered_map<uint32_t, std::vector<long>> faces;
    for (size_t ti = 0; ti < m.triangles.size(); ++ti) {
        auto t = m.triangles[ti];
        const uint32_t shell = find(t[0]);
        if (shellVolume[shell] < 0) std::swap(t[1], t[2]);
        const Vec3 &a = m.vertices[t[0]], &b = m.vertices[t[1]], &c = m.vertices[t[2]];
        const Vec3 n = cross(b - a, c - a);
        const long normal = st.add("DIRECTION(''," + triple(n * (1.0 / norm(n))) + ")");
        // Reference direction: direction of the first edge, which lies in the plane.
        const long axis = st.add("AXIS2_PLACEMENT_3D(''," + ref(point[t[0]]) + "," + ref(normal) + "," +
                                 ref(edges.at(key(t[0], t[1])).dir) + ")");
        const long plane = st.add("PLANE(''," + ref(axis) + ")");
        std::string loop = "EDGE_LOOP('',(";
        for (int e = 0; e < 3; ++e) {
            const uint32_t from = t[e], to = t[(e + 1) % 3];
            const long oe = st.add("ORIENTED_EDGE('',*,*," + ref(edges.at(key(from, to)).curve) +
                                   (from < to ? ",.T.)" : ",.F.)"));
            loop += (e ? "," : "") + ref(oe);
        }
        const long loopId = st.add(loop + "))");
        const long bound = st.add("FACE_OUTER_BOUND(''," + ref(loopId) + ",.T.)");
        const std::string faceName = triangleName ? quoted(triangleName(ti)) : "''";
        faces[shell].push_back(st.add("ADVANCED_FACE(" + faceName + ",(" + ref(bound) + ")," + ref(plane) + ",.T.)"));
    }
    long outer = 0;
    std::vector<long> voids;
    for (const auto& [shell, fs] : faces) {
        const long id = st.add("CLOSED_SHELL(''," + refList(fs) + ")");
        if (shellVolume[shell] > 0) outer = id;
        else voids.push_back(st.add("ORIENTED_CLOSED_SHELL('',*," + ref(id) + ",.F.)"));
    }
    if (voids.empty()) return st.add("MANIFOLD_SOLID_BREP(" + quoted(name) + "," + ref(outer) + ")");
    return st.add("BREP_WITH_VOIDS(" + quoted(name) + "," + ref(outer) + "," + refList(voids) + ")");
}

}  // namespace

bool writeStep(const std::string& path, const std::vector<Mesh>& parts, const std::string& name,
               const std::vector<std::string>& solidNames, const TriangleName& triangleName) {
    return writeStepFile(path, name, [&](StepFile& st) {
        std::vector<long> solids;
        for (size_t k = 0; k < parts.size(); ++k) {
            std::function<std::string(size_t)> faceName;
            if (triangleName) faceName = [&, k](size_t t) { return triangleName(k, t); };
            solids.push_back(writeSolid(st, parts[k], k < solidNames.size() ? solidNames[k] : "", faceName));
        }
        return solids;
    });
}

bool writeStepFile(const std::string& path, const std::string& name,
                   const std::function<std::vector<long>(StepFile&)>& solids, double uncertainty) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::setvbuf(f, nullptr, _IOFBF, 1 << 20);

    char stamp[32];
    const std::time_t now = std::time(nullptr);
    std::strftime(stamp, sizeof stamp, "%Y-%m-%dT%H:%M:%S", std::localtime(&now));
    const std::string q = quoted(name);

    std::fprintf(f,
                 "ISO-10303-21;\nHEADER;\n"
                 "FILE_DESCRIPTION(('TPMS lattice'),'2;1');\n"
                 "FILE_NAME(%s,'%s',(''),(''),'TPMS Forge','TPMS Forge','');\n"
                 "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n"
                 "ENDSEC;\nDATA;\n",
                 q.c_str(), stamp);

    StepFile st(f);
    const long mm = st.add("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))");
    const long rad = st.add("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))");
    const long sr = st.add("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())");
    const long unc = st.add("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(" + real(uncertainty) + ")," + ref(mm) +
                            ",'distance_accuracy_value','confusion accuracy')");
    const long ctx = st.add("(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((" +
                            ref(unc) + "))GLOBAL_UNIT_ASSIGNED_CONTEXT((" + ref(mm) + "," + ref(rad) + "," +
                            ref(sr) + "))REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY'))");

    std::vector<long> items;
    const long origin = st.add("CARTESIAN_POINT('',(0.,0.,0.))");
    const long dz = st.add("DIRECTION('',(0.,0.,1.))");
    const long dx = st.add("DIRECTION('',(1.,0.,0.))");
    items.push_back(st.add("AXIS2_PLACEMENT_3D(''," + ref(origin) + "," + ref(dz) + "," + ref(dx) + ")"));
    for (long id : solids(st)) items.push_back(id);
    const long rep = st.add("ADVANCED_BREP_SHAPE_REPRESENTATION(" + q + "," + refList(items) + "," + ref(ctx) + ")");

    const long app = st.add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')");
    st.add("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000," + ref(app) + ")");
    const long pctx = st.add("PRODUCT_CONTEXT(''," + ref(app) + ",'mechanical')");
    const long product = st.add("PRODUCT(" + q + "," + q + ",''," + "(" + ref(pctx) + "))");
    st.add("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(" + ref(product) + "))");
    const long pdf = st.add("PRODUCT_DEFINITION_FORMATION('',''," + ref(product) + ")");
    const long pdctx = st.add("PRODUCT_DEFINITION_CONTEXT('part definition'," + ref(app) + ",'design')");
    const long pd = st.add("PRODUCT_DEFINITION('design',''," + ref(pdf) + "," + ref(pdctx) + ")");
    const long pds = st.add("PRODUCT_DEFINITION_SHAPE(''," + std::string("''") + "," + ref(pd) + ")");
    st.add("SHAPE_DEFINITION_REPRESENTATION(" + ref(pds) + "," + ref(rep) + ")");

    std::fputs("ENDSEC;\nEND-ISO-10303-21;\n", f);
    const bool ok = !std::ferror(f);
    return std::fclose(f) == 0 && ok;
}

}  // namespace tpms
