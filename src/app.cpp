// tpms-gui: every tpms option as a widget, a live 3D preview, and export on
// a worker thread. The GUI builds a tpms command line and runs it through the
// same run() as the CLI, so both always offer the same options.

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <commdlg.h>
#endif

#include "run.hpp"

#include <GLFW/glfw3.h>  // also includes GL/gl.h: fixed-function OpenGL 1.1, no loader needed

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <future>
#include <mutex>
#include <string>
#include <vector>

using namespace tpms;

namespace {

// Shapes in the order of the GUI list, with their --shape names.
const char* const kShapeArgs[] = {"box", "tube", "torus-shell", "cylinder", "sphere", "torus"};
const char* const kShapeLabels =
    "Box\0Tube (cells follow the shape)\0Torus shell (cells follow the shape)\0"
    "Cylinder (trimmed)\0Sphere (trimmed)\0Torus (trimmed)\0";
enum { kBox, kTube, kTorusShell, kCylinder, kSphere, kTorus };

struct Settings {
    int type = 0;  // index into surfaces()
    int form = 0;  // sheet, network
    int part = 0;  // solid, void
    int thickness = 0;  // by density, by iso level
    float density = 0.2f, level = 0.5f;
    int shape = kBox;
    bool cubic = true;
    float cell[3] = {10, 10, 10};
    int cells[3] = {2, 2, 2};
    float r0 = 15, r1 = 25, height = 30, major = 40;
    int radial = 0;  // log, linear
    int patches = 3;
    float tol = 0;  // 0: automatic
    bool exact = false, faceted = false;
    int res = 24;
    float minFragment = 0.01f;
    char out[1024] = "lattice.step";
};

// Sizes that suit each shape with 10 mm cells.
void shapeDefaults(Settings& s) {
    const int counts[6][3] = {{2, 2, 2}, {1, 12, 3}, {1, 24, 6}, {2, 2, 2}, {2, 2, 2}, {2, 2, 2}};
    for (int k = 0; k < 3; ++k) s.cells[k] = counts[s.shape][k];
    if (s.shape == kTube) s.r0 = 15, s.r1 = 25, s.height = 30;
    if (s.shape == kTorusShell) s.r0 = 5, s.r1 = 15, s.major = 40;
    if (s.shape == kCylinder) s.r1 = 15, s.height = 30;
    if (s.shape == kSphere) s.r1 = 15;
    if (s.shape == kTorus) s.r1 = 10, s.major = 30;
}

std::string num(double v) {
    char b[32];
    std::snprintf(b, sizeof b, "%g", v);
    return b;
}

bool isStl(const char* path) {
    const std::string e = std::filesystem::path(path).extension().string();
    return e == ".stl" || e == ".STL";
}

// Arguments that shape the solid: everything the preview shows.
std::vector<std::string> geometryArgs(const Settings& s) {
    std::vector<std::string> a = {"--type", surfaces()[size_t(s.type)].name, "--form", s.form ? "network" : "sheet",
                                  "--part", s.part ? "void" : "solid"};
    if (s.thickness) a.insert(a.end(), {"--level", num(s.level)});
    else a.insert(a.end(), {"--density", num(s.density)});
    if (s.shape != kBox) a.insert(a.end(), {"--shape", kShapeArgs[s.shape]});
    if (s.shape == kTube || s.shape == kTorusShell) a.insert(a.end(), {"--radius", num(s.r0), num(s.r1)});
    if (s.shape >= kCylinder) a.insert(a.end(), {"--radius", num(s.r1)});
    if (s.shape == kTube || s.shape == kCylinder) a.insert(a.end(), {"--height", num(s.height)});
    if (s.shape == kTorusShell || s.shape == kTorus) a.insert(a.end(), {"--major", num(s.major)});
    if (s.shape != kTube && s.shape != kTorusShell) {
        a.push_back("--cell");
        for (int k = 0; k < (s.cubic ? 1 : 3); ++k) a.push_back(num(s.cell[k]));
    }
    if (s.shape <= kTorusShell) {
        a.push_back("--cells");
        const bool same = s.cells[0] == s.cells[1] && s.cells[1] == s.cells[2];
        for (int k = 0; k < (same ? 1 : 3); ++k) a.push_back(std::to_string(s.cells[k]));
    }
    if ((s.shape == kTube || s.shape == kTorusShell) && s.radial) a.insert(a.end(), {"--radial", "linear"});
    if (s.minFragment != 0.01f) a.insert(a.end(), {"--min-fragment", num(s.minFragment)});
    return a;
}

// The full command: geometry, export options that differ from the defaults,
// and the output file.
std::vector<std::string> exportArgs(const Settings& s) {
    std::vector<std::string> a = geometryArgs(s);
    if (isStl(s.out) || s.faceted) {
        if (s.res != 24) a.insert(a.end(), {"--res", std::to_string(s.res)});
    }
    if (!isStl(s.out)) {
        if (s.faceted) a.push_back("--faceted");
        else {
            if (s.patches != 3) a.insert(a.end(), {"--patches", std::to_string(s.patches)});
            if (s.tol > 0) a.insert(a.end(), {"--tol", num(s.tol)});
            if (s.exact) a.push_back("--exact");
        }
    }
    a.insert(a.end(), {"-o", s.out});
    return a;
}

std::string commandLine(const std::vector<std::string>& args) {
    std::string c = "tpms";
    for (const auto& a : args) c += a.find(' ') == std::string::npos ? " " + a : " \"" + a + "\"";
    return c;
}

bool saveDialog(char* path, size_t size) {
#ifdef _WIN32
    OPENFILENAMEA ofn{};
    ofn.lStructSize = sizeof ofn;
    ofn.lpstrFilter = "STEP (*.step)\0*.step;*.stp\0STL (*.stl)\0*.stl\0";
    ofn.nFilterIndex = isStl(path) ? 2 : 1;
    ofn.lpstrFile = path;
    ofn.nMaxFile = DWORD(size);
    ofn.lpstrDefExt = "step";
    ofn.Flags = OFN_OVERWRITEPROMPT | OFN_NOCHANGEDIR | OFN_PATHMUSTEXIST;
    return GetSaveFileNameA(&ofn);
#else
    (void)path, (void)size;
    return false;
#endif
}

// Preview: the walls and the cut faces, built on a worker thread.
struct Preview {
    std::vector<float> tris;  // per vertex: x y z nx ny nz, in mm
    float lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    bool solid = true, box = true;
    std::string info, error;
};

Preview buildPreview(const std::vector<std::string>& args) {
    Preview p;
    Options o;
    if (!parseArgs(args, o, p.error)) return p;
    const Shape& shape = o.shape;
    const Lattice lat = makeLattice(o);
    Lattice open = lat;  // the lattice without the trimming shape
    open.shape = Shape{};
    // At most about 600 000 grid nodes, 6 to 24 steps per cell: fast enough
    // to follow the inputs.
    const Grid unit = gridOf(o, 1);
    const double units = double(unit.steps[0]) * unit.steps[1] * unit.steps[2];
    const int res = std::clamp(int(std::cbrt(6e5 / units)), 6, 24);
    const Grid grid = gridOf(o, res);
    p.solid = o.part == Part::Solid, p.box = shape.kind == ShapeKind::Box;

    double h = 1e300;
    for (int a = 0; a < 3; ++a) h = std::min(h, grid.size[a] / grid.steps[a]);
    const double d = 1e-3 * h;
    auto gradient = [&](const Field& f, Vec3 q) {
        Vec3 g;
        for (int k = 0; k < 3; ++k) {
            Vec3 a = q, b = q;
            a[k] += d, b[k] -= d;
            g[k] = f(a.x, a.y, a.z) - f(b.x, b.y, b.z);
        }
        return g;
    };
    const Field openField = [&](double x, double y, double z) { return open(x, y, z); };
    const Field shapeField = [&](double x, double y, double z) { return shape.inside(x, y, z); };

    // A point and normal in mm.
    auto emitMm = [&](Vec3 q, Vec3 n) {
        n = n * (1.0 / std::max(norm(n), 1e-30));
        p.tris.insert(p.tris.end(), {float(q.x), float(q.y), float(q.z), float(n.x), float(n.y), float(n.z)});
    };
    // A point and normal in grid space, mapped to mm.
    auto emit = [&](Vec3 q, Vec3 n) {
        if (shape.conformal()) n = shape.normal(q, n), q = shape.map(q);
        emitMm(q, n);
    };
    // The part of a triangle where g < 0, with interpolated normals.
    auto clipEmit = [&](const Vec3 q[3], const Vec3 n[3], const double g[3]) {
        Vec3 pq[4], pn[4];
        int m = 0;
        for (int k = 0; k < 3; ++k) {
            const int l = (k + 1) % 3;
            if (g[k] < 0) pq[m] = q[k], pn[m++] = n[k];
            if ((g[k] < 0) != (g[l] < 0)) {
                const double t = g[k] / (g[k] - g[l]);
                pq[m] = q[k] + (q[l] - q[k]) * t, pn[m++] = n[k] + (n[l] - n[k]) * t;
            }
        }
        for (int k = 1; k + 1 < m; ++k) emit(pq[0], pn[0]), emit(pq[k], pn[k]), emit(pq[k + 1], pn[k + 1]);
    };
    auto onGridFace = [&](const Vec3& a, const Vec3& b, const Vec3& c) {
        for (int k = 0; k < 3; ++k)
            for (double side : {grid.origin[k], grid.origin[k] + grid.size[k]})
                if (std::fabs(a[k] - side) < d && std::fabs(b[k] - side) < d && std::fabs(c[k] - side) < d) return true;
        return false;
    };

    // Walls: each level set f = c meshed alone, as the boundary of f < c.
    // A sheet's two walls can be closer than a grid step, but each wall alone
    // is smooth, so the preview has no holes at any thickness. Normals are the
    // gradient of the exported field: smooth, and pointing out of the solid.
    std::vector<double> levels = {lat.level};
    if (lat.form == Form::Sheet) levels.push_back(-lat.level);
    for (double c : levels) {
        Lattice wall = open;
        wall.form = Form::Network, wall.part = Part::Solid, wall.level = c;
        const Mesh m = meshField([&](double x, double y, double z) { return wall(x, y, z); }, grid);
        std::vector<Vec3> normals, mm;  // per vertex: normal in grid space; position and normal in mm
        for (const Vec3& v : m.vertices) {
            normals.push_back(gradient(openField, v));
            if (!shape.conformal()) continue;
            mm.push_back(shape.map(v));
            mm.push_back(shape.normal(v, normals.back()));
        }
        for (const auto& t : m.triangles) {
            const Vec3 q[3] = {m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]]};
            if (onGridFace(q[0], q[1], q[2])) continue;
            if (shape.conformal()) {  // map each vertex once, not once per triangle
                for (uint32_t v : t) emitMm(mm[2 * v], mm[2 * v + 1]);
                continue;
            }
            const Vec3 n[3] = {normals[t[0]], normals[t[1]], normals[t[2]]};
            double g[3] = {-1, -1, -1};
            if (shape.trimmed())
                for (int k = 0; k < 3; ++k) g[k] = shape.inside(q[k].x, q[k].y, q[k].z);
            clipEmit(q, n, g);
        }
    }

    if (shape.trimmed()) {
        // Cut faces on the shape: its surface, cut to the lattice region.
        // Each triangle is split 4 x 4 first, as thin sheets give thin bands.
        const Mesh surf = meshField(shapeField, grid);
        const int k = 4;
        for (const auto& t : surf.triangles) {
            const Vec3 &A = surf.vertices[t[0]], &B = surf.vertices[t[1]], &C = surf.vertices[t[2]];
            auto at = [&](int i, int j) { return A + (B - A) * (double(i) / k) + (C - A) * (double(j) / k); };
            for (int i = 0; i < k; ++i)
                for (int j = 0; i + j < k; ++j)
                    for (int up = 0; up < (i + j + 1 < k ? 2 : 1); ++up) {
                        const Vec3 q[3] = {up ? at(i + 1, j + 1) : at(i, j), up ? at(i, j + 1) : at(i + 1, j),
                                           up ? at(i + 1, j) : at(i, j + 1)};
                        Vec3 n[3];
                        double g[3];
                        for (int c = 0; c < 3; ++c)
                            n[c] = gradient(shapeField, q[c]), g[c] = open(q[c].x, q[c].y, q[c].z);
                        clipEmit(q, n, g);
                    }
        }
    } else {
        // Cut faces on the grid faces (box faces, or the faces of a conformal
        // parameter box, which map to cylinders, tori and flat rings), on a 2D
        // grid 4 times finer than the walls: thin sheets give thin bands.
        for (int a = 0; a < 3; ++a)
            for (int side = 0; side < 2 && !(grid.periodic >> a & 1); ++side) {
                const int u = (a + 1) % 3, v = (a + 2) % 3, nu = 4 * grid.steps[u], nv = 4 * grid.steps[v];
                Vec3 n;
                n[a] = side ? 1 : -1;
                auto at = [&](int i, int j) {
                    Vec3 q;
                    q[a] = grid.origin[a] + side * grid.size[a];
                    q[u] = grid.origin[u] + grid.size[u] * i / nu, q[v] = grid.origin[v] + grid.size[v] * j / nv;
                    return q;
                };
                std::vector<double> g(size_t(nu + 1) * size_t(nv + 1));
                for (int j = 0; j <= nv; ++j)
                    for (int i = 0; i <= nu; ++i) {
                        const Vec3 q = at(i, j);
                        g[size_t(j) * size_t(nu + 1) + size_t(i)] = open(q.x, q.y, q.z);
                    }
                const Vec3 ns[3] = {n, n, n};
                auto cell = [&](int i0, int j0, int i1, int j1, int i2, int j2) {
                    const Vec3 q[3] = {at(i0, j0), at(i1, j1), at(i2, j2)};
                    const double gq[3] = {g[size_t(j0) * size_t(nu + 1) + size_t(i0)],
                                          g[size_t(j1) * size_t(nu + 1) + size_t(i1)],
                                          g[size_t(j2) * size_t(nu + 1) + size_t(i2)]};
                    clipEmit(q, ns, gq);
                };
                for (int j = 0; j < nv; ++j)
                    for (int i = 0; i < nu; ++i) cell(i, j, i + 1, j, i + 1, j + 1), cell(i, j, i + 1, j + 1, i, j + 1);
            }
    }
    if (p.tris.empty()) {
        p.error = "Nothing to show: no solid at this density / level.";
        return p;
    }
    for (int a = 0; a < 3; ++a) p.lo[a] = 1e30f, p.hi[a] = -1e30f;
    for (size_t i = 0; i < p.tris.size(); i += 6)
        for (int a = 0; a < 3; ++a) p.lo[a] = std::min(p.lo[a], p.tris[i + a]), p.hi[a] = std::max(p.hi[a], p.tris[i + a]);

    const double frac = sampledFraction(o, lat);
    char b[200];
    std::snprintf(b, sizeof b, "Preview: %zu triangles (%d steps per cell), solid fraction %.3f", p.tris.size() / 18, res,
                  p.solid ? frac : 1 - frac);
    p.info = b;
    return p;
}

struct Camera {
    float yaw = -35, pitch = 25;  // degrees; z is up
    float dist = 0, panX = 0, panY = 0;
};

// Half the diagonal of the preview's bounds.
float radius(const Preview& p) {
    float r2 = 0;
    for (int a = 0; a < 3; ++a) r2 += 0.25f * (p.hi[a] - p.lo[a]) * (p.hi[a] - p.lo[a]);
    return std::sqrt(r2);
}

void fit(Camera& cam, const Preview& p) {
    cam = Camera{};
    cam.dist = 1.25f * radius(p) / std::tan(0.2618f);
}

// Compiles the preview triangles into a display list, which the driver keeps
// on the GPU: a million triangles then orbit smoothly. Frees the copy.
GLuint compile(Preview& p, GLuint old) {
    if (old) glDeleteLists(old, 1);
    if (p.tris.empty()) return 0;
    const GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 6 * sizeof(float), p.tris.data());
    glNormalPointer(GL_FLOAT, 6 * sizeof(float), p.tris.data() + 3);
    glDrawArrays(GL_TRIANGLES, 0, GLsizei(p.tris.size() / 6));
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
    glEndList();
    p.tris = {};
    return list;
}

void drawScene(const Preview& p, GLuint list, const Camera& cam, int x, int y, int w, int h) {
    glViewport(x, y, w, h);
    glScissor(x, y, w, h);
    glEnable(GL_SCISSOR_TEST);
    glClearColor(0.12f, 0.13f, 0.15f, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_SCISSOR_TEST);
    if (!list || w <= 0 || h <= 0) return;

    const float r = radius(p);
    const double zNear = std::max(cam.dist - 2 * r, 0.01f * r), zFar = cam.dist + 2 * r;
    const double top = zNear * std::tan(0.2618), right = top * w / h;  // 30 degree field of view
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-right, right, -top, top, zNear, zFar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    // Eye space: a key light from the upper right and a weak fill light from
    // the left, so faces turned away from the key light stay readable.
    const float key[4] = {0.4f, 0.6f, 1, 0}, keyDiffuse[4] = {0.7f, 0.7f, 0.7f, 1};
    const float fill[4] = {-0.8f, -0.2f, 0.5f, 0}, fillDiffuse[4] = {0.3f, 0.3f, 0.32f, 1};
    glLightfv(GL_LIGHT0, GL_POSITION, key);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, keyDiffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, fill);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, fillDiffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, fillDiffuse);
    glTranslatef(cam.panX, cam.panY, -cam.dist);
    glRotatef(cam.pitch - 90, 1, 0, 0);
    glRotatef(cam.yaw, 0, 0, 1);
    glTranslatef(-0.5f * (p.lo[0] + p.hi[0]), -0.5f * (p.lo[1] + p.hi[1]), -0.5f * (p.lo[2] + p.hi[2]));

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    // One-sided: normals point out of the solid, whatever the triangle winding.
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
    const float ambient[4] = {0.22f, 0.22f, 0.22f, 1}, specular[4] = {0.2f, 0.2f, 0.2f, 1};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 40);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    if (p.solid) glColor3f(0.72f, 0.75f, 0.80f);
    else glColor3f(0.40f, 0.66f, 0.82f);
    glCallList(list);
    glDisable(GL_LIGHTING);

    // Box edges (box shape only), and x / y / z axes in red / green / blue at
    // the origin.
    glLineWidth(1.5f);
    glBegin(GL_LINES);
    glColor3f(0.55f, 0.58f, 0.62f);
    for (int a = 0; a < 3 && p.box; ++a)
        for (int k = 0; k < 4; ++k) {
            float c0[3], c1[3];
            const int u = (a + 1) % 3, v = (a + 2) % 3;
            c0[a] = p.lo[a], c1[a] = p.hi[a];
            c0[u] = c1[u] = (k & 1) ? p.hi[u] : p.lo[u];
            c0[v] = c1[v] = (k & 2) ? p.hi[v] : p.lo[v];
            glVertex3fv(c0);
            glVertex3fv(c1);
        }
    glEnd();
    glDisable(GL_DEPTH_TEST);  // axes stay visible through the lattice
    glLineWidth(2.5f);
    glBegin(GL_LINES);
    const float axis = 0.6f * radius(p);
    const float colors[3][3] = {{0.9f, 0.3f, 0.3f}, {0.3f, 0.85f, 0.3f}, {0.35f, 0.5f, 1}};
    for (int a = 0; a < 3; ++a) {
        float end[3] = {0, 0, 0};
        end[a] = axis;
        glColor3fv(colors[a]);
        glVertex3f(0, 0, 0);
        glVertex3fv(end);
    }
    glEnd();
    glDisable(GL_DEPTH_TEST);
}

void savePpm(const char* path, int w, int h) {
    std::vector<unsigned char> px(size_t(w) * size_t(h) * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, px.data());
    if (std::FILE* f = std::fopen(path, "wb")) {
        std::fprintf(f, "P6\n%d %d\n255\n", w, h);
        for (int row = h - 1; row >= 0; --row) std::fwrite(&px[size_t(row) * size_t(w) * 3], 1, size_t(w) * 3, f);
        std::fclose(f);
    }
}

// Settings from tpms options, so `tpms-gui [tpms options]` opens with a
// copied command line filled in.
Settings fromOptions(const Options& o) {
    Settings s;
    for (size_t k = 0; k < surfaces().size(); ++k)
        if (&surfaces()[k] == findSurface(o.type)) s.type = int(k);
    s.form = o.form == Form::Network, s.part = o.part == Part::Void;
    s.thickness = o.hasLevel, s.density = float(o.density);
    if (o.hasLevel) s.level = float(o.level);
    for (int a = 0; a < 3; ++a) s.cell[a] = float(o.cell[a]), s.cells[a] = o.cells[a];
    s.cubic = o.cell[0] == o.cell[1] && o.cell[1] == o.cell[2];
    s.shape = int(o.shape.kind);  // the GUI list follows ShapeKind
    if (s.shape != kBox) {
        s.r0 = float(o.shape.r0), s.r1 = float(o.shape.r1), s.height = float(o.shape.height);
        s.major = float(o.shape.major), s.radial = o.shape.logRadial ? 0 : 1;
    }
    const double autoTol = 0.001 * std::min({o.cell[0], o.cell[1], o.cell[2]});
    s.tol = std::fabs(o.tol - autoTol) < 1e-12 ? 0 : float(o.tol);
    s.patches = o.patches, s.exact = o.exact, s.faceted = o.faceted, s.res = o.res;
    s.minFragment = float(o.minFragment);
    if (!o.outputs.empty()) std::snprintf(s.out, sizeof s.out, "%s", o.outputs[0].c_str());
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);
    // --screenshot FILE.ppm: render once the preview is ready, save, quit (docs and tests).
    const char* screenshot = nullptr;
    if (args.size() >= 2 && args[0] == "--screenshot") screenshot = argv[2], args.erase(args.begin(), args.begin() + 2);
    Settings s;
    Options startOptions;
    std::string startError;
    if (!args.empty() && parseArgs(args, startOptions, startError)) s = fromOptions(startOptions);

    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_SCALE_TO_MONITOR, GLFW_TRUE);
    GLFWwindow* window = glfwCreateWindow(1400, 860, "TPMS lattice", nullptr, nullptr);
    if (!window) return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = nullptr;
    float scale = 1;
    glfwGetWindowContentScale(window, &scale, nullptr);
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 0, style.WindowBorderSize = 0, style.FrameRounding = 4, style.GrabRounding = 4;
    style.ScaleAllSizes(scale);
    const char* font = "C:/Windows/Fonts/segoeui.ttf";
    if (std::filesystem::exists(font)) io.Fonts->AddFontFromFileTTF(font, 17 * scale);
    else io.FontGlobalScale = scale;
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    Camera cam;
    Preview shown;
    GLuint list = 0;  // display list of the shown preview
    std::string shownKey, wantedKey, lastBox;
    std::future<Preview> previewJob;
    std::future<int> exportJob;
    std::mutex logMutex;
    std::string log = startError.empty() ? "" : "error: " + startError + "\n", status;
    auto exportStart = std::chrono::steady_clock::now();
    size_t logSize = 0;
    int frames = 0;

    std::vector<const char*> typeNames;
    for (const Surface& sf : surfaces()) typeNames.push_back(sf.name.c_str());

    while (!glfwWindowShouldClose(window)) {
        const bool busy = previewJob.valid() || exportJob.valid();
        if (busy || screenshot) glfwWaitEventsTimeout(0.05);
        else glfwWaitEventsTimeout(0.5);

        // Jobs.
        const std::vector<std::string> geometry = geometryArgs(s);
        wantedKey = commandLine(geometry);
        if (previewJob.valid() && previewJob.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
            shown = previewJob.get();
            std::string box;  // refit the view when the bounds change
            for (int a = 0; a < 3; ++a) box += num(shown.lo[a]) + " " + num(shown.hi[a]) + " ";
            if (box != lastBox && !shown.tris.empty()) fit(cam, shown), lastBox = box;
            list = compile(shown, list);
        }
        if (!previewJob.valid() && wantedKey != shownKey) {
            shownKey = wantedKey;
            previewJob = std::async(std::launch::async, buildPreview, geometry);
        }
        if (exportJob.valid() && exportJob.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
            const int rc = exportJob.get();
            const double secs =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - exportStart).count();
            char b[128];
            std::snprintf(b, sizeof b, rc == 0 ? "Done in %.1f s." : "Failed after %.1f s: see the log.", secs);
            status = b;
        }

        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Settings panel.
        const float panelW = 380 * scale;
        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize({panelW, io.DisplaySize.y});
        ImGui::Begin("TPMS lattice", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |
                                                  ImGuiWindowFlags_NoCollapse);
        ImGui::PushItemWidth(-140 * scale);

        ImGui::SeparatorText("Surface");
        ImGui::Combo("Type", &s.type, typeNames.data(), int(typeNames.size()));
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("f = %s", surfaces()[size_t(s.type)].formula.c_str());
        ImGui::Combo("Form", &s.form, "Sheet   |f| <= c\0Network   f <= c\0");
        ImGui::Combo("Part", &s.part, "Solid\0Void (fluid domain)\0");
        ImGui::Combo("Thickness by", &s.thickness, "Solid fraction\0Iso level c\0");
        // Numbers are typed, with +/- steps; values are clamped to what tpms accepts.
        if (s.thickness == 0) ImGui::InputFloat("Solid fraction", &s.density, 0.01f, 0.05f, "%.3f");
        else ImGui::InputFloat("Iso level c", &s.level, 0.01f, 0.1f, "%.4f");
        s.density = std::clamp(s.density, 0.01f, 0.99f);

        ImGui::SeparatorText("Shape");
        if (ImGui::Combo("Shape", &s.shape, kShapeLabels)) shapeDefaults(s);
        const bool conformal = s.shape == kTube || s.shape == kTorusShell;
        if (s.shape == kTorusShell || s.shape == kTorus) ImGui::InputFloat("Ring radius (mm)", &s.major, 1, 5, "%.2f");
        if (conformal) {
            ImGui::InputFloat(s.shape == kTube ? "Inner radius (mm)" : "Inner tube radius", &s.r0, 1, 5, "%.2f");
            ImGui::InputFloat(s.shape == kTube ? "Outer radius (mm)" : "Outer tube radius", &s.r1, 1, 5, "%.2f");
        } else if (s.shape != kBox) {
            ImGui::InputFloat(s.shape == kTorus ? "Tube radius (mm)" : "Radius (mm)", &s.r1, 1, 5, "%.2f");
        }
        if (s.shape == kTube || s.shape == kCylinder) ImGui::InputFloat("Height (mm)", &s.height, 1, 5, "%.2f");
        // Keep sizes valid: positive, inner below outer, the ring outside the tube.
        s.r0 = std::max(0.1f, s.r0), s.r1 = std::max(conformal ? s.r0 + 0.1f : 0.1f, s.r1);
        s.height = std::max(0.1f, s.height), s.major = std::max(s.r1 + 0.1f, s.major);
        if (conformal) {
            ImGui::InputInt3(s.shape == kTube ? "Cells: r, around, z" : "Cells: r, ring, tube", s.cells);
            ImGui::Combo("Radial spacing", &s.radial, "Log (square cells)\0Linear (even thickness)\0");
        } else {
            ImGui::Checkbox("Cubic cell", &s.cubic);
            if (s.cubic) {
                ImGui::InputFloat("Cell size (mm)", &s.cell[0], 1, 5, "%.2f");
                s.cell[1] = s.cell[2] = s.cell[0];
            } else {
                ImGui::InputFloat3("Cell size (mm)", s.cell, "%.2f");
            }
            if (s.shape == kBox) ImGui::InputInt3("Cells x y z", s.cells);
        }
        for (int k = 0; k < 3; ++k) s.cells[k] = std::clamp(s.cells[k], 1, 200), s.cell[k] = std::max(0.1f, s.cell[k]);
        const float tau = 6.2831853f;
        if (s.shape == kBox)
            ImGui::TextDisabled("Box %g x %g x %g mm", s.cells[0] * s.cell[0], s.cells[1] * s.cell[1],
                                s.cells[2] * s.cell[2]);
        else if (s.shape == kTube)
            ImGui::TextDisabled("Cells: %.1f mm radial, %.1f-%.1f around, %.1f axial", (s.r1 - s.r0) / s.cells[0],
                                tau * s.r0 / s.cells[1], tau * s.r1 / s.cells[1], s.height / s.cells[2]);
        else if (s.shape == kTorusShell)
            ImGui::TextDisabled("Cells: %.1f mm radial, %.1f-%.1f ring, %.1f-%.1f tube",
                                (s.r1 - s.r0) / s.cells[0], tau * (s.major - s.r1) / s.cells[1],
                                tau * (s.major + s.r1) / s.cells[1], tau * s.r0 / s.cells[2], tau * s.r1 / s.cells[2]);
        ImGui::InputFloat("Min fragment", &s.minFragment, 0.005f, 0.05f, "%.3f");
        s.minFragment = std::clamp(s.minFragment, 0.0f, 1.0f);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Drop pieces smaller than this fraction of the largest piece (cut off by the shape).");

        ImGui::SeparatorText("Export");
        ImGui::InputText("##file", s.out, sizeof s.out);
#ifdef _WIN32
        ImGui::SameLine();
        if (ImGui::Button("Browse...")) saveDialog(s.out, sizeof s.out);
#endif
        if (isStl(s.out)) {
            ImGui::InputInt("Mesh steps per cell", &s.res);
        } else {
            ImGui::Checkbox("Faceted (one face per triangle)", &s.faceted);
            if (s.faceted) {
                ImGui::InputInt("Mesh steps per cell", &s.res);
            } else if (s.shape != kBox) {
                ImGui::PushTextWrapPos(0);
                ImGui::TextColored({1, 0.6f, 0.4f, 1}, "Smooth STEP works for the box only so far: tick Faceted or "
                                                       "export .stl.");
                ImGui::PopTextWrapPos();
            } else {
                ImGui::InputInt("Start patches per cell", &s.patches);
                ImGui::InputFloat("Tolerance (mm)", &s.tol, 0, 0, s.tol > 0 ? "%.4f" : "%.0f (auto)");
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("0 = auto: 0.001 x the cell size.");
                ImGui::Checkbox("Exact edges (more faces)", &s.exact);
            }
        }
        s.res = std::max(4, s.res), s.patches = std::max(2, s.patches), s.tol = std::max(0.0f, s.tol);

        const std::vector<std::string> args = exportArgs(s);
        const bool exporting = exportJob.valid();
        const bool unsupported = s.shape != kBox && !isStl(s.out) && !s.faceted;
        ImGui::BeginDisabled(exporting || !shown.error.empty() || unsupported);
        if (ImGui::Button(exporting ? "Exporting..." : "Export", {-1, 0})) {
            {
                std::lock_guard<std::mutex> lock(logMutex);
                log = "> " + commandLine(args) + "\n";
            }
            status.clear();
            exportStart = std::chrono::steady_clock::now();
            exportJob = std::async(std::launch::async, [args, &log, &logMutex] {
                return run(args, [&](const std::string& text) {
                    std::lock_guard<std::mutex> lock(logMutex);
                    log += text + "\n";
                });
            });
        }
        ImGui::EndDisabled();
        if (exporting)
            ImGui::Text("Working... %.0f s",
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - exportStart).count());
        else if (!status.empty())
            ImGui::TextColored(status[0] == 'D' ? ImVec4(0.5f, 0.9f, 0.5f, 1) : ImVec4(1, 0.5f, 0.4f, 1), "%s",
                               status.c_str());

        ImGui::SeparatorText("Command line");
        const std::string cmd = commandLine(args);
        ImGui::PushTextWrapPos(0);
        ImGui::TextDisabled("%s", cmd.c_str());
        ImGui::PopTextWrapPos();
        if (ImGui::SmallButton("Copy")) ImGui::SetClipboardText(cmd.c_str());

        ImGui::SeparatorText("Log");
        ImGui::BeginChild("log", {0, 0}, ImGuiChildFlags_Borders);
        {
            std::lock_guard<std::mutex> lock(logMutex);
            ImGui::TextUnformatted(log.c_str());
            if (log.size() != logSize) ImGui::SetScrollHereY(1), logSize = log.size();
        }
        ImGui::EndChild();
        ImGui::PopItemWidth();
        ImGui::End();

        // 3D view: drag to orbit, right / middle drag to pan, wheel to zoom,
        // double-click to reset.
        if (!io.WantCaptureMouse) {
            if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
                cam.yaw += io.MouseDelta.x * 0.4f;
                cam.pitch = std::clamp(cam.pitch + io.MouseDelta.y * 0.4f, -89.0f, 89.0f);
            }
            if (ImGui::IsMouseDragging(ImGuiMouseButton_Right) || ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
                const float k = cam.dist * 0.0015f / scale;
                cam.panX += io.MouseDelta.x * k, cam.panY -= io.MouseDelta.y * k;
            }
            if (io.MouseWheel != 0) cam.dist *= std::pow(0.88f, io.MouseWheel);
            if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) fit(cam, shown);
        }
        ImDrawList* fg = ImGui::GetForegroundDrawList();
        const ImVec2 at = {panelW + 12 * scale, 10 * scale};
        if (!shown.error.empty()) fg->AddText(at, IM_COL32(255, 130, 110, 255), shown.error.c_str());
        else fg->AddText(at, IM_COL32(200, 205, 210, 255), shown.info.c_str());
        if (previewJob.valid()) fg->AddText({at.x, at.y + 22 * scale}, IM_COL32(160, 165, 170, 255), "Updating...");
        fg->AddText({at.x, io.DisplaySize.y - 26 * scale}, IM_COL32(130, 135, 140, 255),
                    "Drag: orbit   Right drag: pan   Wheel: zoom   Double-click: reset view");

        ImGui::Render();
        int fbW, fbH;
        glfwGetFramebufferSize(window, &fbW, &fbH);
        const int panelPx = int(panelW * (io.DisplaySize.x > 0 ? fbW / io.DisplaySize.x : 1));
        drawScene(shown, list, cam, panelPx, 0, fbW - panelPx, fbH);
        glViewport(0, 0, fbW, fbH);
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        if (screenshot && ++frames > 3 && !previewJob.valid() && shownKey == wantedKey) {
            savePpm(screenshot, fbW, fbH);
            glfwSetWindowShouldClose(window, GLFW_TRUE);
        }
        glfwSwapBuffers(window);
    }

    if (previewJob.valid()) previewJob.wait();
    if (exportJob.valid()) exportJob.wait();
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

#ifdef _WIN32
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) { return main(__argc, __argv); }
#endif
