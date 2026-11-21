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

struct Settings {
    int type = 0;  // index into surfaces()
    int form = 0;  // sheet, network
    int part = 0;  // solid, void
    int thickness = 0;  // by density, by iso level
    float density = 0.2f, level = 0.5f;
    bool cubic = true;
    float cell[3] = {10, 10, 10};
    int cells[3] = {2, 2, 2};
    int patches = 3;
    float tol = 0;  // 0: automatic
    bool exact = false, faceted = false;
    int res = 24;
    float minFragment = 0.01f;
    char out[1024] = "lattice.step";
};

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
    a.push_back("--cell");
    for (int k = 0; k < (s.cubic ? 1 : 3); ++k) a.push_back(num(s.cell[k]));
    a.push_back("--cells");
    const bool same = s.cells[0] == s.cells[1] && s.cells[1] == s.cells[2];
    for (int k = 0; k < (same ? 1 : 3); ++k) a.push_back(std::to_string(s.cells[k]));
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

// Preview: the faceted mesh at a resolution that keeps about 2 million grid
// nodes, so it follows the sliders. Built on a worker thread.
struct Preview {
    std::vector<float> tris;  // per vertex: x y z nx ny nz
    float box[3] = {0, 0, 0};
    bool solid = true;
    std::string info, error;
};

Preview buildPreview(const std::vector<std::string>& args) {
    Preview p;
    Options o;
    if (!parseArgs(args, o, p.error)) return p;
    const double nCells = double(o.cells[0]) * o.cells[1] * o.cells[2];
    const int res = std::clamp(int(std::cbrt(2e6 / nCells)), 6, 20);
    const Lattice lat = makeLattice(o);
    const Pieces pieces = keepPieces(meshLattice(lat, o.cells, res), o.minFragment);
    double boxVol = 1;
    for (int a = 0; a < 3; ++a) p.box[a] = float(o.cells[a] * o.cell[a]), boxVol *= p.box[a];
    p.solid = o.part == Part::Solid;
    // Walls get the exact normal (gradient of the field) at each vertex, so
    // they shade smoothly; caps on the box faces stay flat.
    const double d = 1e-4 * std::min({o.cell[0], o.cell[1], o.cell[2]});
    auto onBoxFace = [&](const Vec3& a, const Vec3& b, const Vec3& c) {
        for (int k = 0; k < 3; ++k)
            for (double side : {0.0, double(p.box[k])})
                if (std::fabs(a[k] - side) < d && std::fabs(b[k] - side) < d && std::fabs(c[k] - side) < d) return true;
        return false;
    };
    for (const Mesh& m : pieces.kept) {
        std::vector<Vec3> normals;
        for (const Vec3& v : m.vertices) {
            Vec3 g;
            for (int k = 0; k < 3; ++k) {
                Vec3 q = v, r = v;
                q[k] += d, r[k] -= d;
                g[k] = lat(q.x, q.y, q.z) - lat(r.x, r.y, r.z);
            }
            normals.push_back(g * (1.0 / std::max(norm(g), 1e-30)));
        }
        for (const auto& t : m.triangles) {
            const Vec3 &a = m.vertices[t[0]], &b = m.vertices[t[1]], &c = m.vertices[t[2]];
            Vec3 flat = cross(b - a, c - a);
            flat = flat * (1.0 / std::max(norm(flat), 1e-30));
            const bool cap = onBoxFace(a, b, c);
            for (int k = 0; k < 3; ++k) {
                const Vec3& v = m.vertices[t[size_t(k)]];
                const Vec3& n = cap ? flat : normals[t[size_t(k)]];
                p.tris.insert(p.tris.end(), {float(v.x), float(v.y), float(v.z), float(n.x), float(n.y), float(n.z)});
            }
        }
    }
    if (pieces.kept.empty()) {
        p.error = "Nothing to show: no solid at this density / level.";
        return p;
    }
    const double frac = p.solid ? pieces.volume / boxVol : 1 - pieces.volume / boxVol;
    char b[256];
    std::snprintf(b, sizeof b, "Preview: %zu triangles (%d steps per cell), %zu solid(s), solid fraction %.3f",
                  pieces.triangles, res, pieces.kept.size(), frac);
    p.info = b;
    if (pieces.fragments) p.info += ", " + std::to_string(pieces.fragments) + " fragment(s) dropped";
    if (pieces.cavities) p.info += ", " + std::to_string(pieces.cavities) + " cavity shell(s) dropped";
    return p;
}

struct Camera {
    float yaw = -35, pitch = 25;  // degrees; z is up
    float dist = 0, panX = 0, panY = 0;
};

void fit(Camera& cam, const float box[3]) {
    cam = Camera{};
    cam.dist = 1.25f * 0.5f * std::sqrt(box[0] * box[0] + box[1] * box[1] + box[2] * box[2]) / std::tan(0.2618f);
}

void drawScene(const Preview& p, const Camera& cam, int x, int y, int w, int h) {
    glViewport(x, y, w, h);
    glScissor(x, y, w, h);
    glEnable(GL_SCISSOR_TEST);
    glClearColor(0.12f, 0.13f, 0.15f, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_SCISSOR_TEST);
    if (p.tris.empty() || w <= 0 || h <= 0) return;

    const float r = 0.5f * std::sqrt(p.box[0] * p.box[0] + p.box[1] * p.box[1] + p.box[2] * p.box[2]);
    const double zNear = std::max(cam.dist - 2 * r, 0.01f * r), zFar = cam.dist + 2 * r;
    const double top = zNear * std::tan(0.2618), right = top * w / h;  // 30 degree field of view
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-right, right, -top, top, zNear, zFar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    const float light[4] = {0.4f, 0.6f, 1, 0};  // eye space: a headlight from the upper right
    const float diffuse[4] = {0.8f, 0.8f, 0.8f, 1};
    glLightfv(GL_LIGHT0, GL_POSITION, light);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glTranslatef(cam.panX, cam.panY, -cam.dist);
    glRotatef(cam.pitch - 90, 1, 0, 0);
    glRotatef(cam.yaw, 0, 0, 1);
    glTranslatef(-0.5f * p.box[0], -0.5f * p.box[1], -0.5f * p.box[2]);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    const float ambient[4] = {0.22f, 0.22f, 0.22f, 1}, specular[4] = {0.2f, 0.2f, 0.2f, 1};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 40);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    if (p.solid) glColor3f(0.72f, 0.75f, 0.80f);
    else glColor3f(0.40f, 0.66f, 0.82f);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 6 * sizeof(float), p.tris.data());
    glNormalPointer(GL_FLOAT, 6 * sizeof(float), p.tris.data() + 3);
    glDrawArrays(GL_TRIANGLES, 0, GLsizei(p.tris.size() / 6));
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisable(GL_LIGHTING);

    // Box edges, and x / y / z axes in red / green / blue at the origin.
    glLineWidth(1.5f);
    glBegin(GL_LINES);
    glColor3f(0.55f, 0.58f, 0.62f);
    for (int a = 0; a < 3; ++a)
        for (int k = 0; k < 4; ++k) {
            float c0[3], c1[3];
            const int u = (a + 1) % 3, v = (a + 2) % 3;
            c0[a] = 0, c1[a] = p.box[a];
            c0[u] = c1[u] = (k & 1) ? p.box[u] : 0;
            c0[v] = c1[v] = (k & 2) ? p.box[v] : 0;
            glVertex3fv(c0);
            glVertex3fv(c1);
        }
    glEnd();
    glDisable(GL_DEPTH_TEST);  // axes stay visible through the lattice
    glLineWidth(2.5f);
    glBegin(GL_LINES);
    const float axis = 0.3f * std::max({p.box[0], p.box[1], p.box[2]});
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

}  // namespace

int main(int argc, char** argv) {
    // --screenshot FILE.ppm: render once the preview is ready, save, quit (docs and tests).
    const char* screenshot = argc == 3 && std::strcmp(argv[1], "--screenshot") == 0 ? argv[2] : nullptr;

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

    Settings s;
    Camera cam;
    Preview shown;
    std::string shownKey, wantedKey, lastBox;
    std::future<Preview> previewJob;
    std::future<int> exportJob;
    std::mutex logMutex;
    std::string log, status;
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
            const std::string box = num(shown.box[0]) + " " + num(shown.box[1]) + " " + num(shown.box[2]);
            if (box != lastBox && !shown.tris.empty()) fit(cam, shown.box), lastBox = box;
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
        if (s.thickness == 0) ImGui::SliderFloat("Solid fraction", &s.density, 0.02f, 0.95f, "%.3f");
        else ImGui::DragFloat("Iso level c", &s.level, 0.005f, -10, 10, "%.4f");

        ImGui::SeparatorText("Box");
        ImGui::Checkbox("Cubic cell", &s.cubic);
        if (s.cubic) {
            ImGui::DragFloat("Cell size (mm)", &s.cell[0], 0.1f, 0.1f, 1000, "%.2f");
            s.cell[1] = s.cell[2] = s.cell[0];
        } else {
            ImGui::DragFloat3("Cell size (mm)", s.cell, 0.1f, 0.1f, 1000, "%.2f");
        }
        ImGui::DragInt3("Cells x y z", s.cells, 0.05f, 1, 50);
        for (int k = 0; k < 3; ++k) s.cells[k] = std::max(1, s.cells[k]), s.cell[k] = std::max(0.1f, s.cell[k]);
        ImGui::TextDisabled("Box %g x %g x %g mm", s.cells[0] * s.cell[0], s.cells[1] * s.cell[1],
                            s.cells[2] * s.cell[2]);
        ImGui::DragFloat("Min fragment", &s.minFragment, 0.001f, 0, 0.5f, "%.3f");
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Drop pieces smaller than this fraction of the largest piece (cut off by the box).");

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
        ImGui::BeginDisabled(exporting || !shown.error.empty());
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
            if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) fit(cam, shown.box);
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
        drawScene(shown, cam, panelPx, 0, fbW - panelPx, fbH);
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
