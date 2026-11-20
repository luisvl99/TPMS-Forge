#pragma once

#include "mesher.hpp"

#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace tpms {

// Binary STL with all parts.
bool writeStl(const std::string& path, const std::vector<Mesh>& parts);

// Faceted STEP AP214. Each part becomes one MANIFOLD_SOLID_BREP with one
// planar ADVANCED_FACE per triangle. Units: mm.
bool writeStep(const std::string& path, const std::vector<Mesh>& parts, const std::string& name);

struct SmoothReport {
    size_t faces = 0, solids = 0, fragments = 0, cavities = 0;
    double volume = 0;    // mm^3 of the kept solids
    double maxError = 0;  // max distance of the faces from the exact surface, mm
    int patches = 0;      // patch grid steps per unit cell edge that was used
    int samples = 0;      // samples per patch edge
    std::string error;    // empty on success
};

// Smooth STEP AP214: one bicubic B-spline face per patch, planar caps on the
// box faces. `patches` = first patch grid steps per unit cell edge to try.
// Doubles the samples per patch until the fit error is below `tol` (mm).
SmoothReport writeSmoothStep(const std::string& path, const std::string& name, const Lattice& lat,
                             const int cells[3], int patches, double tol, double minFragment);

// Minimal ISO 10303-21 writer shared by both STEP writers.
class StepFile {
public:
    explicit StepFile(std::FILE* f) : f_(f) {}
    long add(const std::string& body) {
        std::fprintf(f_, "#%ld=%s;\n", next_, body.c_str());
        return next_++;
    }

private:
    std::FILE* f_;
    long next_ = 1;
};

std::string real(double v);  // STEP REAL, always with a decimal point
std::string ref(long id);
std::string refList(const std::vector<long>& ids);
std::string triple(Vec3 p);

// Writes header, units and product data. `solids` writes the
// MANIFOLD_SOLID_BREPs and returns their ids.
bool writeStepFile(const std::string& path, const std::string& name,
                   const std::function<std::vector<long>(StepFile&)>& solids);

}  // namespace tpms
