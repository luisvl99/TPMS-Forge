#pragma once

#include "mesher.hpp"

#include <charconv>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace tpms {

// A path with a suffix before its extension: ("a/b.step", "_A") -> "a/b_A.step".
inline std::string withSuffix(const std::string& path, const std::string& suffix) {
    std::filesystem::path p(path);
    return p.replace_filename(p.stem().string() + suffix + p.extension().string()).string();
}

// Binary STL with all parts.
bool writeStl(const std::string& path, const std::vector<Mesh>& parts);

// Faceted STEP AP214. Each part becomes one MANIFOLD_SOLID_BREP with one
// planar ADVANCED_FACE per triangle. Units: mm. Optional names: one per
// part (solid), and per triangle (part index, triangle index).
using TriangleName = std::function<std::string(size_t part, size_t triangle)>;
bool writeStep(const std::string& path, const std::vector<Mesh>& parts, const std::string& name,
               const std::vector<std::string>& solidNames = {}, const TriangleName& triangleName = {});

struct SmoothReport {
    size_t faces = 0, solids = 0, fragments = 0, cavities = 0;
    size_t patchFaces = 0;  // faces before merging
    bool merged = false;    // patches merged into large trimmed faces
    std::string warning;    // why merging was skipped
    double volume = 0;    // mm^3 of the kept solids
    double maxError = 0;  // max distance of the faces from the exact surface, mm
    int patches = 0;      // patch grid steps per unit cell edge that was used
    int samples = 0;      // samples per patch edge
    std::string error;    // empty on success
    std::vector<std::string> files;  // written: the path, or its _A and _B files (splitChannels)
    std::string note;                // why the channels stayed in one file
};

// Smooth STEP AP214 of a box, tube or torus shell: bicubic B-spline faces,
// caps on the shape's faces (planes, cylinders, tori). `patches` = first
// patch grid steps per unit cell edge to try. Doubles the samples per patch
// until the fit error is below `tol` (mm). Then merges the patches into large
// trimmed faces with edges within `tol`, unless `exact`. `fieldVolume` (mm^3,
// from samples of the field) checks the result against folded patches; it
// must match within 2% of `shapeVolume`. With `splitChannels` (a sheet's
// fluid), the solids of channel A (f > c) and B (f < -c) go to two files,
// path with _A and _B, when every solid lies in one channel.
SmoothReport writeSmoothStep(const std::string& path, const std::string& name, const Lattice& lat,
                             const int cells[3], int patches, double tol, double minFragment, bool exact,
                             double fieldVolume, double shapeVolume, bool splitChannels = false);

// Minimal ISO 10303-21 writer shared by both STEP writers.
class StepFile {
public:
    explicit StepFile(std::FILE* f) : f_(f) {}
    long add(const std::string& body) {
        // fwrite, not fprintf("%s"): about 3 times faster on big files.
        char id[24] = "#";
        char* end = std::to_chars(id + 1, id + sizeof id - 1, next_).ptr;
        *end++ = '=';
        std::fwrite(id, 1, size_t(end - id), f_);
        std::fwrite(body.data(), 1, body.size(), f_);
        std::fwrite(";\n", 1, 2, f_);
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
// MANIFOLD_SOLID_BREPs and returns their ids. `uncertainty` (mm) is the
// distance below which the reader treats points as equal.
bool writeStepFile(const std::string& path, const std::string& name,
                   const std::function<std::vector<long>(StepFile&)>& solids, double uncertainty = 1e-7);

}  // namespace tpms
