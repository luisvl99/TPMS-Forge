#pragma once
// The tpms command line as a function, shared by the CLI and the GUI.

#include "mesher.hpp"
#include "tpms.hpp"

#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace tpms {

// Receives output text without a trailing newline. Errors start with "error:".
using Print = std::function<void(const std::string&)>;

struct Options {
    std::string type = "gyroid";
    Form form = Form::Sheet;
    Part part = Part::Solid;
    double density = 0.2, density1 = 0;  // density1: at the end of a graded lattice
    double level = 0, level1 = 0;
    bool hasLevel = false;
    bool graded = false;  // density or level changes linearly along axis `grade`
    int grade = 0;
    double cell[3] = {10, 10, 10};  // box and trimmed shapes
    int cells[3] = {2, 2, 2};        // box, tube and torus shell (along the parameter axes)
    Shape shape;
    int res = 24;
    int patches = 3;
    double tol = 0;  // 0: autoTol()
    bool faceted = false;
    bool exact = false;
    double minFragment = 0.01;
    std::vector<std::string> outputs;
    bool report = false;    // print the property report
    bool unitCell = false;  // box: a periodic cell (exact edges, every piece kept)
    bool plates = false, buffers = false;  // box: slabs (shape.slab) on a solid, or on the fluid
    double shift[3] = {0, 0, 0};  // offset of the lattice, in cells
};

// The --shape name of a shape kind.
const char* shapeName(ShapeKind k);

// Parses and checks arguments (without the program name).
bool parseArgs(const std::vector<std::string>& args, Options& o, std::string& error);

// Default --tol: 0.001 x the smallest cell edge in mm.
double autoTol(const Options& o);

// The lattice the options describe; solves the iso level for --density.
Lattice makeLattice(const Options& o);

// The grid to mesh at `res` steps per cell: the box, a mapped
// parameter box (periodic around), or the trimmed shape's bounds plus one
// step of margin.
Grid gridOf(const Options& o, int res);

// The grid of gridOf, with whole steps added for plates or buffers.
Grid meshGrid(const Options& o, int res);

// Volume of the plates or buffers, mm^3.
double slabVolume(const Options& o);

// Closed mesh of the exported region in mm (mapped shapes mapped to mm),
// with plates or buffers.
Mesh meshSolid(const Options& o, const Lattice& lat, int res, size_t* droppedDegenerate = nullptr);

// Volume of the part's shape, mm^3.
double shapeVolume(const Options& o);

// Fraction of the shape volume inside the exported region, from
// quasi-random samples of the field. Independent of any mesh.
double sampledFraction(const Options& o, const Lattice& lat);

// Properties of the solid and the fluid (void), from closed meshes of both
// at `res` steps per cell: volumes, wetted area, hydraulic diameters, wall
// thickness, and which faces of the shape each piece touches.
struct ReportPiece {
    double volume = 0;               // mm^3
    int channel = 0;                 // fluid of a sheet: 0 = side f > c, 1 = side f < -c
    std::vector<std::string> faces;  // faces of the shape it touches
};
struct Report {
    int res = 0;
    double shapeVolume = 0, solidVolume = 0, fluidVolume = 0;  // mm^3
    double wettedArea = 0, cutArea = 0;  // mm^2: walls between solid and fluid; solid on the shape's faces
    int channels = 1;                    // fluid channels: 2 for a sheet (one on each side), else 1
    double channelVolume[2] = {0, 0}, channelArea[2] = {0, 0};
    double thickness[2] = {0, 0};        // sheet: local wall thickness, min and max (mm)
    double step = 0;                     // typical mesh grid step, mm
    size_t cavities = 0;                 // inner cavity shells (solid and fluid)
    std::vector<ReportPiece> solids, fluids;  // largest first
    std::vector<std::pair<std::string, std::string>> opposite;  // pairs of opposite faces of the shape
};
Report makeReport(const Options& o, const Lattice& lat, int res);
// The report as (label, value) rows, for the command line and the GUI.
std::vector<std::pair<std::string, std::string>> reportRows(const Report& r, const Options& o);
void printReport(const Report& r, const Options& o, const Print& print);
bool writeReportJson(const std::string& path, const Report& r, const Options& o, const Lattice& lat);

// The exported pieces of a closed mesh: without cavities and fragments
// smaller than minFragment x the largest piece.
struct Pieces {
    std::vector<Mesh> kept;
    size_t fragments = 0, cavities = 0, triangles = 0;
    double volume = 0;
};
Pieces keepPieces(const Mesh& mesh, double minFragment);

// Runs a tpms command line. Returns the exit code: 0 ok, 1 a check or write
// failed, 2 usage error.
int run(const std::vector<std::string>& args, const Print& print);

}  // namespace tpms
