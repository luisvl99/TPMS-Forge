#pragma once
// The tpms command line as a function, shared by the CLI and the GUI.

#include "mesher.hpp"
#include "tpms.hpp"

#include <functional>
#include <string>
#include <vector>

namespace tpms {

// Receives output text without a trailing newline. Errors start with "error:".
using Print = std::function<void(const std::string&)>;

struct Options {
    std::string type = "gyroid";
    Form form = Form::Sheet;
    Part part = Part::Solid;
    double density = 0.2;
    double level = 0;
    bool hasLevel = false;
    double cell[3] = {10, 10, 10};  // box and trimmed shapes
    int cells[3] = {2, 2, 2};        // box and conformal shapes (along the parameter axes)
    Shape shape;
    int res = 24;
    int patches = 3;
    double tol = 0;  // 0: autoTol()
    bool faceted = false;
    bool exact = false;
    double minFragment = 0.01;
    std::vector<std::string> outputs;
};

// Parses and checks arguments (without the program name).
bool parseArgs(const std::vector<std::string>& args, Options& o, std::string& error);

// Default --tol: 0.001 x the smallest cell edge in mm.
double autoTol(const Options& o);

// The lattice the options describe; solves the iso level for --density.
Lattice makeLattice(const Options& o);

// The grid to mesh at `res` steps per cell: the box, the conformal
// parameter box (periodic around), or the trimmed shape's bounds plus one
// step of margin.
Grid gridOf(const Options& o, int res);

// Closed mesh of the exported region in mm (conformal shapes mapped).
Mesh meshSolid(const Options& o, const Lattice& lat, int res, size_t* droppedDegenerate = nullptr);

// Volume of the part's shape, mm^3.
double shapeVolume(const Options& o);

// Fraction of the shape volume inside the exported region, from
// quasi-random samples of the field. Independent of any mesh.
double sampledFraction(const Options& o, const Lattice& lat);

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
