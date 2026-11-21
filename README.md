# tpms

`tpms` turns a TPMS lattice into a STEP or STL file. TPMS stands for triply periodic minimal surface. There is a command-line tool, `tpms`, and a desktop app with a live 3D preview, `tpms-gui`.

**Status: v2.** STEP output has smooth bicubic B-spline faces that follow the exact surface within a set tolerance. STL output, and STEP with `--faceted`, are triangle meshes (the v0 path). The GUI offers every command-line option. Besides boxes, the lattice can fill a tube or torus shell with cells that follow the curvature (smooth STEP too), or a cylinder, sphere or torus cut from a straight lattice (STL and faceted STEP so far).

![tpms-gui](docs/gui.png)

## GUI

Run `build/tpms-gui`.

- **Left panel:** every option of the command line, grouped as Surface, Shape and Export. Numbers are typed (with − / + buttons) and clamped to valid values. Options that do not apply to the chosen output file are hidden: mesh resolution for STEP, patches and tolerance for STL.
- **Preview:** shows the solid as you change the settings. Changes to the surface, thickness, cell size, cell count or minimum fragment rebuild it on a worker thread in about half a second. Each wall (level set) is meshed on its own, so even very thin sheets show without holes; the cuts on the box faces come from a 2D grid 4 times finer. The line above the view gives the triangle count and the solid fraction (from quasi-random samples of the field, not from the preview mesh). Pieces that the export drops as fragments are still drawn. Drag to orbit, right-drag to pan, use the wheel to zoom, and double-click to reset the view. Axes: x red, y green, z blue (z up).
- **Export:** writes the file on a worker thread. The log shows the same report as the command line.
- **Command line:** the equivalent `tpms` command, with a Copy button. Use it to repeat an export or to script one.

`tpms-gui` also takes the `tpms` options, so `tpms-gui --type diamond --density 0.3 -o d.step` (for example a copied command line) opens with those settings. `tpms-gui --screenshot FILE.ppm [options]` renders once the preview is final, saves the window as an image and quits (for docs and tests).

## Build

You need a C++17 compiler and CMake 3.20 or newer. OpenMP is used when the compiler has it. The command-line tool has no other dependencies.

`tpms-gui` uses Dear ImGui 1.91.5 and GLFW 3.4. CMake downloads them at configure time (pinned versions, SHA-256 checked). Turn the GUI off with `-DTPMS_GUI=OFF` to build only the command-line tool without downloads.

```bash
cmake -S . -B build -G Ninja
```

```bash
cmake --build build
```

With Visual Studio, leave out `-G Ninja` and build with `cmake --build build --config Release`.

## Examples

Gyroid sheet, 25 % solid, 4 x 4 x 4 cells of 8 mm, written as STEP and STL:

```bash
build/tpms --type gyroid --form sheet --density 0.25 --cell 8 --cells 4 -o gyroid.step -o gyroid.stl
```

Fluid domain for CFD (the box minus the solid) of the same lattice:

```bash
build/tpms --type gyroid --form sheet --density 0.25 --cell 8 --cells 4 --part void -o gyroid_fluid.step
```

Schwarz P network with an explicit iso level and non-cubic cells:

```bash
build/tpms --type schwarz-p --form network --level 0.3 --cell 10 10 15 --cells 3 3 2 -o p.step
```

Run `build/tpms --help` for all options. Run `build/tpms --list` for all surface types.

## Surface types

All surfaces are level-set (nodal) approximations with period 2π. The CLI maps one period onto one unit cell.

| Name | Aliases | f(x, y, z) |
|---|---|---|
| gyroid | g | sin x cos y + sin y cos z + sin z cos x |
| schwarz-p | p, primitive | cos x + cos y + cos z |
| diamond | d, schwarz-d | sin x sin y sin z + sin x cos y cos z + cos x sin y cos z + cos x cos y sin z |
| neovius | | 3 (cos x + cos y + cos z) + 4 cos x cos y cos z |
| iwp | schoen-iwp | 2 (cos x cos y + cos y cos z + cos z cos x) − (cos 2x + cos 2y + cos 2z) |
| frd | schoen-frd | 4 cos x cos y cos z − (cos 2x cos 2y + cos 2y cos 2z + cos 2z cos 2x) |
| fischer-koch-s | fks | cos 2x sin y cos z + cos x cos 2y sin z + sin x cos y cos 2z |

## Forms and parts

- `--form sheet` is the region |f| ≤ c: a wall around the surface. The wall thickness comes from the level set, so it is not exactly constant.
- `--form network` is the region f ≤ c: one of the two channel networks.
- `--density R` finds c so that one unit cell has solid volume fraction R. It always means the solid fraction, also with `--part void`.
- `--part void` exports the shape minus the solid. Use it as the fluid domain for CFD. A sheet lattice in a box has two separate fluid regions, so you get two solids.

## Shapes

| `--shape` | Size options | Cells | Kind |
|---|---|---|---|
| `box` (default) | — | `--cell L`, `--cells nx ny nz` | box `[0, n x L]` |
| `tube` | `--radius R0 R1 --height H` | `--cells radial around axial` | cells follow the shape |
| `torus-shell` | `--radius R0 R1 --major R` | `--cells radial ring tube` | cells follow the shape |
| `cylinder` | `--radius R --height H` | `--cell L` | trimmed |
| `sphere` | `--radius R` | `--cell L` | trimmed |
| `torus` | `--radius r --major R` | `--cell L` | trimmed |

- **Cells follow the shape** (tube, torus shell): the lattice is built in a parameter box with one unit per cell and bent into the shape. No cell is cut on the curved walls, and the cell count around must be whole, so the lattice closes on itself. `--radial log` (default) keeps cells square: they grow with the radius. `--radial linear` keeps the radial cell size constant, so cells stretch around on the outside. The GUI shows the resulting cell sizes.
- **Trimmed** (cylinder, sphere, torus): a straight lattice of `--cell` cells, cut by the shape. Use this for solid shapes: bending cells around an axis or a pole would squeeze them to nothing there.
- Round shapes are centered at the origin with their axis along z.
- Bending stretches cells, and stretched walls are thicker, so the solid fraction differs from `--density` (which is per unit cell). The tool samples the field and prints the real value.
- Smooth STEP works for the box, tube and torus shell. The faces where the lattice is cut are exact: cylinders and end planes for the tube, tori for the torus shell. For the trimmed shapes, use `.stl` or `--faceted` (or cut a box lattice in CAD).

Example: a gyroid tube, 15 to 25 mm radius, 30 mm long, 12 cells around:

```bash
build/tpms --shape tube --radius 15 25 --height 30 --cells 1 12 3 -o tube.step
```

## STEP options

- `--patches N` is the patch grid steps per unit cell edge to start from (default 3). Each step gives about one B-spline face per wall crossing. If the result fails a check, the tool tries N+1, N+2, … (up to 16, or N+7 if larger) and prints the value it used. FRD and Fischer-Koch S sheets need 8–11.
- `--tol T` is the max distance in mm between a face and the exact surface (default 0.001 × the smallest cell edge; for a tube or torus shell, the smallest bent cell).
- `--exact` keeps one face per patch. Neighbor faces then share their edges exactly. Without it (the default), patches are merged into large trimmed faces, which gives 3–8 times fewer faces. Their edges lie within `--tol` of the faces instead of exactly on them (tolerant edges, as in most industrial STEP files).
- `--faceted` writes the v0 STEP (one planar face per triangle) instead.

## What the tool checks

Smooth STEP runs print the patch count used, face count, solids, solid fraction and max fit error. The tool exits with an error when the closed-shell, convergence or volume check still fails at the highest patch count.

Mesh runs (STL, `--faceted`) print the triangle count, solids, whether the mesh is closed and the solid fraction.

Both warn when:

- the solid fraction of the result is more than 5 % away from what the field gives (sampled at quasi-random points). This usually means walls thinner than the grid: raise `--res` or `--patches`.
- pieces were dropped. Small fragments cut off by the shape are dropped by default (`--min-fragment`), because a floating piece breaks an FEA or CFD setup.
- inner cavities were found. The tool cannot write a solid with an inner cavity yet, so it drops the cavity shell and warns. This only happens at extreme levels, for example when a network breaks into separate blobs.

`ctest` in the build folder runs each surface type (smooth STEP) and each shape (STL) once, and fails if the tool's own checks fail.

## Checking a STEP file with OpenCASCADE

[tools/check_step.py](tools/check_step.py) reads STEP files with OpenCASCADE. It reports the solid count, closed shells, face count, `BRepCheck_Analyzer` validity, volume and read time. It needs `cadquery-ocp`:

```bash
python -m pip install cadquery-ocp
```

```bash
python tools/check_step.py lattice.step
```

Results for 2 x 2 x 2 cells at density 0.3 (sheet and network, solid and void):

- With `--exact`, all 7 surface types read as valid closed solids. The volumes match the solid fraction the tool prints within 0.0003.
- With merged faces (default), 26 of 28 cases are valid. The Neovius sheet fails: about 30 faces along the box faces get edges that OpenCASCADE splits on reading. Use `--exact` for Neovius sheets until this is fixed.

Curved shapes at density 0.3, merged faces: a tube (radius 10 to 20 mm, 10 mm long, 1 x 6 x 1 cells; all 28 cases) and a torus shell (radius 4 to 10 mm, major radius 25 mm, 1 x 12 x 4 cells; sheet solid and network void of each type):

- Everything reads as valid closed solids except the Neovius sheets, as in the box. The volumes match the solid fraction the tool prints: for the gyroid sheet tube, solid 2 817.5 mm³ + void 6 607.3 mm³ = 9 424.8 mm³ = π (20² − 10²) × 10.
- The Neovius network at this density is isolated blobs, so its void part has inner cavities. The tool drops them with a warning, as in the box.

| Faces, merged (`--exact`) | tube, 1 x 6 x 1 cells | torus shell, 1 x 12 x 4 cells |
|---|---|---|
| gyroid sheet | 354 (1 262) | 3 302 (8 582) |
| Schwarz P sheet | 189 (952) | 1 739 (6 880) |
| FRD sheet | 1 637 (10 578) | 13 354 (80 868) |

A torus shell takes 10–160 s, most of it in merging.

| Faces, 2 x 2 x 2 cells, density 0.3 | `--exact` | merged |
|---|---|---|
| gyroid sheet / network | 1 518 / 756 | 378 / 180 |
| Schwarz P sheet / network | 1 176 / 552 | 241 / 144 |
| diamond sheet / network | 1 998 / 996 | 714 / 358 |
| I-WP sheet / network | 4 398 / 642 | 524 / 162 |
| FRD sheet / network | 13 632 / 2 868 | 1 894 / 1 228 |
| Fischer-Koch S sheet / network | 2 478 / 4 932 | 2 016 / 955 |

| Default run (gyroid sheet, 2 x 2 x 2 cells) | Faces | File | OpenCASCADE read |
|---|---|---|---|
| smooth STEP, merged | 373 | 10 MB | 1.2 s |
| smooth STEP, `--exact` | 1 506 | 5 MB | 0.7 s |
| `--faceted --res 24` | 261 396 | 196 MB | 117 s |

Merged files have fewer faces but can be larger: each merged face is a dense control net (up to 60 x 60 points) over a rectangle that includes the parts outside the face.

## Notes for NX and simulation

- Face count and file size grow with the cell count. FRD and Fischer-Koch S need finer patch grids, so they give the most faces per cell.
- Merged faces are checked at the patch samples and cell mid points only. Between them a steep face can drift up to a few times `--tol`. Use a smaller `--tol` if that matters.
- For the faceted path, `--res` sets the grid steps per unit cell edge. Triangle count grows with `res²`. Thin sheets need a higher `res`; the density warning tells you when it is too low.
- Units are mm.

## How it works

### Smooth STEP (default for .step)

1. Split the solid boundary into walls. A network has one wall (f = c). A sheet has two (f = c and f = −c). Each wall is a smooth level set.
2. Build a quad mesh of each wall by dual contouring on a coarse grid (`--patches` steps per cell). Around every grid edge that the wall crosses, the cells around it give one quad. At the box, boundary squares and box edges join in, so quads end exactly on the box faces. A cell can hold several pieces of the wall. Its crossing edges are grouped through its faces, and a face with 4 crossings is split by the asymptotic decider. Both cells that share the face get the same answer, so the mesh stays manifold.
3. Project every quad corner onto the exact surface (Newton steps along the gradient). Corners on a box face stay in that face.
4. Sample each quad edge and project the samples onto the surface. Neighbor quads use the same samples for a shared edge.
5. Fill each quad with a Coons grid, project it onto the surface, and interpolate it with a bicubic B-spline. The boundary rows of the control net equal the edge curve, so neighbor faces meet exactly and no sewing is needed.
6. Check the distance of each face from the exact surface at the middle of every sample cell. Double the samples (4, 8, 16, 32) until it is below `--tol`.
7. Close the solid on each box face with planar faces. Their loops are the wall boundary curves plus the inside parts of the box edges.
8. Check the result: every edge is used by exactly two faces in opposite directions, all corners converged, and the enclosed volume matches a quasi-random sample of the field within 2 % of the shape volume. If a check fails, the patch grid is too coarse for some feature, so the tool tries the next `--patches` value.
9. Merge patches into large faces (skipped with `--exact`):
   - Group neighbor patches of one wall into regions. The tool keeps merging the pair with the narrowest combined normal cone, while all normals stay within 82° of the region's mean normal. A region is then a height field over the plane normal to that mean.
   - Fit one bicubic B-spline height field per region by least squares (banded Cholesky), adding control points until every sample is within `tol` / 2. A small second-difference penalty fills the parts of the rectangle outside the region.
   - A region that does not fit is grouped again with a 15° smaller cone. A region fails if it folds over itself, touches itself at a vertex, or has more than one outer loop. A single patch that is too curved for a height field keeps its own surface.
   - Join the patch edges between the same two faces into one B-spline curve within `tol` / 4. Chains end at kinks sharper than 30°: one curve cannot follow a corner. Zero-length edges (two corners on one point, e.g. on a mirror plane) are collapsed.
   - The STEP uncertainty is set to `tol`, so readers accept the edge gaps.
10. Write STEP AP214: `B_SPLINE_SURFACE_WITH_KNOTS` faces, `B_SPLINE_CURVE_WITH_KNOTS` edges, `PLANE` caps with `LINE` edges on the box.

For a tube or torus shell, the same steps run in the parameter box (one unit per cell):

- Corners and samples are projected onto the surface in the parameter box, then mapped to mm. The B-splines are fitted to the mapped samples, and the fit error is measured in mm.
- The faces where the parameter box wraps around (the seam) get no caps. The field repeats there, so both copies of the seam get the same patches. After fitting, the far copy of every seam vertex and edge is welded to the near one.
- Caps lie on exact `CYLINDRICAL_SURFACE`s (tube walls), `TOROIDAL_SURFACE`s (torus shell walls) or `PLANE`s (tube ends). Their edges along the shape's box edges (circles and straight lines) are B-spline curves within `tol` / 4, so the STEP uncertainty is `tol`. A curved cap keeps the seam as an edge it uses twice, as CAD seams are. A tube end is joined across the seam into an annulus with holes.
- Volumes use a flux in the parameter box: Phi(u) with dPhi/du = det J, the volume stretch of the mapping. The flux through the closed patch and cap surfaces is the volume in mm, exact up to the face samples. On the seam the flux is zero, so no cap is needed there.

### Faceted STEP and STL

1. Sample the field on a regular grid. The field is |f| − c (sheet) or f − c (network); `--part void` uses its negative. The grid covers the box; for a tube or torus shell, the parameter box (one unit per cell); for a trimmed shape, its bounds plus one grid step, with the field cut by the shape (the larger of the two values).
2. Extract the surface with marching tetrahedra. Each grid cube is split into 6 tetrahedra along its main diagonal. The split is the same in every cube, so no ambiguous cases occur and neighboring cubes share edges exactly.
3. Close the solid with planar caps on the box faces. Each face square is split along the same diagonal as the tetrahedra and clipped to the inside region. The caps reuse the surface vertices, so the mesh is watertight. Axes that wrap around (around a tube; around and along a torus shell) get no caps: their two ends are welded instead.
4. Keep each surface vertex at least 2 % of a grid edge away from grid nodes. This prevents near-zero edges and does not change the mesh topology.
5. For a tube or torus shell, map every vertex from the parameter box to mm. The mapping keeps orientation, so the mesh stays closed and outward.
6. Split the mesh into connected solids and drop small fragments.
7. Write STEP AP214 directly: `MANIFOLD_SOLID_BREP` → `CLOSED_SHELL` → `ADVANCED_FACE` on a `PLANE` per triangle, with shared `EDGE_CURVE`s. STL is binary.

## Roadmap

- **v1 (done):** smooth STEP with B-spline faces, error checked against the exact implicit function, shared boundary curves.
- **v1 merged faces (done):** 3–8 times fewer faces through trimmed height-field regions.
- **Next for v1:** test import in NX; fix merged Neovius sheets; smaller merged files (adaptive knots); inner cavities (`BREP_WITH_VOIDS`).
- **v2 (done):** GUI with every command-line option and a live preview.
- **Curved shapes, phase 1 (done):** tube and torus shell with cells that follow the shape; trimmed cylinder, sphere and torus; preview, STL and faceted STEP.
- **Curved shapes, phase 2 (done):** smooth STEP for the tube and torus shell, with exact cylinder, torus and plane caps.
- **Next:** faster merge and smaller merged files on curved shapes; graded density and cell size.
- **v3:** infill of an arbitrary part, cylindrical and spherical lattices, exact minimal surfaces (Weierstrass) as an option.
