"""Read STEP files with OpenCASCADE and report validity, closure and volume.

Needs: pip install cadquery-ocp
Usage: python tools/check_step.py [--periodic] file1.step [file2.step ...]

--periodic (for --unit-cell files): also compare the faces on each pair of
opposite sides of the bounding box: their area, and their vertices moved
across the box. Periodic boundary conditions in FEA need them to match.
"""
import sys
import time

from OCP.Bnd import Bnd_Box
from OCP.BRep import BRep_Tool
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS


def items(shape, kind):
    exp = TopExp_Explorer(shape, kind)
    while exp.More():
        yield exp.Current()
        exp.Next()


def periodic(shape):
    """Per axis: the faces lying on the low and high side of the bounding box
    must have the same area and the same vertices (moved across the box)."""
    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box)
    lo, hi = box.CornerMin(), box.CornerMax()
    lo, hi = (lo.X(), lo.Y(), lo.Z()), (hi.X(), hi.Y(), hi.Z())
    size = max(h - l for l, h in zip(lo, hi))
    eps = 1e-6 * size
    out = []
    for a in range(3):
        sides = []
        for at in (lo[a], hi[a]):
            area, verts = 0.0, set()
            for f in items(shape, TopAbs_FACE):
                pts = [BRep_Tool.Pnt_s(TopoDS.Vertex(v)) for v in items(f, TopAbs_VERTEX)]
                coords = [(p.X(), p.Y(), p.Z()) for p in pts]
                if not coords or any(abs(c[a] - at) > eps for c in coords):
                    continue
                props = GProp_GProps()
                BRepGProp.SurfaceProperties_s(f, props, 1e-10)  # adaptive: the default can differ by 1e-5
                area += props.Mass()
                verts |= {tuple(round(c[k] / eps) for k in range(3) if k != a) for c in coords}
            sides.append((area, verts))
        (a0, v0), (a1, v1) = sides
        out.append(f"{'xyz'[a]}: area {a0:.4f} / {a1:.4f} mm2, vertices {'match' if v0 == v1 else 'DIFFER'}")
    return "; ".join(out)


args = [a for a in sys.argv[1:] if a != "--periodic"]
for path in args:
    t0 = time.time()
    reader = STEPControl_Reader()
    if reader.ReadFile(path) != IFSelect_RetDone:
        print(f"{path}: READ FAILED")
        continue
    reader.TransferRoots()
    shape = reader.OneShape()
    t_read = time.time() - t0

    solids = list(items(shape, TopAbs_SOLID))
    shells = list(items(shape, TopAbs_SHELL))
    faces = sum(1 for _ in items(shape, TopAbs_FACE))
    closed = sum(BRep_Tool.IsClosed_s(TopoDS.Shell(s)) for s in shells)
    t1 = time.time()
    valid = BRepCheck_Analyzer(shape).IsValid()
    t_check = time.time() - t1
    vol = 0.0
    for s in solids:
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(s, props)
        vol += props.Mass()
    print(f"{path.split('/')[-1]}: solids={len(solids)} shells={len(shells)} closed_shells={closed} "
          f"faces={faces} valid={valid} volume={vol:.4f} mm3 read={t_read:.1f}s check={t_check:.1f}s")
    if "--periodic" in sys.argv:
        print("  periodic:", periodic(shape))
