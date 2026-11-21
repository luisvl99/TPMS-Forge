"""Read STEP files with OpenCASCADE and report validity, closure and volume.

Needs: pip install cadquery-ocp
Usage: python tools/check_step.py [--periodic] [--names] file1.step [file2.step ...]

--periodic (for --unit-cell files): also compare the faces on each pair of
opposite sides of the bounding box: their area, and their vertices moved
across the box. Periodic boundary conditions in FEA need them to match.
--names: the names OpenCASCADE's document reader finds, per solid: its name
and how many faces carry each face name.
"""
import collections
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


def names(path):
    """Solid and face names, read as a CAD reader with sub-shape names on."""
    from OCP.Interface import Interface_Static
    from OCP.OCP.collections import Sequence_TDF_Label
    from OCP.STEPCAFControl import STEPCAFControl_Reader
    from OCP.TCollection import TCollection_ExtendedString
    from OCP.TDataStd import TDataStd_Name
    from OCP.TDocStd import TDocStd_Document
    from OCP.XCAFDoc import XCAFDoc_DocumentTool

    doc = TDocStd_Document(TCollection_ExtendedString("doc"))
    reader = STEPCAFControl_Reader()
    reader.SetNameMode(True)
    Interface_Static.SetIVal_s("read.stepcaf.subshapes.name", 1)  # after the reader exists
    reader.ReadFile(path)
    reader.Transfer(doc)
    tool = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())

    def name(label):
        n = TDataStd_Name()
        return n.Get().ToExtString() if label.FindAttribute(TDataStd_Name.GetID_s(), n) else ""

    top = Sequence_TDF_Label()
    tool.GetFreeShapes(top)
    whole = tool.GetShape_s(top.Value(1))
    subs = Sequence_TDF_Label()
    tool.GetSubShapes_s(top.Value(1), subs)
    named = {whole.TShape(): name(top.Value(1))}  # a single solid takes the product name
    for i in range(1, subs.Size() + 1):
        named[tool.GetShape_s(subs.Value(i)).TShape()] = name(subs.Value(i))
    out = []
    for solid in items(whole, TopAbs_SOLID):
        count = collections.Counter(named.get(f.TShape(), "-") for f in items(solid, TopAbs_FACE))
        out.append(f"{named.get(solid.TShape(), '-')!r}: " + ", ".join(f"{k} {v}" for k, v in sorted(count.items())))
    return "\n  ".join(out)


args = [a for a in sys.argv[1:] if a not in ("--periodic", "--names")]
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
    if "--names" in sys.argv:
        print("  names:", names(path))
