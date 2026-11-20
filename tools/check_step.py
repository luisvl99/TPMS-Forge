"""Read STEP files with OpenCASCADE and report validity, closure and volume.

Needs: pip install cadquery-ocp
Usage: python tools/check_step.py file1.step [file2.step ...]
"""
import sys
import time

from OCP.BRep import BRep_Tool
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS


def items(shape, kind):
    exp = TopExp_Explorer(shape, kind)
    while exp.More():
        yield exp.Current()
        exp.Next()


for path in sys.argv[1:]:
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
