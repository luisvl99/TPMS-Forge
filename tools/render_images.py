"""Renders the README banner, docs/banner.png (also GitHub's social preview
image), and the program icon, src/tpms.ico, from models that build/tpms-forge makes.

Needs the project's .venv (VTK comes with cadquery-ocp) and a built tpms-forge.
Usage: .venv/Scripts/python tools/render_images.py   (.venv/bin/python elsewhere)

The models are kept as STL files in build/images; delete them to remake them.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import vtk
from PIL import Image, ImageDraw, ImageFont
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy

ROOT = Path(__file__).resolve().parent.parent
TPMS = ROOT / "build" / ("tpms-forge.exe" if os.name == "nt" else "tpms-forge")
WORK = ROOT / "build" / "images"
BOTTOM, TOP = (0.02, 0.03, 0.07), (0.10, 0.13, 0.24)  # background gradient
PALETTE = [(0.13, 0.83, 0.93), (0.39, 0.40, 0.95), (0.93, 0.28, 0.60)]  # cyan, indigo, pink
FLUIDS = [(0.22, 0.74, 0.97), (0.98, 0.45, 0.09)]  # channel A, channel B


def model(name, args):
    """A model made by tpms-forge, with its solids numbered (RegionId) and smooth normals."""
    stl = WORK / f"{name}.stl"
    if not stl.exists():
        WORK.mkdir(parents=True, exist_ok=True)
        subprocess.run([str(TPMS), *args.split(), "-o", str(stl)], check=True, stdout=subprocess.DEVNULL)
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(stl))
    solids = vtk.vtkPolyDataConnectivityFilter()  # before the normals split points at sharp edges
    solids.SetInputConnection(reader.GetOutputPort())
    solids.SetExtractionModeToAllRegions()
    solids.ColorRegionsOn()
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(solids.GetOutputPort())
    normals.SetFeatureAngle(50)
    normals.Update()
    return normals.GetOutput()


def paint(mesh, direction, palette=PALETTE, per_solid=False):
    """Point colors: a ramp through the palette along the direction, or one
    palette color per solid, darker towards the start of the direction."""
    pts = vtk_to_numpy(mesh.GetPoints().GetData())
    t = pts @ np.asarray(direction, float)
    t = (t - t.min()) / max(np.ptp(t), 1e-12)
    if per_solid:
        solid = vtk_to_numpy(mesh.GetPointData().GetArray("RegionId")) % len(palette)
        rgb = np.asarray(palette)[solid] * (0.6 + 0.4 * t)[:, None]
    else:
        stops = np.linspace(0, 1, len(palette))
        rgb = np.stack([np.interp(t, stops, [c[k] for c in palette]) for k in range(3)], axis=1)
    colors = numpy_to_vtk((rgb * 255).astype(np.uint8), deep=1)
    colors.SetName("colors")
    mesh.GetPointData().SetScalars(colors)


def actor(mesh, place, size, turn):
    """The mesh scaled so its largest extent is size, turned (degrees about
    x, then y, then z) about its center, and centered at place."""
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(mesh)
    mapper.SetColorModeToDirectScalars()
    a = vtk.vtkActor()
    a.SetMapper(mapper)
    b = mesh.GetBounds()
    center = [(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, (b[4] + b[5]) / 2]
    a.SetOrigin(center)
    a.SetPosition(np.subtract(place, center))
    a.SetScale(size / max(b[1] - b[0], b[3] - b[2], b[5] - b[4]))
    a.RotateX(turn[0])
    a.RotateY(turn[1])
    a.RotateZ(turn[2])
    p = a.GetProperty()
    p.SetAmbient(0.12)
    p.SetDiffuse(0.85)
    p.SetSpecular(0.35)
    p.SetSpecularPower(30)
    return a


def render(actors, size, eye, target, angle, background=True, marks=()):
    """An RGBA picture of the actors (transparent without background), and
    the pixels (from the top left) of the world points in marks."""
    ren = vtk.vtkRenderer()
    for a in actors:
        ren.AddActor(a)
    if background:
        ren.GradientBackgroundOn()
        ren.SetBackground(*BOTTOM)
        ren.SetBackground2(*TOP)
    else:
        ren.SetBackgroundAlpha(0)
    vtk.vtkLightKit().AddLightsToRenderer(ren)
    occlusion = vtk.vtkSSAOPass()  # darkens the creases and the channels
    occlusion.SetRadius(0.12)
    occlusion.SetBias(0.004)
    occlusion.SetKernelSize(128)
    occlusion.BlurOn()
    occlusion.SetDelegatePass(vtk.vtkRenderStepsPass())
    ren.SetPass(occlusion)
    win = vtk.vtkRenderWindow()
    win.SetOffScreenRendering(1)
    win.SetAlphaBitPlanes(1)
    win.SetSize(*size)
    win.AddRenderer(ren)
    cam = ren.GetActiveCamera()
    cam.SetPosition(*eye)
    cam.SetFocalPoint(*target)
    cam.SetViewUp(0, 0, 1)
    cam.SetViewAngle(angle)
    ren.ResetCameraClippingRange()
    win.Render()
    grab = vtk.vtkWindowToImageFilter()
    grab.SetInput(win)
    grab.SetInputBufferTypeToRGBA()
    grab.ReadFrontBufferOff()
    grab.Update()
    w, h, _ = grab.GetOutput().GetDimensions()
    px = vtk_to_numpy(grab.GetOutput().GetPointData().GetScalars()).reshape(h, w, 4)[::-1]
    pixels = []
    for p in marks:
        ren.SetWorldPoint(*p, 1)
        ren.WorldToDisplay()
        x, y, _ = ren.GetDisplayPoint()
        pixels.append((x, h - y))
    return Image.fromarray(np.ascontiguousarray(px), "RGBA"), pixels


def font(names, size):
    fonts = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"
    for name in names:
        if (fonts / name).exists():
            return ImageFont.truetype(str(fonts / name), size)
    return ImageFont.load_default(size)


def banner():
    # name, tpms-forge options, paint options, label, size, turn
    models = [
        ("fluids", "--part void --density 0.45 --res 32", dict(direction=(0, 0, 1), palette=FLUIDS, per_solid=True),
         "two fluid channels", 1.75, (0, 0, -30)),
        ("sphere", "--shape sphere-shell --radius 15 25 --cells 1 12 5 --res 32", dict(direction=(0.3, 0, 1)),
         "sphere shell", 2.25, (0, 0, 0)),
        ("torus", "--shape torus-shell --radius 5 15 --major 40 --cells 1 24 6", dict(direction=(1, 0.4, 0.2)),
         "torus shell", 2.3, (62, 0, 0)),
        ("graded", "--density 0.1 0.4 --grade x --cells 4 2 2 --res 32", dict(direction=(1, 0, 0.2)),
         "graded density", 2.3, (0, 0, 28)),
    ]
    W, H = 2560, 1280  # twice the final size: smooth edges after scaling down
    actors, marks = [], []
    for k, (name, args, colors, label, size, turn) in enumerate(models):
        mesh = model(name, args)
        paint(mesh, **colors)
        x = (k - 1.5) * 2.6
        actors.append(actor(mesh, (x, 0, 0), size, turn))
        marks.append((x, 0, -1.45))
    img, spots = render(actors, (W, H), (0, -13.5, 4.6), (0, 0, 0.95), 25, marks=marks)
    draw = ImageDraw.Draw(img)
    title = font(["segoeuib.ttf", "DejaVuSans-Bold.ttf"], 190)
    body = font(["segoeui.ttf", "DejaVuSans.ttf"], 58)
    small = font(["segoeui.ttf", "DejaVuSans.ttf"], 42)
    draw.text((150, 50), "TPMS Forge", font=title, fill=(255, 255, 255))
    draw.text((156, 320), "Lattices of triply periodic minimal surfaces as smooth STEP solids, for CAD, CFD and FEA",
              font=body, fill=(196, 205, 230))
    for (px, py), m in zip(spots, models):
        draw.text((px, py), m[3], font=small, fill=(150, 160, 190), anchor="mt")
    img.convert("RGB").resize((W // 2, H // 2), Image.LANCZOS).save(ROOT / "docs" / "banner.png", optimize=True)


def icon():
    mesh = model("cell", "--type gyroid --cells 1 --density 0.35 --res 48")
    paint(mesh, (0.3, 0.2, 1))
    img, _ = render([actor(mesh, (0, 0, 0), 2, (0, 0, 0))], (1024, 1024), (6, -6, 5), (0, 0, 0), 22, False)
    img = img.crop(img.getbbox())
    img = img.resize((round(img.width * 212 / max(img.size)), round(img.height * 212 / max(img.size))), Image.LANCZOS)
    # On a dark rounded tile, as the banner, so it reads on light and dark backgrounds.
    rows = np.linspace(0, 1, 256)[:, None, None]
    shade = np.asarray(TOP) * (1 - rows) + np.asarray(BOTTOM) * rows
    tile = Image.fromarray((np.broadcast_to(shade, (256, 256, 3)) * 255).astype(np.uint8)).convert("RGBA")
    mask = Image.new("L", (256, 256))
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, 255, 255), radius=52, fill=255)
    tile.putalpha(mask)
    tile.alpha_composite(img, ((256 - img.width) // 2, (256 - img.height) // 2))
    tile.save(WORK / "icon.png")
    tile.save(ROOT / "src" / "tpms.ico", sizes=[(s, s) for s in (16, 24, 32, 48, 64, 128, 256)])


if __name__ == "__main__":
    banner()
    icon()
    print("wrote docs/banner.png and src/tpms.ico")
