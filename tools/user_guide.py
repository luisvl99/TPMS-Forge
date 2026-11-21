"""Builds the user guide, docs/tpms-user-guide.pdf, with screenshots of tpms-gui.

Needs: python -m pip install reportlab pillow, and a built tpms-gui.
Usage: python tools/user_guide.py [--no-shots]

The screenshots go to build/guide; --no-shots reuses them. Fonts: Segoe UI
and Consolas (Windows), else DejaVu (matplotlib's copy).
"""
import os
import subprocess
import sys
from pathlib import Path

from PIL import Image as PILImage
from PIL import ImageChops
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, CondPageBreak, Frame, Image, KeepTogether, NextPageTemplate, PageBreak,
                                PageTemplate, Paragraph, Preformatted, Spacer, Table, TableStyle)
from reportlab.platypus.tableofcontents import TableOfContents

ROOT = Path(__file__).resolve().parent.parent
SHOTS = ROOT / "build" / "guide"
OUT = ROOT / "docs" / "tpms-user-guide.pdf"
GUI = ROOT / "build" / ("tpms-gui.exe" if os.name == "nt" else "tpms-gui")

# Screenshots: name, tpms options.
VIEWS = [
    ("gui", "--type gyroid --density 0.2"),
    ("t-gyroid", "--type gyroid --cells 1"),
    ("t-schwarz-p", "--type schwarz-p --cells 1"),
    ("t-diamond", "--type diamond --cells 1"),
    ("t-neovius", "--type neovius --cells 1"),
    ("t-iwp", "--type iwp --cells 1"),
    ("t-frd", "--type frd --cells 1"),
    ("t-fks", "--type fischer-koch-s --cells 1"),
    ("f-network", "--type gyroid --form network --density 0.3"),
    ("f-void", "--type gyroid --part void"),
    ("s-tube", "--shape tube --radius 15 25 --height 30 --cells 1 12 3"),
    ("s-torus-shell", "--shape torus-shell --radius 5 15 --major 40 --cells 1 24 6"),
    ("s-cylinder", "--shape cylinder --radius 15 --height 30"),
    ("s-sphere", "--shape sphere --radius 15"),
    ("s-torus", "--shape torus --radius 10 --major 30"),
    ("s-sphere-shell", "--shape sphere-shell --radius 15 25 --cells 1 12 5"),
    ("g-density", "--density 0.1 0.4 --grade x --cells 4 2 2"),
    ("g-cell", "--cell 5 --cells 4 2 2 --cell-grade x 15"),
    ("u-cell", "--unit-cell --type schwarz-p --form network --density 0.3 --shift 0.5"),
    ("p-plates", "--plates z 1.5"),
    ("p-buffers", "--part void --buffers x 10"),
    ("hx-void", "--part void --channel-buffers x y 5 1"),
    ("hx-solid", "--channel-buffers x y 5 1"),
    ("report", "--report --material 4430 114 --flow 0.1 --channel-buffers x y 5 1"),
]
BACKGROUND = (31, 33, 38)  # the preview's clear color


def take_shots():
    SHOTS.mkdir(parents=True, exist_ok=True)
    for name, args in VIEWS:
        ppm = SHOTS / f"{name}.ppm"
        subprocess.run([str(GUI), "--screenshot", str(ppm), *args.split()], check=True, timeout=300)
        PILImage.open(ppm).save(SHOTS / f"{name}.png")
        ppm.unlink()
        print("screenshot", name)


def view(name):
    """The 3D view of a screenshot, cropped to the model."""
    out = SHOTS / f"{name}-view.jpg"
    im = PILImage.open(SHOTS / f"{name}.png").convert("RGB")
    w, h = im.size
    s = w / 1400  # the window is 1400 x 860 at the display's scale
    area = im.crop((int(395 * s), int(45 * s), w, h - int(35 * s)))  # without the panel and the text lines
    differs = ImageChops.difference(area, PILImage.new("RGB", area.size, BACKGROUND)).convert("L")
    x0, y0, x1, y1 = differs.point(lambda d: 255 if d > 8 else 0).getbbox()
    pad = int(14 * s)
    box = (max(0, x0 - pad), max(0, y0 - pad), min(area.size[0], x1 + pad), min(area.size[1], y1 + pad))
    area.crop(box).save(out, quality=88)
    return out


def report_window(name):
    """The property report window of a screenshot."""
    out = SHOTS / f"{name}-window.jpg"
    im = PILImage.open(SHOTS / f"{name}.png").convert("RGB")
    s = im.size[0] / 1400
    x0, y0 = (380 + 30) * s, 50 * s
    im.crop((int(x0), int(y0), int(x0 + 620 * s), int(y0 + 460 * s))).save(out, quality=90)
    return out


def full(name):
    out = SHOTS / f"{name}-full.jpg"
    PILImage.open(SHOTS / f"{name}.png").convert("RGB").save(out, quality=88)
    return out


# ---------------------------------------------------------------- styles

def register_fonts():
    win = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"
    if (win / "segoeui.ttf").exists():
        faces = {"Body": "segoeui.ttf", "Body-Bold": "segoeuib.ttf", "Body-Italic": "segoeuii.ttf",
                 "Body-Semi": "seguisb.ttf", "Mono": "consola.ttf"}
        for name, file in faces.items():
            path = win / file
            if not path.exists():
                path = win / "segoeuib.ttf"
            pdfmetrics.registerFont(TTFont(name, str(path)))
    else:
        import matplotlib
        ttf = Path(matplotlib.__file__).parent / "mpl-data" / "fonts" / "ttf"
        for name, file in {"Body": "DejaVuSans.ttf", "Body-Bold": "DejaVuSans-Bold.ttf",
                           "Body-Italic": "DejaVuSans-Oblique.ttf", "Body-Semi": "DejaVuSans-Bold.ttf",
                           "Mono": "DejaVuSansMono.ttf"}.items():
            pdfmetrics.registerFont(TTFont(name, str(ttf / file)))
    pdfmetrics.registerFontFamily("Body", normal="Body", bold="Body-Bold", italic="Body-Italic",
                                  boldItalic="Body-Bold")


ACCENT = colors.HexColor("#1F4E79")
CODE_BG = colors.HexColor("#F2F4F7")
RULE = colors.HexColor("#C9D3DF")
HEAD_BG = colors.HexColor("#E6EDF5")
NOTE_BG = colors.HexColor("#EEF5FB")

register_fonts()
ST = {
    "title": ParagraphStyle("title", fontName="Body-Bold", fontSize=34, leading=40, textColor=ACCENT),
    "subtitle": ParagraphStyle("subtitle", fontName="Body", fontSize=15, leading=20, textColor=colors.HexColor("#33475B")),
    "h1": ParagraphStyle("h1", fontName="Body-Bold", fontSize=18, leading=23, textColor=ACCENT, spaceBefore=6, keepWithNext=1,
                         spaceAfter=8),
    "contents": ParagraphStyle("contents", fontName="Body-Bold", fontSize=18, leading=23, textColor=ACCENT,
                               spaceAfter=10),
    "h2": ParagraphStyle("h2", fontName="Body-Semi", fontSize=12.5, leading=16, textColor=ACCENT, spaceBefore=10, keepWithNext=1,
                         spaceAfter=4),
    "body": ParagraphStyle("body", fontName="Body", fontSize=9.6, leading=13.6, spaceAfter=5),
    "bullet": ParagraphStyle("bullet", fontName="Body", fontSize=9.6, leading=13.4, leftIndent=13, bulletIndent=3,
                             spaceAfter=2.5),
    "cell": ParagraphStyle("cell", fontName="Body", fontSize=8.6, leading=11.2),
    "cellhead": ParagraphStyle("cellhead", fontName="Body-Bold", fontSize=8.6, leading=11.2, textColor=ACCENT),
    "code": ParagraphStyle("code", fontName="Mono", fontSize=8.3, leading=10.6),
    "caption": ParagraphStyle("caption", fontName="Body-Italic", fontSize=8.4, leading=10.5,
                              textColor=colors.HexColor("#4A5A6A"), alignment=TA_CENTER, spaceBefore=3),
    "toc1": ParagraphStyle("toc1", fontName="Body-Semi", fontSize=10.2, leading=13, leftIndent=0, spaceBefore=3),
    "toc2": ParagraphStyle("toc2", fontName="Body", fontSize=9, leading=11.5, leftIndent=14),
}
WIDTH = A4[0] - 2 * 20 * mm


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def c(s):
    """Inline code."""
    return f'<font face="Mono" size="8.8" color="#8B1E3F">{esc(s)}</font>'


def P(text, style="body"):
    return Paragraph(text, ST[style])


def bullets(items):
    return [Paragraph(t, ST["bullet"], bulletText="•") for t in items]


def code(text):
    t = Table([[Preformatted(text.strip("\n"), ST["code"])]], colWidths=[WIDTH])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), CODE_BG), ("BOX", (0, 0), (-1, -1), 0.4, RULE),
                           ("LEFTPADDING", (0, 0), (-1, -1), 7), ("RIGHTPADDING", (0, 0), (-1, -1), 7),
                           ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
    return [t, Spacer(1, 5)]


def note(text):
    t = Table([[Paragraph(text, ST["cell"])]], colWidths=[WIDTH])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), NOTE_BG), ("LINEBEFORE", (0, 0), (0, -1), 2.2, ACCENT),
                           ("LEFTPADDING", (0, 0), (-1, -1), 8), ("TOPPADDING", (0, 0), (-1, -1), 5),
                           ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
    return [t, Spacer(1, 6)]


def table(rows, widths, head=True):
    """rows: lists of markup strings; widths: fractions of the text width."""
    data = [[Paragraph(cell, ST["cellhead" if head and i == 0 else "cell"]) for cell in row]
            for i, row in enumerate(rows)]
    t = Table(data, colWidths=[w * WIDTH for w in widths], repeatRows=1 if head else 0)
    style = [("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, -1), 0.3, RULE),
             ("TOPPADDING", (0, 0), (-1, -1), 3.2), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.2),
             ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4)]
    if head:
        style += [("BACKGROUND", (0, 0), (-1, 0), HEAD_BG), ("LINEBELOW", (0, 0), (-1, 0), 0.8, ACCENT)]
    t.setStyle(TableStyle(style))
    return [t, Spacer(1, 7)]


def picture(path, width, caption=None):
    w, h = PILImage.open(path).size
    img = Image(str(path), width=width, height=width * h / w)
    parts = [img] + ([Paragraph(caption, ST["caption"])] if caption else [])
    return KeepTogether(parts + [Spacer(1, 8)])


def gallery(items, columns, height=None):
    """items: (image path, caption). Equal columns; images fit a common height."""
    colw = WIDTH / columns
    cells = []
    for path, caption in items:
        w, h = PILImage.open(path).size
        iw = colw - 8
        ih = iw * h / w
        if height and ih > height:
            ih, iw = height, height * w / h
        cells.append([Image(str(path), width=iw, height=ih), Paragraph(caption, ST["caption"])])
    rows = []
    for i in range(0, len(cells), columns):
        row = cells[i:i + columns]
        rows.append(row + [""] * (columns - len(row)))
    t = Table(rows, colWidths=[colw] * columns)
    t.setStyle(TableStyle([("ALIGN", (0, 0), (-1, -1), "CENTER"), ("VALIGN", (0, 0), (-1, -1), "BOTTOM"),
                           ("BOTTOMPADDING", (0, 0), (-1, -1), 8)]))
    return [t, Spacer(1, 4)]


# ---------------------------------------------------------------- document

class Guide(BaseDocTemplate):
    def __init__(self, path):
        super().__init__(str(path), pagesize=A4, leftMargin=20 * mm, rightMargin=20 * mm, topMargin=20 * mm,
                         bottomMargin=20 * mm, title="tpms user guide",
                         subject="TPMS lattices as smooth STEP and STL", creator="tools/user_guide.py")
        frame = Frame(self.leftMargin, self.bottomMargin, self.width, self.height, id="body")
        self.addPageTemplates([PageTemplate("cover", [frame]),
                               PageTemplate("page", [frame], onPage=self.decorate)])

    @staticmethod
    def decorate(canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.5)
        y = 12 * mm
        canvas.line(20 * mm, y + 4 * mm, A4[0] - 20 * mm, y + 4 * mm)
        canvas.setFont("Body", 8)
        canvas.setFillColor(colors.HexColor("#5A6B7C"))
        canvas.drawString(20 * mm, y, "tpms user guide")
        canvas.drawRightString(A4[0] - 20 * mm, y, str(doc.page))
        canvas.restoreState()

    def afterFlowable(self, flowable):
        """Headings go to the contents and the PDF outline."""
        if isinstance(flowable, Paragraph) and flowable.style.name in ("h1", "h2"):
            level = 0 if flowable.style.name == "h1" else 1
            text = flowable.getPlainText()
            key = "h-" + text
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(text, key, level=level, closed=level > 0)
            self.notify("TOCEntry", (level, text, self.page, key))


NUMBER = [0, 0]  # current section and subsection


def H1(text):
    NUMBER[0], NUMBER[1] = NUMBER[0] + 1, 0
    return [PageBreak(), Paragraph(f"{NUMBER[0]}&nbsp;&nbsp;{text}", ST["h1"])]


def H2(text):
    """A subsection heading; it starts a new page if less than 45 mm are left."""
    NUMBER[1] += 1
    return [CondPageBreak(45 * mm), Paragraph(f"{NUMBER[0]}.{NUMBER[1]}&nbsp;&nbsp;{text}", ST["h2"])]


def build():
    v = {name: view(name) for name, _ in VIEWS if name not in ("gui", "report")}
    s = []

    # Cover
    s += [Spacer(1, 30 * mm), P("tpms", "title"), Spacer(1, 4),
          P("TPMS lattices as smooth STEP and STL<br/>for CAD, CFD and FEA", "subtitle"), Spacer(1, 4),
          P('<font color="#5A6B7C">User guide · November 2026</font>'), Spacer(1, 14 * mm)]
    s += gallery([(v["hx-void"], "Two-fluid heat exchanger: each channel with its own inlet and outlet"),
                  (v["s-torus-shell"], "Gyroid sheet bent into a torus shell")], 2, height=70 * mm)
    s += [Spacer(1, 10 * mm),
          P("<b>tpms</b> turns triply periodic minimal surface (TPMS) lattices into CAD files. Its STEP files "
            "have smooth B-spline faces that follow the exact surface within a set tolerance, so CAD systems "
            "such as Siemens NX import them as valid closed solids, ready for meshing. It also writes STL "
            "and faceted STEP meshes. A command-line tool, <b>tpms</b>, and a desktop app with a live 3D "
            "preview, <b>tpms-gui</b>, offer the same options."),
          NextPageTemplate("page"), PageBreak()]

    # Contents
    toc = TableOfContents()
    toc.levelStyles = [ST["toc1"], ST["toc2"]]
    toc.dotsMinLevel = 0
    toc.tableStyle = TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                                 ("RIGHTPADDING", (0, 0), (-1, -1), 0), ("TOPPADDING", (0, 0), (-1, -1), 0.5),
                                 ("BOTTOMPADDING", (0, 0), (-1, -1), 0.5)])
    s += [Paragraph("Contents", ST["contents"]), toc]

    # 1 Overview
    s += H1("What tpms does")
    s.append(P("A TPMS lattice is built from a periodic surface f(x, y, z) = 0, repeated cell by cell. tpms "
               "makes the solid around that surface (a sheet) or on one side of it (a network), fills a box or "
               "a curved shape with it, and writes it as a file that CAD, CFD and FEA tools read directly."))
    s += bullets([
        "<b>7 surface types</b>: gyroid, Schwarz P, diamond, Neovius, I-WP, FRD and Fischer-Koch S.",
        "<b>Sheet or network</b>, exported as the <b>solid</b> or as the <b>fluid</b> (the void around it). "
        "A sheet's fluid is two separate channels, written one file each.",
        "<b>Thickness</b> by solid fraction, by iso level or as a <b>wall thickness in mm</b>, constant or "
        "<b>graded</b> along an axis; <b>graded cell size</b> in a box.",
        "<b>Shapes</b>: box; tube, torus shell and sphere shell with cells that follow the curvature; cylinder, "
        "sphere and torus cut from a straight lattice.",
        "<b>Smooth STEP</b> with few, large faces (or one face per patch with exact edges), <b>faceted STEP</b> "
        "and binary <b>STL</b>. Solids and faces carry <b>names</b> (walls, inlets, outlets).",
        "<b>Periodic unit cells</b> for homogenization, with matching opposite faces.",
        "<b>Plates</b> for sandwich panels, <b>buffers</b> as CFD inlet and outlet zones, and <b>channel "
        "buffers</b> that give each fluid of a heat exchanger its own inlet and outlet.",
        "<b>Property report</b>: volumes, wetted area, hydraulic diameters, wall thickness, connectivity, "
        "mass, stiffness bounds, permeability and pressure drop estimates; as text, JSON or CSV.",
        "<b>Self-checks</b> on every file: closed solids, fit error against the exact surface, and the volume "
        "against the field.",
    ])
    s += H2("Typical uses")
    s += table([
        ["Task", "What to export", "Section"],
        ["FEA of a lattice part or panel", "the solid as smooth STEP; plates for sandwich panels", "9.1"],
        ["CFD of flow through a lattice", "the fluid (void) as smooth STEP, with inlet and outlet buffers", "9.2"],
        ["Two-fluid heat exchanger (conjugate heat transfer)", "channel buffers: the solid and both fluids, "
         "three files that fill the domain exactly", "9.3"],
        ["Effective stiffness or conductivity", "a periodic unit cell", "9.4"],
        ["3D printing", "STL", "9.5"],
        ["Comparing designs", "the property report as CSV, one row per run", "9.6"],
    ], [0.36, 0.52, 0.12])
    s += H2("Units and coordinates")
    s += bullets([
        "All lengths are in mm.",
        f"A box lattice fills [0, n × L] along each axis ({c('--cells')} n, {c('--cell')} L).",
        "Round shapes are centered at the origin with their axis along z.",
    ])

    # 2 Getting started
    s += H1("Getting started")
    s += H2("Build")
    s.append(P("You need a C++17 compiler and CMake 3.20 or newer. OpenMP is used when the compiler has it, so "
               "fitting and meshing run on all cores. The GUI downloads Dear ImGui and GLFW at configure time "
               "(pinned versions, checked by SHA-256)."))
    s += code("cmake -S . -B build -G Ninja\ncmake --build build")
    s += bullets([
        f"The programs are {c('build/tpms')} and {c('build/tpms-gui')} ({c('.exe')} on Windows). With MinGW they "
        "are linked statically and run without extra DLLs.",
        f"With Visual Studio, leave out {c('-G Ninja')} and build with {c('cmake --build build --config Release')}.",
        f"{c('-DTPMS_GUI=OFF')} builds only the command-line tool, without downloads.",
        f"{c('ctest')} in the build folder runs the test suite.",
    ])
    s += H2("Your first lattice in the GUI")
    s += bullets([
        f"Start {c('tpms-gui')}. It opens with a gyroid sheet, 20 % solid, 2 × 2 × 2 cells of 10 mm.",
        "Change the surface, thickness, shape and cells in the left panel; the preview follows in about half a "
        "second.",
        f"Type an output file name under Export ({c('.step')} or {c('.stl')}), or pick one with Browse.",
        "Press <b>Export</b>. The log shows the checks and the file written.",
        "Press <b>Property report</b> for the volumes, areas and estimates in a window of their own.",
    ])
    s += H2("Your first lattice on the command line")
    s.append(P(f"Without other options, {c('tpms')} writes the same default lattice:"))
    s += code("""tpms -o lattice.step

Surface   : gyroid, sheet, solid part
Level c   : 0.30973 (from target solid fraction 0.200)
Box       : 20 x 20 x 20 mm (2 x 2 x 2 cells of 10 x 10 x 10 mm)
Smooth    : 4 patches per cell, 371 faces (merged from 1518), 1 solid(s), solid fraction 0.1998
Fit error : 0.00928 mm max (tolerance 0.01 mm, 4 samples per patch edge)
note: dropped 2 small fragment(s) or filled tiny cavities (< 0.01 x largest).
Wrote     : lattice.step (4.8 MB)
Time      : 0.12 s""")
    s.append(P("The lines tell you what was built and that it passed the checks: the iso level found for the "
               "target solid fraction, the patch grid used, the face count, the solid fraction of the result "
               "and the largest distance between the faces and the exact surface. Small pieces cut off at "
               f"the box faces are dropped, because a floating piece breaks an FEA or CFD setup "
               f"({c('--min-fragment')}, section 7.5). The exit code is 0 when everything passed, 1 when a "
               "check or a write failed, and 2 for a usage error."))
    s.append(P(f"{c('tpms --help')} lists every option; {c('tpms --list')} lists the surface types. Section 10 "
               f"has the full reference. The examples write {c('tpms')} for {c('build/tpms')}: add the build "
               "folder to the PATH, or write the path."))

    # 3 GUI
    s += H1("The desktop app")
    s.append(picture(full("gui"), WIDTH, "tpms-gui: settings on the left, the live preview on the right"))
    s += table([
        ["Area", "What it does"],
        ["Surface", "type, form (sheet or network), part (solid or void), thickness by solid fraction or iso "
                    "level, grading"],
        ["Shape", "box, tube, torus shell, cylinder, sphere or torus with their sizes; cell size and counts; "
                  "unit cell, graded cell size, plates or buffers, channel buffers; offset and minimum fragment"],
        ["Export", "output file, and the options that apply to it: mesh steps for STL and faceted STEP; start "
                   "patches, tolerance and exact edges for smooth STEP"],
        ["Export button", "writes the file on a worker thread; the log shows the same lines as the command line"],
        ["Property report", "computes the report and shows it in a window, with Copy, Save JSON and Save CSV; "
                            "<i>Report inputs</i> holds the optional material and flow values"],
        ["Command line", "the equivalent tpms command, with a Copy button: repeat an export, or script it"],
        ["Preview", "drag to orbit, right-drag to pan, wheel to zoom, double-click to reset. Axes: x red, "
                    "y green, z blue (up). The solid is gray; the fluid is blue, and a sheet's second fluid "
                    "channel red"],
    ], [0.2, 0.8])
    s += bullets([
        "Numbers are typed (with − / + buttons) and kept within valid values. Options that exclude each other "
        "are greyed out; hover over an option for a short explanation.",
        "The line above the preview gives the triangle count and the solid fraction, sampled from the field.",
        f"{c('tpms-gui')} takes the same options as {c('tpms')}: {c('tpms-gui --type diamond --density 0.3')} "
        f"opens with those settings, so a copied command line opens where it left off. Add {c('--report')} to "
        "open with the report.",
    ])

    # 4 Lattice
    s += H1("Designing the lattice")
    s += H2("Surface types")
    s.append(P(f"Choose with {c('--type NAME')} (or the Type list). The surfaces are the standard level-set "
               "approximations with period 2π; one period fills one cell."))
    s += gallery([(v["t-gyroid"], "gyroid"), (v["t-schwarz-p"], "schwarz-p"), (v["t-diamond"], "diamond"),
                  (v["t-neovius"], "neovius"), (v["t-iwp"], "iwp"), (v["t-frd"], "frd"),
                  (v["t-fks"], "fischer-koch-s")], 4, height=34 * mm)
    s += table([
        ["Name", "Aliases", "f(x, y, z)"],
        ["gyroid", "g", "sin x cos y + sin y cos z + sin z cos x"],
        ["schwarz-p", "p, primitive", "cos x + cos y + cos z"],
        ["diamond", "d, schwarz-d", "sin x sin y sin z + sin x cos y cos z + cos x sin y cos z + cos x cos y sin z"],
        ["neovius", "", "3 (cos x + cos y + cos z) + 4 cos x cos y cos z"],
        ["iwp", "schoen-iwp", "2 (cos x cos y + cos y cos z + cos z cos x) − (cos 2x + cos 2y + cos 2z)"],
        ["frd", "schoen-frd", "4 cos x cos y cos z − (cos 2x cos 2y + cos 2y cos 2z + cos 2z cos 2x)"],
        ["fischer-koch-s", "fks", "cos 2x sin y cos z + cos x cos 2y sin z + sin x cos y cos 2z"],
    ], [0.17, 0.15, 0.68])
    s += note("FRD and Fischer-Koch S have the finest detail: they need finer patch grids, so they take longer "
              "and give the most faces per cell. Neovius networks at low solid fractions break into separate "
              "blobs.")
    s += H2("Form and part")
    s += bullets([
        f"{c('--form sheet')} (default) is the region |f| ≤ c: a wall around the surface, with fluid on both "
        "sides. The wall thickness follows from the level set, so it is not exactly constant (the report gives "
        "its range).",
        f"{c('--form network')} is the region f ≤ c: one of the two channel networks, like thick struts.",
        f"{c('--part solid')} (default) exports the lattice; {c('--part void')} exports the shape minus the "
        "lattice: the fluid domain for CFD.",
        "A sheet's fluid is two channels that never meet: <b>A</b> where f > c and <b>B</b> where f &lt; −c. "
        f"In STEP each goes to a file of its own: {c('-o x.step')} writes {c('x_A.step')} and {c('x_B.step')}. "
        "STL keeps both in one file.",
    ])
    s += gallery([(v["f-network"], "Gyroid network, solid (--form network --density 0.3)"),
                  (v["f-void"], "Gyroid sheet, fluid: channel A blue, channel B red (--part void)")], 2,
                 height=58 * mm)
    s += H2("Thickness")
    s += bullets([
        f"{c('--density R')} (default 0.2) finds the iso level c that gives solid fraction R in one unit cell. "
        f"It always means the <i>solid</i> fraction, also with {c('--part void')}.",
        f"{c('--level C')} sets the iso level directly instead.",
        f"{c('--wall T')} (sheet) sets the wall thickness in mm. The walls then lie about T/2 on each side of the "
        "surface, so the thickness is nearly constant: 0.48 to 0.49 mm for a gyroid at 0.5 mm, against 0.58 to "
        "0.71 mm by solid fraction 0.2. The solid fraction follows from it and is printed. Smooth STEP takes "
        "several times longer this way.",
        "On bent or trimmed shapes, cells are stretched or cut, so the real solid fraction differs from R. "
        "tpms samples the field and prints the value it gives.",
    ])
    s += H2("Grading")
    s += bullets([
        f"{c('--density R0 R1')} changes the solid fraction linearly from R0 at the start of an axis to R1 at its "
        f"end; {c('--grade x|y|z')} picks the axis (default x). {c('--level C0 C1')} grades the iso level and "
        f"{c('--wall T0 T1')} the wall thickness. On a tube, torus or sphere shell the axes are the cell axes: x "
        "radial, y around, z along the axis, around the tube, or from south to north.",
        f"{c('--cell-grade x|y|z L1')} (box) changes the cell length along that axis from {c('--cell')} at the "
        "start to L1 mm at the end; cells stretch along the axis. The box length along it is n × (L0 + L1) / 2.",
        "Both can be combined and work with every output. tpms prints the overall solid fraction.",
    ])
    s += gallery([(v["g-density"], "Solid fraction 0.1 to 0.4 along x (--density 0.1 0.4)"),
                  (v["g-cell"], "Cells from 5 to 15 mm along x (--cell 5 --cell-grade x 15)")], 2,
                 height=50 * mm)
    s += H2("Offset")
    s.append(P(f"{c('--shift S')} (or {c('Sx Sy Sz')}, in cells) moves the lattice, so the faces of the shape cut "
               "it elsewhere: away from thin slivers at the box faces, or with a node at the center of a unit "
               f"cell ({c('--shift 0.5')} for Schwarz P)."))

    # 5 Shapes
    s += H1("Shapes")
    s += table([
        ["--shape", "Size options", "Cells", "Kind"],
        [c("box") + " (default)", "—", f"{c('--cell L')}, {c('--cells nx ny nz')}", "box [0, n × L]"],
        [c("tube"), c("--radius R0 R1 --height H"), c("--cells radial around axial"), "cells follow the shape"],
        [c("torus-shell"), c("--radius R0 R1 --major R"), c("--cells radial ring tube"), "cells follow the shape"],
        [c("sphere-shell"), c("--radius R0 R1 --polar A0 A1"), c("--cells radial around S-N"), "cells follow the shape"],
        [c("cylinder"), c("--radius R --height H"), c("--cell L"), "trimmed"],
        [c("sphere"), c("--radius R"), c("--cell L"), "trimmed"],
        [c("torus"), c("--radius r --major R"), c("--cell L"), "trimmed"],
    ], [0.19, 0.3, 0.3, 0.21])
    s += code("tpms --shape tube --radius 15 25 --height 30 --cells 1 12 3 -o tube.step\n"
              "tpms --shape sphere-shell --radius 15 25 --cells 1 12 5 --polar 20 160 -o shell.step\n"
              "tpms --shape sphere --radius 15 --cell 10 -o sphere.step")
    s += gallery([(v["s-tube"], "tube"), (v["s-torus-shell"], "torus-shell"), (v["s-sphere-shell"], "sphere-shell"),
                  (v["s-cylinder"], "cylinder"), (v["s-sphere"], "sphere"), (v["s-torus"], "torus")], 3, height=40 * mm)
    s += bullets([
        "<b>Cells follow the shape</b> (tube, torus shell, sphere shell): the lattice is bent into the shape, so no cell is cut "
        "on the curved walls. The cell count around must be whole, so the lattice closes on itself. "
        f"{c('--radial log')} (default) keeps cells square: they grow with the radius; {c('--radial linear')} "
        "keeps their radial size constant. The GUI shows the resulting cell sizes.",
        "<b>Sphere shell</b>: cells that follow a sphere cannot reach its poles (one point would get a different "
        f"lattice for each longitude), so the shell runs from A0 to A1 degrees from the +z axis ({c('--polar')}, "
        "default 20 160) and stays open around the poles. Its ends are cones from the center, or the plane z = 0 "
        f"at 90 degrees: {c('--polar 10 90')} is a dome. Cells narrow around toward the ends; for a closed round "
        "part use the trimmed sphere.",
        "<b>Trimmed</b> (cylinder, sphere, torus): a straight lattice of --cell cells, cut by the shape. Use "
        "these for solid round parts: bending cells around an axis or a pole would squeeze them to nothing.",
        "Smooth STEP works for every shape, and the cut faces are exact: cylinders and end planes, tori, or a "
        "sphere. Trimmed shapes need finer patch grids, so they take from about a second (cylinder) to a minute "
        "(sphere, torus).",
        f"If a trimmed shape fails because a wall passes very close to its axis or center, move the lattice a "
        f"little: for example {c('--shift 0.13 0.21 0.07')}.",
    ])

    # 6 Special domains
    s += H1("Unit cells, plates and buffers")
    s += H2("Unit cell for homogenization")
    s.append(P(f"{c('--unit-cell')} writes a periodic cell (RVE) for FEA with periodic boundary conditions:"))
    s += bullets([
        f"One cell (or {c('--cells N')} for a supercell), in a box only, without grading.",
        "Opposite faces match exactly: point for point and edge for edge. tpms checks the mesh points on each "
        "pair of faces and fails if they differ.",
        f"Exact edges and every piece kept ({c('--exact')} and {c('--min-fragment 0')} are implied): a small piece "
        "at a corner continues in the next cell.",
    ])
    s.append(picture(v["u-cell"], 62 * mm, "Schwarz P network cell with a node at its center "
                                           "(--unit-cell --form network --shift 0.5)"))
    s += H2("Plates and buffers")
    s += bullets([
        f"{c('--plates x|y|z T [T1]')} (solid): plates T mm thick on the faces at the start and end of an axis "
        "(T1 at the end; 0 for none), fused with the lattice: sandwich panels, walls.",
        f"{c('--buffers x|y|z L [L1]')} (void): empty zones L mm long before and after the lattice, fused with "
        "the fluid, so the flow can develop: inlet and outlet regions for CFD. A sheet's two channels meet in "
        "a buffer, so the fluid is one solid.",
        "The slabs lie outside the lattice box, which keeps its size. Solid fractions printed by tpms and the "
        "report are of the lattice box.",
    ])
    s += gallery([(v["p-plates"], "--plates z 1.5"), (v["p-buffers"], "--part void --buffers x 10")], 2,
                 height=52 * mm)
    s += H2("Channel buffers for heat exchangers")
    s.append(P(f"{c('--channel-buffers A B L [T]')} gives each fluid channel of a sheet its own inlet and outlet: "
               "a two-fluid heat exchanger in cross flow."))
    s += bullets([
        "Channel A flows along axis A, channel B along axis B (x, y or z, different).",
        "On the two faces across its axis, each channel gets a header plate T mm thick (default 1) and, beyond "
        "it, a buffer L mm long. The plate has holes exactly where its channel opens on the face, and closes "
        "the other channel there.",
        f"{c('--part void')} writes channel A with its holes and buffers to {c('NAME_A.step')} and channel B to "
        f"{c('NAME_B.step')}. {c('--part solid')} writes the lattice fused with the four plates.",
        "The three files fill the box, plates and buffers without gaps or overlaps: the solid's and both "
        "fluids' volumes add up to the whole.",
        "Box with a sheet only; not with --unit-cell, --cell-grade, --plates or --buffers.",
    ])
    s += gallery([(v["hx-void"], "Fluids: channel A (blue) along x, B (red) along y"),
                  (v["hx-solid"], "Solid: the lattice with its four header plates")], 2, height=58 * mm)
    s += code("tpms --part void --channel-buffers x y 10 1 -o hx_fluid.step   # hx_fluid_A.step and _B.step\n"
              "tpms --channel-buffers x y 10 1 -o hx_solid.step")
    s += note("Ports and walls are named: channel A's buffer ends are its inlet and outlet (x- and x+ here), the "
              "plates and hole walls are <i>wall A</i> or <i>wall B</i> after the channel they face. The two "
              "faces across the third axis (z- and z+ here) stay open cut faces: set them as walls, symmetry or "
              "periodic.")

    # 7 Output files
    s += H1("Output files")
    s.append(P(f"Give one or more outputs with {c('-o FILE')}; the extension picks the format."))
    s += table([
        ["Output", "What you get", "Use it for"],
        [c(".step") + " (default)", "smooth B-spline faces within --tol of the exact surface; patches merged "
                                    "into large trimmed faces (3–8 times fewer faces)", "CAD, meshing for FEA/CFD"],
        [c(".step --exact"), "one face per patch; neighbor faces share their edges exactly", "strict CAD checks, "
                                                                                              "unit cells"],
        [c(".step --faceted"), "one planar face per triangle of the mesh", "tools that need a faceted B-rep"],
        [c(".stl"), "binary triangle mesh (all solids in one file)", "3D printing, mesh-based tools"],
        [c(".json"), "the property report, structured", "scripts"],
        [c(".csv"), "the property report, one header line and one line of values", "spreadsheets, sweeps"],
    ], [0.2, 0.55, 0.25])
    s += H2("Smooth STEP options")
    s += bullets([
        f"{c('--tol T')}: the largest distance in mm between a face and the exact surface. Default: 0.001 × the "
        "smallest cell edge (0.01 mm for 10 mm cells).",
        f"{c('--patches N')}: the patch grid per cell edge to start from (default 3). If a check fails, tpms tries "
        "finer grids by itself and prints the one it used; an ungraded lattice is tried on a single cell "
        "first, which is fast.",
        f"{c('--exact')}: one face per patch with exactly shared edges. Without it, merged faces have edges within "
        "--tol of the faces (tolerant edges, as in most industrial STEP files; the file says so).",
    ])
    s += table([
        ["Gyroid sheet, 2 × 2 × 2 cells of 10 mm", "Faces", "File", "Read time (OpenCASCADE)"],
        ["smooth STEP, merged (default)", "371", "5 MB", "1.6 s"],
        ["smooth STEP, --exact", "1 506", "5 MB", "2.0 s"],
        ["faceted STEP, --res 24", "261 396", "196 MB", "117 s"],
    ], [0.46, 0.14, 0.14, 0.26])
    s += H2("Names in the STEP file")
    s += table([
        ["Name", "What"],
        [f"{c('solid')}, {c('fluid')}, {c('fluid A')}, {c('fluid B')}", "solids: the part; a sheet's fluid "
         "channels A (f > c) and B (f &lt; −c). More pieces get 2, 3, …"],
        [f"{c('wall')}, {c('wall A')}, {c('wall B')}", "walls between solid and fluid; on a sheet, A and B tell "
         "which channel they face (also header plates and their holes)"],
        [f"{c('x-')} {c('x+')} {c('y-')} {c('y+')} {c('z-')} {c('z+')}", "box faces, and the outer faces of "
         "plates and buffers: the inlet of --buffers x is x-"],
        [f"{c('inner')} {c('outer')} {c('bottom')} {c('top')}", "tube, torus shell and sphere shell faces"],
        [f"{c('side')} {c('bottom')} {c('top')} {c('surface')}", "cylinder side and ends; sphere and torus surface"],
    ], [0.36, 0.64])
    s.append(P("Use the names to pick boundaries in the CAD, CFD or FEA setup. Readers differ in what they "
               "import: some show face names only with an option, and a file with one solid often shows the "
               "file name for it."))
    s += H2("Meshes: STL and faceted STEP")
    s.append(P(f"{c('--res N')} sets the grid steps per cell edge (default 24). The triangle count grows with N²; "
               "thin sheets need a higher N. tpms warns when the solid fraction of the mesh is more than 5 % away "
               "from the field's, which usually means walls thinner than the grid."))
    s += H2("Size and speed")
    s += bullets([
        "Face count and file size grow with the number of cells. Start with a few cells, check the import, "
        "then scale up.",
        "FRD and Fischer-Koch S give the most faces per cell; Schwarz P and gyroid the fewest.",
        "A tighter --tol gives more faces. The default suits meshing for simulation.",
    ])
    s += H2("Fragments and cavities")
    s.append(P(f"Pieces smaller than {c('--min-fragment F')} × the largest piece (default 0.01) are dropped: tiny "
               "slivers cut off at the faces of the shape would break a simulation setup. Use 0 to keep all. "
               "Closed cavities inside a solid (for example isolated network blobs in the fluid) are written as "
               "voids of their solid; tiny ones are filled."))

    # 8 Report
    s += H1("Property report")
    s.append(P(f"{c('--report')} prints the properties of the solid and the fluid, whatever {c('--part')} exports. "
               f"An output ending in {c('.json')} or {c('.csv')} writes the report as a file; the GUI shows it in "
               "a window."))
    s.append(picture(report_window("report"), WIDTH * 0.8, "The property report window in tpms-gui"))
    s += table([
        ["Value", "Meaning"],
        ["Solid / fluid volume", "mm³, with solid fraction and porosity"],
        ["Wetted area", "area of the walls between solid and fluid, mm² (without the cut faces)"],
        ["Specific surface", "wetted area / shape volume, 1/mm and m²/m³"],
        ["Hydraulic diameter", "D<sub>h</sub> = 4 × fluid volume / wetted area, overall and per fluid channel"],
        ["Wall thickness (sheet)", "mean 2 × solid volume / wetted area, and the local range on the walls"],
        ["Solid size (network)", "4 × solid volume / wetted area, like a strut diameter"],
        ["Mass", "with --material: solid volume × density, g"],
        ["Stiffness E/E<sub>s</sub>", "Gibson–Ashby bounds at solid fraction φ: φ² (bending-dominated) to φ/3 "
                                      "(stretching-dominated); sheets usually lie near the upper value, networks "
                                      "near the lower. In GPa with --material ... E"],
        ["Permeability", "per fluid channel, Kozeny–Carman: K = ε D<sub>h</sub>² / 80 (ε the channel's porosity)"],
        ["Pressure drop", "with --flow, per channel: Reynolds number, and the pressure gradient μu/K + 7/6 ρu² / "
                          "(ε² D<sub>h</sub>) (Darcy plus Ergun's inertial term), in Pa/m and over the box along "
                          "the flow axis"],
        ["Pieces", "each solid and fluid piece, its volume and the faces of the shape it touches"],
        ["Connects", "for each pair of opposite faces, whether the solid and each fluid channel connect them: flow "
                     "through for CFD, a load path for FEA"],
    ], [0.24, 0.76])
    s += H2("Inputs")
    s += bullets([
        f"{c('--material RHO [E]')}: the solid's density in kg/m³, and its Young's modulus in GPa.",
        f"{c('--flow U [MU RHO]')}: the superficial (inlet) velocity in m/s, and the fluid's viscosity in Pa s and "
        "density in kg/m³ (default water: 0.001 and 998). The flow axis is the channel's --channel-buffers or "
        "--buffers axis, else x.",
    ])
    s += note("The values come from closed meshes at --res steps per cell: volumes and areas within about 1 % at "
              "the default 24 (checked against CAD mass properties). Mass is exact for that volume. Stiffness, "
              "permeability and pressure drop are rough estimates from correlations, for comparing designs: use "
              "FEA on a unit cell and CFD for real values.")
    s += H2("Parameter sweeps")
    s.append(P("The CSV columns are the same for every run (empty where a value does not apply), so the lines of "
               "several runs can be put together into one table. PowerShell:"))
    s += code("""foreach ($d in 0.1, 0.2, 0.3) {
    tpms --density $d --flow 0.1 -o "sweep_$d.csv"
}
$rows = Get-ChildItem sweep_*.csv | ForEach-Object { Get-Content $_ }
$rows[0], ($rows | Where-Object { $_ -notlike 'surface,*' }) | Set-Content sweep.csv""")

    # 9 Workflows
    s += H1("Workflows")
    s += H2("FEA of a lattice part")
    s += bullets([
        "Export the solid as smooth STEP (default). For a sandwich panel add plates: "
        f"{c('--plates z 1.5')}.",
        "Import into the CAD or FEA tool and mesh the solid. Use the face names: load and support on the box "
        f"faces ({c('z-')}, {c('z+')}), and {c('wall')} faces for pressure.",
        "Check with the report that the solid connects the loaded faces (<i>Connects</i>).",
    ])
    s += code("tpms --type gyroid --density 0.25 --cell 8 --cells 5 5 2 --plates z 1.5 -o panel.step")
    s += H2("CFD through a lattice")
    s += bullets([
        f"Export the fluid with {c('--part void')} and inlet and outlet buffers along the flow: {c('--buffers x 10')}.",
        f"Inlet {c('x-')}, outlet {c('x+')}, walls {c('wall')}; the other box faces as walls or symmetry.",
        "The report gives D<sub>h</sub>, porosity, permeability and, with --flow, a first pressure drop estimate "
        "to compare with.",
    ])
    s += code("tpms --part void --buffers x 10 --cells 4 2 2 -o fluid.step --flow 0.2 --report")
    s += H2("Two-fluid heat exchanger (conjugate heat transfer)")
    s += bullets([
        "Write the solid and both fluids with the same lattice options and channel buffers (three files).",
        "Import all three; they touch exactly, so the CFD tool can couple their shared faces as interfaces.",
        "Inlets and outlets are the buffer ends: x- and x+ for channel A along x, y- and y+ for B along y.",
    ])
    s += code("""tpms --cells 4 4 2 --part void --channel-buffers x y 10 1 -o hx_fluid.step
tpms --cells 4 4 2 --channel-buffers x y 10 1 -o hx_solid.step""")
    s += H2("Effective properties from a unit cell")
    s += bullets([
        f"Export one cell with {c('--unit-cell')}; its opposite faces match, so periodic boundary conditions "
        "apply directly.",
        f"Check it with {c('python tools/check_step.py --periodic cell.step')} (section 11.2).",
    ])
    s += code("tpms --unit-cell --type diamond --density 0.3 --cell 8 -o cell.step")
    s += H2("3D printing")
    s.append(P(f"Write STL. Raise {c('--res')} for thin walls or large cells; the warning about the solid fraction "
               "tells you when the grid is too coarse."))
    s += code("tpms --density 0.15 --cell 6 --cells 6 --res 32 -o print.stl")
    s += H2("Importing into CAD")
    s += bullets([
        "STEP files are AP214, in mm, with one solid per piece (and names). Import them as solids, not as "
        "surfaces.",
        "Merged files have tolerant edges (within --tol); the STEP file states that tolerance, and CAD systems "
        "use it when they sew. If a tool is strict about edges, export with --exact.",
        "Before a large model, try 2 × 2 × 2 cells, check volume and validity in the CAD tool, then scale up.",
        "No CAD at hand? A free Onshape account (Parasolid kernel, the same as Siemens NX) imports the files "
        "and gives mass properties to compare with the report.",
    ])

    # 10 Reference
    s += H1("Command-line reference")
    s.append(P(f"Usage: {c('tpms [options] -o FILE [-o FILE ...]')}. {c('tpms-gui')} takes the same options and "
               "opens with them."))

    def ref(title, rows):
        return [KeepTogether(H2(title) + table([["Option", "Meaning"]] + [[c(a), b] for a, b in rows], [0.31, 0.69]))]

    s += ref("Lattice", [
        ("--type NAME", "surface type (default gyroid); see --list"),
        ("--form sheet|network", "sheet: |f| ≤ c, a wall around the surface; network: f ≤ c (default sheet)"),
        ("--density R [R1]", "target solid fraction of one cell, 0 &lt; R &lt; 1 (default 0.2); with R1 graded "
                             "from R to R1 along --grade"),
        ("--level C [C1]", "iso level c (or graded from C to C1), instead of --density"),
        ("--wall T [T1]", "sheet: wall thickness in mm (or graded from T to T1), instead of --density"),
        ("--grade x|y|z", "axis of a graded density or level (default x)"),
        ("--part solid|void", "the solid, or the shape minus the solid (fluid domain)"),
    ])
    s += ref("Shape", [
        ("--shape NAME", "box, tube, torus-shell, sphere-shell (cells follow the shape); cylinder, sphere, torus "
                         "(trimmed)"),
        ("--cell L [Ly Lz]", "unit cell size in mm (box and trimmed shapes; default 10)"),
        ("--cells N [Ny Nz]", "cell counts (default 2): box x y z; tube radial, around, along; torus shell "
                              "radial, around the ring, around the tube; sphere shell radial, around, south to "
                              "north"),
        ("--radius R | R0 R1", "tube, torus or sphere shell: inner and outer radius; cylinder, sphere: radius; "
                               "torus: tube radius"),
        ("--polar A0 A1", "sphere shell: where it starts and ends, degrees from +z (default 20 160)"),
        ("--height H", "tube, cylinder: length along z"),
        ("--major R", "torus shell, torus: from the axis to the tube center"),
        ("--radial log|linear", "tube, torus or sphere shell: radial spacing (default log: square cells)"),
        ("--cell-grade x|y|z L1", "box: cell length from --cell to L1 mm along the axis"),
        ("--shift S [Sy Sz]", "offset of the lattice in cells (default 0)"),
        ("--unit-cell", "box: a periodic cell for homogenization (--cells defaults to 1)"),
        ("--plates x|y|z T [T1]", "box, solid: plates T mm thick at the start and end of the axis"),
        ("--buffers x|y|z L [L1]", "box, void: buffer zones L mm long before and after the lattice"),
        ("--channel-buffers A B L [T]", "box, sheet: channel A along axis A, B along B, each with buffers L mm "
                                        "long behind header plates T mm thick (default 1)"),
    ])
    s += ref("Smooth STEP", [
        ("--patches N", "patch grid per cell edge to start from (default 3); stepped up as needed"),
        ("--tol T", "largest fit error in mm (default 0.001 × the smallest cell edge)"),
        ("--exact", "one face per patch, exactly shared edges (default: merged faces)"),
        ("--faceted", "one planar face per triangle instead"),
    ])
    s += ref("Meshes, output and report", [
        ("--res N", "STL and faceted STEP: grid steps per cell edge (default 24)"),
        ("--min-fragment F", "drop pieces smaller than F × the largest (default 0.01; 0 keeps all)"),
        ("-o, --out FILE", ".step / .stp, .stl, .json or .csv; repeat for several files"),
        ("--report", "print the property report"),
        ("--material RHO [E]", "report: solid density kg/m³, Young's modulus GPa"),
        ("--flow U [MU RHO]", "report: velocity m/s, fluid viscosity Pa s and density kg/m³ (default water)"),
        ("--list", "list the surface types"),
        ("-h, --help", "show the help"),
    ])

    # 11 Troubleshooting
    s += H1("Checks and troubleshooting")
    s.append(P("Every smooth STEP is checked before it is written: each edge is shared by exactly two faces in "
               "opposite directions (closed solids), every corner lies on the surface, the faces are within "
               "--tol of the exact surface, and the enclosed volume matches the field. If a check fails, tpms "
               "tries finer patch grids by itself; it reports an error only when the finest one fails too. Mesh "
               "outputs are checked for closed solids."))
    s += table([
        ["Message", "What to do"],
        ["warning: solid fraction ... is far from the field's", "walls are thinner than the grid: raise --res (STL, "
         "faceted) or --patches (smooth STEP), or use larger cells or a higher solid fraction"],
        ["note: dropped N small fragment(s)", "slivers cut off at the shape's faces were removed. Use --shift to move "
         "the cut, or --min-fragment 0 to keep them"],
        ["error: ... at N patches per cell. Try a higher --patches", "raise --patches; for trimmed shapes also try "
         "a small --shift (walls through the axis or center are hard to fit)"],
        ["warning: mesh not closed", "try another --res"],
        ["note: a fluid solid lies in both channels", "the two fluid channels meet (buffers join them, or a sheet "
         "is thinner than the grid), so they stay in one file"],
        ["error: (usage)", "the message names the option; tpms refuses conflicting options instead of guessing"],
    ], [0.42, 0.58])
    s += H2("Known limits")
    s += bullets([
        "The surfaces are level-set approximations of the minimal surfaces, and sheet walls are not exactly "
        "constant in thickness.",
        "Merged Neovius sheets on tubes and torus shells can be read as invalid; use --exact for those.",
        "A sphere shell stays open around its poles; --wall makes smooth STEP several times slower.",
        "Large lattices give large files: face count grows with the number of cells.",
        "Plates and buffers work on one axis of a box; channel buffers need a box and a sheet; grading is "
        "linear along one axis.",
        "The report's stiffness, permeability and pressure drop are correlations, not simulations.",
    ])

    # 12 OpenCASCADE checker
    s += H2("Checking files with OpenCASCADE")
    s.append(P(f"{c('tools/check_step.py')} reads STEP files with OpenCASCADE, an independent CAD kernel, and "
               "reports solids, closed shells, faces, validity and volume. It exits with 1 if a file is "
               "unreadable, invalid or open."))
    s += code("""python -m venv .venv
.venv/Scripts/python -m pip install cadquery-ocp      (.venv/bin/python on Linux and macOS)
.venv/Scripts/python tools/check_step.py lattice.step
.venv/Scripts/python tools/check_step.py --periodic cell.step     # opposite faces of a unit cell
.venv/Scripts/python tools/check_step.py --names hx_fluid_A.step  # solid and face names""")
    s.append(P(f"With this environment in the project folder, {c('ctest')} also reads every test file back with "
               "OpenCASCADE."))
    s.append(Spacer(1, 10))
    s.append(P('<font color="#5A6B7C">tpms is open source under the MIT license.</font>'))

    doc = Guide(OUT)
    doc.multiBuild(s)
    print("wrote", OUT)


if __name__ == "__main__":
    if "--no-shots" not in sys.argv:
        take_shots()
    build()
