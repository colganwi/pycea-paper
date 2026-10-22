"""Shared paths, plotting style, and figure-assembly helpers."""

from pathlib import Path

import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import colors  # noqa: E402  (project-wide default colors)

DATA = ROOT / "data"
PLOTS = ROOT / "plots"
FIGURES = ROOT / "manuscript" / "figures"

# Cluster colors shared by all panels (first colors of the project palette)
CLUSTER_PALETTE = colors


def set_style():
    plt.rcParams.update(
        {
            "svg.fonttype": "none",  # keep text as text in SVGs
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7,
            "axes.titlesize": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.linewidth": 0.5,
            "pdf.fonttype": 42,
        }
    )


def save_panel(fig, name, subdir):
    """Save a panel as SVG in plots/<subdir>/<name>.svg."""
    out = PLOTS / subdir
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.svg", bbox_inches="tight", transparent=True)
    plt.close(fig)
    return out / f"{name}.svg"


def save_panel_full(fig, name, subdir, pad_inches=0.02, width=None):
    """Like save_panel, but the crop box also covers legends and inset (colorbar) axes added by pycea.

    ``bbox_inches="tight"`` can clip colorbar legends, which are inset axes drawn outside the main axes.
    """
    from matplotlib.legend import Legend
    from matplotlib.transforms import Bbox

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = []
    for ax in fig.axes:
        boxes.append(ax.get_tightbbox(renderer))
        stack = list(ax.child_axes) + [a for a in ax.artists if isinstance(a, Legend)]
        if ax.get_legend() is not None:
            stack.append(ax.get_legend())
        while stack:
            artist = stack.pop()
            boxes.append(artist.get_window_extent(renderer))
            stack.extend(getattr(artist, "child_axes", []))
    for artist in list(fig.texts) + list(fig.patches) + list(fig.artists):  # figure-level text and shapes
        boxes.append(artist.get_window_extent(renderer))
    box = Bbox.union([b for b in boxes if b is not None]).transformed(fig.dpi_scale_trans.inverted())
    box = box.padded(pad_inches)
    if width is not None:  # fixed crop width (inches), centered on the first axes, so panels share one scale
        center = fig.axes[0].get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted()).x0
        center += fig.axes[0].get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted()).width / 2
        box = Bbox.from_extents(center - width / 2, box.y0, center + width / 2, box.y1)
    out = PLOTS / subdir
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.svg", bbox_inches=box, transparent=True)
    plt.close(fig)
    return out / f"{name}.svg"


MM = 1 / 25.4  # inches per mm


def mm_figure(width, height):
    """Figure of `width` x `height` mm, for panels drawn at their final printed size."""
    return plt.figure(figsize=(width * MM, height * MM))


def mm_axes(fig, x, y, width, height, **kwargs):
    """Axes of `width` x `height` mm whose top-left corner is (x, y) mm from the top-left of the figure."""
    fw, fh = fig.get_size_inches() / MM
    return fig.add_axes([x / fw, 1 - (y + height) / fh, width / fw, height / fh], **kwargs)


def mm_text(fig, x, y, s, **kwargs):
    """Figure text at (x, y) mm from the top-left of the figure."""
    fw, fh = fig.get_size_inches() / MM
    return fig.text(x / fw, 1 - y / fh, s, **kwargs)


def save_panel_fixed(fig, name, subdir):
    """Save the whole figure without cropping, so the SVG has exactly the figure's size."""
    out = PLOTS / subdir
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.svg", transparent=True)
    plt.close(fig)
    return out / f"{name}.svg"


def visible_extent(svg_path, dpi_scale=2):
    """Horizontal extent (left, right) of the drawn content of an SVG as fractions of its width."""
    import io

    import cairosvg
    import numpy as np
    from PIL import Image

    png = cairosvg.svg2png(url=str(svg_path), scale=dpi_scale)
    alpha = np.asarray(Image.open(io.BytesIO(png)).convert("RGBA"))[:, :, 3]
    cols = np.where(alpha.max(axis=0) > 0)[0]
    return cols.min() / alpha.shape[1], (cols.max() + 1) / alpha.shape[1]


def visible_top(svg_path, dpi_scale=2):
    """Vertical position of the top of the drawn content of an SVG, as a fraction of its height."""
    import io

    import cairosvg
    import numpy as np
    from PIL import Image

    png = cairosvg.svg2png(url=str(svg_path), scale=dpi_scale)
    alpha = np.asarray(Image.open(io.BytesIO(png)).convert("RGBA"))[:, :, 3]
    rows = np.where(alpha.max(axis=1) > 0)[0]
    return rows.min() / alpha.shape[0]


class SvgFigure:
    """Compose SVG panels into a single figure (units: mm) and export to SVG/PDF.

    Panels are placed by their top-left corner and scaled to a given width or
    height (aspect ratio preserved). Everything is plain SVG so figures are
    fully reproducible from the panel files.
    """

    def __init__(self, width_mm=183, height_mm=120):
        self.w, self.h = width_mm, height_mm
        self.parts = []

    def image(self, path, x, y, width=None, height=None):
        """Place a raster image (PNG) with top-left at (x, y) mm; returns (width, height) in mm."""
        import base64
        import struct

        data = Path(path).read_bytes()
        pw, ph = struct.unpack(">II", data[16:24])  # PNG IHDR
        if width is None:
            width = height * pw / ph
        height = width * ph / pw
        b64 = base64.b64encode(data).decode()
        self.parts.append(
            f'<image x="{x}mm" y="{y}mm" width="{width}mm" height="{height}mm" '
            f'xlink:href="data:image/png;base64,{b64}" preserveAspectRatio="xMidYMid meet"/>'
        )
        return width, height

    def line(self, x1, y1, x2, y2, color="#555555", width=0.4):
        self.parts.append(
            f'<line x1="{x1}mm" y1="{y1}mm" x2="{x2}mm" y2="{y2}mm" stroke="{color}" stroke-width="{width}mm"/>'
        )

    def arrow(self, x1, y, x2, color="#555555", width=0.8, head=2.5):
        """Horizontal arrow from x1 to x2 at height y (mm) with a filled head."""
        self.line(x1, y, x2 - head * 0.8, y, color=color, width=width)
        self.parts.append(
            f'<polygon points="{x2},{y} {x2 - head},{y - head / 2} {x2 - head},{y + head / 2}" fill="{color}"/>'
        )

    @staticmethod
    def _size(path):
        import re

        head = Path(path).read_text()[:2000]
        vb = re.search(r'viewBox="([\d.\-\s]+)"', head)
        if vb:
            _, _, w, h = map(float, vb.group(1).split())
            return w, h
        w = float(re.search(r'width="([\d.]+)', head).group(1))
        h = float(re.search(r'height="([\d.]+)', head).group(1))
        return w, h

    def place(self, path, x, y, width=None, height=None):
        """Place an SVG with top-left at (x, y) mm; returns (width, height) in mm."""
        import re

        pw, ph = self._size(path)
        if width is None:
            width = height * pw / ph
        height = width * ph / pw
        body = Path(path).read_text()
        body = re.sub(r"<\?xml[^>]*\?>", "", body)
        body = re.sub(r"<!DOCTYPE[^>]*>", "", body)
        # Re-root the SVG: nested <svg> with explicit position, size and viewBox
        body = re.sub(
            r"<svg\b[^>]*>",
            f'<svg x="{x}mm" y="{y}mm" width="{width}mm" height="{height}mm" viewBox="0 0 {pw} {ph}" '
            f'preserveAspectRatio="xMidYMid meet" xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink">',
            body,
            count=1,
        )
        self.parts.append(body)
        return width, height

    def text(self, s, x, y, size=7, weight="normal", anchor="start", color="black", mono=False):
        family = "DejaVu Sans Mono, Courier New, monospace" if mono else "Arial, Helvetica, DejaVu Sans, sans-serif"
        self.parts.append(
            f'<text x="{x}mm" y="{y}mm" font-family="{family}" '
            f'font-size="{size * 0.3528}mm" font-weight="{weight}" text-anchor="{anchor}" fill="{color}">{s}</text>'
        )

    def letter(self, s, x, y):
        """Bold panel letter with baseline at (x, y) mm."""
        self.text(s, x, y, size=10, weight="bold")

    def save(self, path_stem):
        path_stem = Path(path_stem)
        path_stem.parent.mkdir(parents=True, exist_ok=True)
        svg = (
            f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{self.w}mm" height="{self.h}mm" viewBox="0 0 {self.w} {self.h}">'
            '<rect width="100%" height="100%" fill="white"/>'
            # The viewBox is in mm, so parts must use plain user units (1 unit = 1 mm): strip their "mm" suffixes.
            + "".join(self.parts).replace('mm"', '"')
            + "</svg>"
        )
        svg_path = path_stem.with_suffix(".svg")
        svg_path.write_text(svg)
        import cairosvg

        cairosvg.svg2pdf(url=str(svg_path), write_to=str(path_stem.with_suffix(".pdf")))
        return svg_path
