"""Shared style for the figures (matplotlib; widths follow common journal column sizes).

Import in every figure script:  from style import *
Widths: single column 89 mm, 1.5 column 120 mm, double column 183 mm; max height 247 mm.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.frameon": False,
    "lines.linewidth": 1.2,
})

PROJ = Path(__file__).resolve().parents[1]  # repository root
FIGDIR = PROJ / "figures" / "output"
DATADIR = PROJ / "figures" / "data"
MM = 1 / 25.4
W1, W15, W2 = 89 * MM, 120 * MM, 183 * MM

# disease order follows the evidence tiers used in the manuscript
DISEASES = ["gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs"]
DISEASE_LABEL = {"gallstone": "Gallstone disease", "peptic_ulcer": "Peptic ulcer", "gerd": "GERD",
                 "diverticular": "Diverticular disease", "ibs": "IBS",
                 "gallstone_surgery": "Gallstone disease requiring cholecystectomy",
                 "nco_cyst": "Follicular cyst (negative control)", "nco_earwax": "Ear wax (negative control)"}
# Okabe-Ito-derived, colour-blind safe; negative controls in grey
DISEASE_COLOR = {"gallstone": "#D55E00", "peptic_ulcer": "#CC79A7", "gerd": "#0072B2",
                 "diverticular": "#009E73", "ibs": "#E69F00",
                 "nco_cyst": "#7F7F7F", "nco_earwax": "#B0B0B0"}
BEHAVIOUR_COLOR = {"sleep": "#B4C0E4", "inactivity": "#CFCECE", "light": "#42949E", "mvpa": "#0F4D92"}
BEHAVIOUR_LABEL = {"sleep": "Sleep", "inactivity": "Inactivity", "light": "Light activity", "mvpa": "MVPA"}
GROUP_LABEL = {"short": "<7 h", "normal": "7–9 h", "long": ">9 h"}
GROUP_MARKER = {"short": "^", "normal": "o", "long": "s"}
GREY = "#5A5A5A"


def panel_label(ax, letter, x=-0.12, y=1.04):
    """Bold lower-case panel letter (Nature style)."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    """Write svg, pdf, 600-dpi tiff, 300-dpi png and a small preview png to figures/output/."""
    FIGDIR.mkdir(parents=True, exist_ok=True)
    for ext, kw in [("svg", {}), ("pdf", {}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}}),
                    ("png", {"dpi": 300})]:
        fig.savefig(FIGDIR / f"{name}.{ext}", bbox_inches="tight", **kw)
    fig.savefig(FIGDIR / f"{name}_preview.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
