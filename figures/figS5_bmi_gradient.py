"""Supplementary Fig. S5 | Strength of association and adiposity gradient of each disease (exploratory
meta-regression): 30-min log RR (full cohort) against log HR per SD of BMI, a sleep -> MVPA, b MVPA -> sleep.

Data: figures/data/fig4_figS5_data.R -> figS5_points.csv, figS5_fits.csv
"""
from style import *  # noqa: F403

# italic P via mathtext must also be Arial (local override; style.py untouched)
plt.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic"})

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
from matplotlib.transforms import blended_transform_factory

pts = pd.read_csv(DATADIR / "figS5_points.csv")
fits = pd.read_csv(DATADIR / "figS5_fits.csv")

FW, FH = 112.0, 56.0  # mm
fig = plt.figure(figsize=(FW * MM, FH * MM))


def ax_mm(left, top, w, h):
    return fig.add_axes([left / FW, 1 - (top + h) / FH, w / FW, h / FH])


MINUS = "−"


def fmt_num(x, nd=2):
    s = f"{x:.{nd}f}"
    return s.replace("-", MINUS)


def ratio_fmt(v, _pos=None):
    s = f"{v:.2f}".rstrip("0").rstrip(".") if v not in (1.0,) else "1"
    return s.replace("-", MINUS)


def log_ratio_axis(axis, ticks):
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_minor_locator(NullLocator())
    axis.set_major_formatter(FuncFormatter(lambda v, p: ratio_fmt(v)))


# ----------------------------------------------------------------------------------------------
# a: bubble meta-regression
# ----------------------------------------------------------------------------------------------
TEXT = "#1A1A1A"
w_all = 1 / pts["se_log"] ** 2
S_MAX = 170.0  # marker area (pt^2) for the most precise estimate; area proportional to 1/SE^2
s_scale = S_MAX / w_all.max()

a_specs = {
    "sleep_to_mvpa": dict(left=12.5, title="Sleep→MVPA, 30 min/day", ylim=(0.765, 1.085),
                          yticks=[0.8, 0.85, 0.9, 0.95, 1.0]),
    "mvpa_to_sleep": dict(left=59.5, title="MVPA→sleep, 30 min/day", ylim=(0.925, 1.45),
                          yticks=[1.0, 1.1, 1.2, 1.3, 1.4]),
}
# label placement: (dx, dy) in points from bubble centre, ha, va
a_labels = {
    "sleep_to_mvpa": {
        "gallstone": (-6, -4, "right", "center"),
        "peptic_ulcer": (-3, -5, "right", "top"),
        "gerd": (-5, -7, "right", "top"),
        "diverticular": (4, 6, "left", "bottom"),
        "ibs": (0, -4, "center", "top"),
    },
    "mvpa_to_sleep": {
        "gallstone": (3, 5, "right", "bottom"),
        "peptic_ulcer": (-4, 0, "right", "center"),
        "gerd": (-5, 3, "right", "bottom"),
        "diverticular": (5, -4, "left", "top"),
        "ibs": (0, 4, "center", "bottom"),
    },
}
A_TOP, A_H, A_W = 8.0, 38.0, 35.0
axes_a = []
for direction, sp in a_specs.items():
    ax = ax_mm(sp["left"], A_TOP, A_W, A_H)
    axes_a.append(ax)
    d = pts[pts["direction"] == direction].set_index("disease").loc[DISEASES]
    x = np.exp(d["log_hr_bmi"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(1, color="#8C8C8C", lw=0.6, ls=(0, (3, 2)), zorder=0)
    f = fits[fits["direction"] == direction].iloc[0]
    xx = np.linspace(1.095, 1.56, 50)
    ax.plot(xx, np.exp(f["intercept"] + f["slope"] * np.log(xx)), color="#4D4D4D", lw=0.9, zorder=1)
    order = np.argsort(-(1 / d["se_log"] ** 2).values)  # draw large bubbles first
    for i in order:
        o = d.index[i]
        ax.scatter(x.iloc[i], d["rr"].iloc[i], s=s_scale / d["se_log"].iloc[i] ** 2,
                   color=DISEASE_COLOR[o], alpha=0.8, edgecolor="white", linewidth=0.5, zorder=3)
    for o in DISEASES:
        dx, dy, ha, va = a_labels[direction][o]
        ax.annotate(DISEASE_LABEL[o], (x[o], d.loc[o, "rr"]), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va=va, fontsize=5.5, color=TEXT, zorder=4)
    ax.set_xlim(1.07, 1.62)
    ax.set_ylim(*sp["ylim"])
    log_ratio_axis(ax.xaxis, [1.1, 1.2, 1.3, 1.4, 1.5])
    log_ratio_axis(ax.yaxis, sp["yticks"])
    ax.set_xlabel("HR per SD of BMI")
    ax.set_ylabel("RR", labelpad=2)
    ax.set_title(sp["title"], pad=4)
    ptxt = f"{f['p_dispersion']:.3f}"
    note = f"Slope {fmt_num(f['slope'], 3)}\nDispersion-adjusted $P$ = {ptxt}"
    if direction == "sleep_to_mvpa":
        ax.text(0.98, 0.985, note, transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=TEXT)
    else:
        ax.text(0.98, 0.015, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5, color=TEXT)



# panel letters and a key to the circle size
for ax, lab in zip(axes_a, "ab"):
    ax.text(-0.30, 1.05, lab, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")
fig.text(0.5, 1 - (A_TOP + A_H + 9.5) / FH, "Circle area is proportional to the precision (1/SE$^2$) of each "
         "disease's log RR; lines are the weighted regressions.", ha="center", va="top", fontsize=5.5, color=GREY)
save(fig, "FigS5_bmi_gradient")
print("saved")
