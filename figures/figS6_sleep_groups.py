"""Supplementary Fig. S6 | Reallocating 30 min/day from sleep to MVPA, by accelerometer-measured sleep-duration
group: full-cohort model vs separate models per sleep-duration group.

Data: figures/data/figS6_groups.csv, figS6_ratios.csv, figS6_heterogeneity.csv
      (written by figures/data/figS6_data.R from point_estimate.rds + boots/*.rds).
Run:  python3 figures/figS6_sleep_groups.py
"""
from style import *  # noqa: F403

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

# italic P via mathtext must also be Arial (local override; style.py untouched)
plt.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic"})

grp = pd.read_csv(DATADIR / "figS6_groups.csv")
rat = pd.read_csv(DATADIR / "figS6_ratios.csv")
het = pd.read_csv(DATADIR / "figS6_heterogeneity.csv").iloc[0]

NCOS = ["nco_cyst", "nco_earwax"]
ROWS = DISEASES + NCOS
ROW_LABEL = {"gallstone": "Gallstone\ndisease", "peptic_ulcer": "Peptic ulcer", "gerd": "GERD",
             "diverticular": "Diverticular\ndisease", "ibs": "IBS",
             "nco_cyst": "Follicular cyst\n(L72)", "nco_earwax": "External ear\ndisorders (H61)"}
NCO_GAP = 1.0
YPOS = {d: i for i, d in enumerate(DISEASES)}
YPOS.update({d: len(DISEASES) + NCO_GAP + i for i, d in enumerate(NCOS)})
YLIM = (YPOS["nco_earwax"] + 0.55, -0.55)          # inverted: first disease on top
GROUP_OFF = {"short": -0.27, "normal": 0.0, "long": 0.27}
MODEL_OFF = {"full_cohort": -0.15, "separate": 0.15}

MS, LW, CAP = 3.3, 0.8, 0           # marker size (pt), CI line width
TXT = 6.0                            # annotation text size (pt)
TXT_S = 5.5
DARK = "#262626"

# ---- figure canvas in mm. save() uses bbox_inches="tight" (pad 0.1 in per side), so the canvas
# is drawn 2 x 2.54 mm smaller than the target and an invisible frame fixes the tight bbox to it.
PAD = 2.54
FW, FH = 183 - 2 * PAD, 139 - 2 * PAD
fig = plt.figure(figsize=(FW * MM, FH * MM))
fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure, fill=False, lw=0))


def ax_at(left, top, width, height):
    """Axes from mm coordinates measured from the top-left corner of the canvas."""
    return fig.add_axes([left / FW, (FH - top - height) / FH, width / FW, height / FH])


def fig_text(x, y_top, s, **kw):
    fig.text(x / FW, (FH - y_top) / FH, s, **kw)


def letter(x, y_top, s):
    fig_text(x, y_top, s, fontsize=8, fontweight="bold", va="baseline", ha="left")


def log_axis(ax, lo, hi, ticks):
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.axvline(1, color="#8C8C8C", lw=0.6, ls=(0, (3, 2)), zorder=0)


def colour(d):
    return DISEASE_COLOR[d]


def forest_point(ax, x, lo, hi, y, c, marker, filled=True, ms=MS):
    ax.plot([lo, hi], [y, y], color=c, lw=LW, solid_capstyle="butt", zorder=2)
    # the pale H61 grey (style.py) gets a darker marker edge, as in Fig. 3, so it stays legible
    edge = "#7F7F7F" if c == DISEASE_COLOR["nco_earwax"] else c
    ax.plot([x], [y], marker=marker, ms=ms, color=c, mec=edge, mew=0.8,
            mfc=c if filled else "white", ls="none", zorder=3)


def fmt_p(p, n):
    floor = 2 / n
    if p < floor - 1e-12:
        return f"<{floor:g}"
    return f"{p:.3f}" if n >= 500 else f"{p:.2f}"


# ---- layout (mm, from top-left) ----------------------------------------------------------------
LAB_R = 20.5                      # right edge of the row-label column
A_L, B_L, C_L = 22.0, 71.0, 127.0
AB_W, C_W = 45.0, 38.0
TOP_Y, TOP_H = 15.0, 66.0         # top block axes: top edge, height
D_TOP, D_H = 103.0, 19.0

ax_a = ax_at(A_L, TOP_Y, AB_W, TOP_H)
ax_b = ax_at(B_L, TOP_Y, AB_W, TOP_H)
ax_c = ax_at(C_L, TOP_Y, C_W, TOP_H)

# ---- panels a and b: group-specific RRs -------------------------------------------------------
AB_LIM, AB_TICKS = (0.40, 1.60), [0.5, 0.7, 1, 1.4]
for ax, model in [(ax_a, "full_cohort"), (ax_b, "separate")]:
    log_axis(ax, *AB_LIM, AB_TICKS)
    sub = grp[grp.model == model]
    for _, r in sub.iterrows():
        y = YPOS[r.disease] + GROUP_OFF[r.cohort]
        forest_point(ax, r.rr, r.lower, r.upper, y, colour(r.disease), GROUP_MARKER[r.cohort])
    ax.set_ylim(*YLIM)
    ax.set_xlabel("RR, sleep→MVPA, 30 min/day", labelpad=2)

for ax in (ax_a, ax_b, ax_c):
    ax.set_yticks([YPOS[d] for d in ROWS])
    ax.tick_params(axis="y", length=0)
for ax in (ax_b, ax_c):
    ax.set_yticklabels([])
ax_a.set_yticklabels([ROW_LABEL[d] for d in ROWS], fontsize=6.5, linespacing=1.0)
for t in ax_a.get_yticklabels():  # negative-control names in the same colour as the diseases (Figs 1c, 3)
    t.set_color(DARK)
ax_a.tick_params(axis="y", pad=1.5)

# negative-control block header (label column) and the 'not fitted' note in b
y_gap = (YPOS["ibs"] + YPOS["nco_cyst"]) / 2
ax_a.text(-1.5 / AB_W, y_gap, "Negative control\noutcomes", transform=ax_a.get_yaxis_transform(),
          ha="right", va="center", fontsize=TXT_S, style="italic", color=GREY, linespacing=1.1)
ax_b.text(0.5, (YPOS["nco_cyst"] + YPOS["nco_earwax"]) / 2,
          "Not fitted for negative control outcomes", transform=ax_b.get_yaxis_transform(),
          ha="center", va="center", fontsize=TXT_S, style="italic", color=GREY)

# ---- panel c: long vs normal ratio, full-cohort and separate models ---------------------------
C_LIM, C_TICKS = (0.78, 1.30), [0.8, 0.9, 1, 1.1, 1.25]
log_axis(ax_c, *C_LIM, C_TICKS)
ln = rat[rat.contrast == "long_vs_normal"]
P_DEC = 1.0 + 4.8 / C_W            # P column decimal point, axes fraction (2.5 mm gap before "<0")
for _, r in ln[ln.model.isin(["full_cohort", "separate"])].iterrows():
    y = YPOS[r.disease] + MODEL_OFF[r.model]
    forest_point(ax_c, r.ratio, r.lower, r.upper, y, colour(r.disease), "o",
                 filled=(r.model == "full_cohort"))
    p_int, p_frac = fmt_p(r.p_boot, r.n_boot).split(".")
    for s_, ha_ in [(p_int, "right"), ("." + p_frac, "left")]:
        ax_c.text(P_DEC, y, s_, transform=ax_c.get_yaxis_transform(), ha=ha_, va="center",
                  fontsize=TXT_S, color=DARK, clip_on=False)
ax_c.text(P_DEC + 1.0 / C_W, YLIM[1], "$P$", transform=ax_c.get_yaxis_transform(), ha="center", va="bottom",
          fontsize=6, color=DARK, clip_on=False)
ax_c.set_ylim(*YLIM)
ax_c.set_xlabel("Ratio of RRs, >9 h vs 7–9 h", labelpad=2)
# degrees of freedom written out: "df" is not defined in the Fig. 4 legend
ax_c.text(0.5, y_gap, f"Between-disease heterogeneity\n(full-cohort model):\n$Q$ = {het.Q:.2f} "
          f"({int(het.df)} degrees of freedom),\n$P$ = {het.p:.2f}", transform=ax_c.get_yaxis_transform(),
          ha="center", va="center",
          fontsize=TXT_S, color=DARK, linespacing=1.15, zorder=4,
          bbox=dict(facecolor="white", edgecolor="none", pad=0.6))

# direction cues under the ratio axes
def cues(ax, left, right, y=-0.30):
    ax.text(0.0, y, left, transform=ax.transAxes, ha="left", va="top", fontsize=TXT_S, color=GREY)
    ax.text(1.0, y, right, transform=ax.transAxes, ha="right", va="top", fontsize=TXT_S, color=GREY)


# ---- titles, letters and in-panel legends -------------------------------------------------------
TITLE_Y = 5.0
letter(0.0, TITLE_Y, "a")
fig_text(4.0, TITLE_Y, "Full-cohort model", fontsize=7, va="baseline", ha="left")
letter(B_L - 4.0, TITLE_Y, "b")
fig_text(B_L, TITLE_Y, "Separate models per sleep-duration group", fontsize=7, va="baseline", ha="left")
letter(C_L - 4.0, TITLE_Y, "c")
fig_text(C_L, TITLE_Y, "Long vs normal sleep", fontsize=7, va="baseline", ha="left")

grp_handles = [Line2D([], [], marker=GROUP_MARKER[g], ms=MS, color=GREY, ls="-", lw=LW)
               for g in ["short", "normal", "long"]]
grp_labels = ["Short sleep (<7 h/day)", "Normal sleep (7–9 h/day)", "Long sleep (>9 h/day)"]
fig.legend(grp_handles, grp_labels, loc="center", ncol=3, fontsize=TXT,
           bbox_to_anchor=((A_L + B_L + AB_W) / 2 / FW, (FH - 10.0) / FH),
           handlelength=1.6, handletextpad=0.4, columnspacing=1.6, borderaxespad=0, borderpad=0)
mod_handles = [Line2D([], [], marker="o", ms=MS, color=GREY, mfc=GREY, lw=LW),
               Line2D([], [], marker="o", ms=MS, color=GREY, mfc="white", mew=0.8, lw=LW)]
fig.legend(mod_handles, ["Full-cohort model", "Separate models"], loc="center left", ncol=2,
           fontsize=TXT, bbox_to_anchor=(C_L / FW, (FH - 10.0) / FH), handlelength=1.6,
           handletextpad=0.4, columnspacing=1.2, borderaxespad=0, borderpad=0)

# ---- panel d: gallstone long vs short, three settings ------------------------------------------
D_L, D_W = 34.0, 95.0
ax_d = ax_at(D_L, D_TOP, D_W, D_H)
D_LIM, D_TICKS = (0.80, 2.0), [0.8, 1, 1.25, 1.6, 2]
log_axis(ax_d, *D_LIM, D_TICKS)
settings = [("full_cohort", "Full-cohort model\n(primary)", "o", True),
            ("main_effects", "Main-effects model\n(no interaction terms)", "D", True),
            ("separate", "Separate models per\nsleep-duration group", "o", False)]
ls_ = rat[(rat.contrast == "long_vs_short") & (rat.disease == "gallstone")].set_index("model")
c_gs = colour("gallstone")
for i, (m, lab, mk, filled) in enumerate(settings):
    r = ls_.loc[m]
    forest_point(ax_d, r.ratio, r.lower, r.upper, i, c_gs, mk, filled=filled,
                 ms=MS if mk == "o" else MS - 0.6)
    ax_d.text(1.0 + 3.0 / D_W, i, f"{r.ratio:.2f} ({r.lower:.2f}–{r.upper:.2f})",
              transform=ax_d.get_yaxis_transform(), ha="left", va="center", fontsize=TXT, color=DARK)
ax_d.text(1.0 + 3.0 / D_W, -0.62, "Ratio (95% CI)", transform=ax_d.get_yaxis_transform(),
          ha="left", va="bottom", fontsize=TXT, color=DARK)
ax_d.set_ylim(len(settings) - 0.5, -0.5)
ax_d.set_yticks(range(len(settings)))
ax_d.set_yticklabels([s[1] for s in settings], fontsize=6.5, linespacing=1.0)
ax_d.tick_params(axis="y", length=0, pad=1.5)
ax_d.set_xlabel("Ratio of RRs, >9 h vs <7 h", labelpad=2)
for ax in (ax_a, ax_b, ax_c, ax_d):  # forest-style axes without a left spine, as in Figs 2a, 3 and 5c
    ax.spines["left"].set_visible(False)
D_TITLE_Y = D_TOP - 6.5
letter(0.0, D_TITLE_Y, "d")
fig_text(4.0, D_TITLE_Y, "Gallstone disease: long vs short sleep", fontsize=7, va="baseline", ha="left")

cues(ax_c, "← Lower in >9 h", "Lower in 7–9 h →", y=-8.8 / TOP_H)
cues(ax_d, "← Lower in >9 h", "Lower in <7 h →", y=-8.8 / D_H)

save(fig, "FigS6_sleep_groups")

# Flatten the raster exports to RGB on white (production expects RGB TIFF); style.save() is shared and left as is.
from PIL import Image

for ext, kw in [("tiff", {"compression": "tiff_lzw", "dpi": (600, 600)}), ("png", {"dpi": (300, 300)})]:
    p = FIGDIR / f"FigS6_sleep_groups.{ext}"
    with Image.open(p) as im:
        im.load()
        rgba = im.convert("RGBA")
        rgb = Image.new("RGB", im.size, (255, 255, 255))
        rgb.paste(rgba, mask=rgba.split()[-1])
    rgb.save(p, **kw)
print("saved")
