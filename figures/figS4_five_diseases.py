"""Supplementary Fig. S4 | Risk ratios, risk differences and dose-response across the five diseases.

Primary analyses, full cohort (results/primary_range_thr for gallstone disease,
results/outcomes/<o>/primary_range_thr for the others; 500 bootstrap replicates).
Data: figures/data/fig2_full_cohort.csv, written by figures/data/fig2_data.R.

a  30-min RR (95% CI), five diseases x six reallocation directions, shared log x-axis.
b  RD per 1,000 (95% CI) for sleep->MVPA and MVPA->sleep at 30 min, reference risk annotated.
c  Sleep<->MVPA dose-response, -60 to +60 min/day, log y-axis.
"""
from style import *  # noqa: F401

import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
from matplotlib.transforms import blended_transform_factory

df = pd.read_csv(DATADIR / "fig2_full_cohort.csv")
assert (df["cohort"] == "full").all() and df["reported"].all() and (df["n_boot"] == 500).all()

DIRS = [("sleep_to_mvpa", "Sleep→MVPA"), ("mvpa_to_sleep", "MVPA→sleep"),
        ("sleep_to_light", "Sleep→light"), ("light_to_sleep", "Light→sleep"),
        ("sleep_to_inactivity", "Sleep→inactivity"), ("inactivity_to_sleep", "Inactivity→sleep")]
REF_LS = (0, (3, 2))          # dashed reference line
REF_COL = "#8C8C8C"
TXT = "#272727"


def get(dis, direction, minutes):
    r = df[(df.disease == dis) & (df.direction == direction) & (df.minutes == minutes)]
    assert len(r) == 1, (dis, direction, minutes)
    return r.iloc[0]


def tint(c, f=0.55):
    """Mix colour c with white (f = share of white)."""
    r, g, b = mcolors.to_rgb(c)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


# ---------------------------------------------------------------- layout (mm)
FW, FH = 183.0, 169.0
fig = plt.figure(figsize=(FW * MM, FH * MM))


def add_ax(left, top, width, height):
    """Axes from mm coordinates measured from the figure's top-left corner."""
    return fig.add_axes([left / FW, (FH - top - height) / FH, width / FW, height / FH])


def letter(ax, ch, left_mm, top_mm, ax_left, ax_top, ax_w, ax_h):
    """Place a panel letter at an absolute mm position via style.panel_label."""
    panel_label(ax, ch, x=(left_mm - ax_left) / ax_w, y=(ax_top - top_mm) / ax_h + 1)


L, R = 26.0, 178.5            # common left/right edge of the data area (export stays <= 183 mm with tight-bbox pad)

# ---------------------------------------------------------------- a: RR grid
A_TOP, A_H = 10.0, 33.0
g_in, g_out = 2.2, 5.0
wa = (R - L - 3 * g_in - 2 * g_out) / 6
a_lefts, x = [], L
for j in range(6):
    a_lefts.append(x)
    x += wa + (g_in if j % 2 == 0 else g_out)

XLIM = (0.72, 1.62)
XT = [0.8, 1.0, 1.25, 1.5]
XTL = ["0.8", "1", "1.25", "1.5"]
ypos = {d: 4 - i for i, d in enumerate(DISEASES)}
axes_a = []
for j, (key, lab) in enumerate(DIRS):
    ax = add_ax(a_lefts[j], A_TOP, wa, A_H)
    axes_a.append(ax)
    ax.set_xscale("log")
    ax.set_xlim(*XLIM)
    ax.xaxis.set_major_locator(FixedLocator(XT))
    ax.xaxis.set_major_formatter(FixedFormatter(XTL))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis="x", pad=1.5)
    ax.set_ylim(-0.6, 4.6)
    ax.axvline(1, color=REF_COL, lw=0.6, ls=REF_LS, zorder=0)
    for d in DISEASES:
        r = get(d, key, 30)
        y = ypos[d]
        # CI as a thick tinted bar, RR as a thin vertical tick, so narrow CIs are not hidden by a marker
        ax.plot([r.lower, r.upper], [y, y], color=tint(DISEASE_COLOR[d], 0.35), lw=2.2,
                solid_capstyle="butt", zorder=2)
        ax.plot(r.rr, y, "|", ms=5.0, mec=DISEASE_COLOR[d], mew=1.1, zorder=3)
    ax.spines["left"].set_visible(False)
    ax.set_title(lab, fontsize=7, pad=4)
    if j == 0:
        ax.set_yticks([ypos[d] for d in DISEASES])
        ax.set_yticklabels([DISEASE_LABEL[d] for d in DISEASES], fontsize=6.5)
        ax.tick_params(axis="y", length=0, pad=2)
    else:
        ax.set_yticks([])
fig.text((L + R) / 2 / FW, 1 - (A_TOP + A_H + 7.2) / FH, "RR (95% CI), 30 min/day",
         ha="center", va="center", fontsize=7)
letter(axes_a[0], "a", 1.0, 3.5, a_lefts[0], A_TOP, wa, A_H)

# ---------------------------------------------------------------- c geometry (b is aligned to it)
C_TOP, C_H = 120.0, 33.0
g_c = 4.0
wc = (R - L - 4 * g_c) / 5
c_lefts = [L + i * (wc + g_c) for i in range(5)]
c_centres = [cl + wc / 2 - L for cl in c_lefts]      # in b's data units (mm from L)

# ---------------------------------------------------------------- b: RD per 1,000
B_TOP, B_H = 66.0, 38.0
axb = add_ax(L, B_TOP, R - L, B_H)
axb.set_xlim(0, R - L)
axb.set_ylim(-9.0, 13.0)
axb.set_yticks([-5, 0, 5, 10])
axb.axhline(0, color="#4D4D4D", lw=0.6, zorder=1)
axb.spines["bottom"].set_visible(False)
bw, off = 5.2, 3.0
for d, cx in zip(DISEASES, c_centres):
    col = DISEASE_COLOR[d]
    for key, dx, fc in [("sleep_to_mvpa", -off, col), ("mvpa_to_sleep", +off, tint(col))]:
        r = get(d, key, 30)
        axb.bar(cx + dx, r.rd_per_1000, width=bw, color=fc, edgecolor=col, linewidth=0.6, zorder=2)
        axb.errorbar(cx + dx, r.rd_per_1000,
                     yerr=[[r.rd_per_1000 - r.rd_lower], [r.rd_upper - r.rd_per_1000]],
                     fmt="none", ecolor=TXT, elinewidth=0.6, capsize=1.3, capthick=0.6, zorder=3)
    ref = get(d, "sleep_to_mvpa", 30).ref_per_1000
    assert abs(ref - get(d, "mvpa_to_sleep", 30).ref_per_1000) < 1e-9
    axb.text(cx, 1.03, f"{ref:.1f}", transform=blended_transform_factory(axb.transData, axb.transAxes),
             ha="center", va="bottom", fontsize=6, color=GREY)
axb.set_xticks(c_centres)
axb.set_xticklabels([DISEASE_LABEL[d] for d in DISEASES], fontsize=6.5)
axb.tick_params(axis="x", length=0, pad=3)
axb.set_ylabel("RD per 1,000 over\n9.8 years (95% CI)", labelpad=2, linespacing=1.1)
fig.text((L - 2.0) / FW, 1 - (B_TOP - 3.2) / FH, "Reference risk\nper 1,000", ha="right", va="center",
         fontsize=6, color=GREY, linespacing=1.1)
leg = axb.legend(handles=[Patch(fc="#6E6E6E", ec="#6E6E6E", lw=0.6, label="Sleep→MVPA"),
                          Patch(fc=tint("#6E6E6E"), ec="#6E6E6E", lw=0.6, label="MVPA→sleep")],
                 loc="upper right", bbox_to_anchor=(1.0, 1.0), fontsize=6.5, handlelength=1.2,
                 handleheight=0.9, handletextpad=0.5, borderaxespad=0.3, labelspacing=0.35)
letter(axb, "b", 1.0, B_TOP - 8.0, L, B_TOP, R - L, B_H)

# ---------------------------------------------------------------- c: dose-response
YLIM = (0.55, 2.5)
YT = [0.6, 0.8, 1.0, 1.5, 2.0]
YTL = ["0.6", "0.8", "1", "1.5", "2"]
axes_c = []
for i, d in enumerate(DISEASES):
    ax = add_ax(c_lefts[i], C_TOP, wc, C_H)
    axes_c.append(ax)
    col = DISEASE_COLOR[d]
    # duration < 0: sleep->MVPA; duration > 0: MVPA->sleep; 0 = natural course (RR 1)
    sub = df[(df.disease == d) & (df["from"] == "mvpa")].sort_values("duration")
    neg, pos = sub[sub.duration < 0], sub[sub.duration > 0]
    xs = np.concatenate([neg.duration, [0], pos.duration])
    rr = np.concatenate([neg.rr, [1.0], pos.rr])
    lo = np.concatenate([neg.lower, [1.0], pos.lower])
    hi = np.concatenate([neg.upper, [1.0], pos.upper])
    assert list(xs) == [-60, -45, -30, -15, 0, 15, 30, 45, 60]
    ax.set_yscale("log")
    ax.set_ylim(*YLIM)
    ax.set_xlim(-68, 68)
    ax.axhline(1, color=REF_COL, lw=0.6, ls=REF_LS, zorder=0)
    ax.fill_between(xs, lo, hi, color=col, alpha=0.22, lw=0, zorder=1)
    ax.plot(xs, rr, "-", color=col, lw=1.0, zorder=2)
    ax.plot(xs, rr, "o", ms=2.8, mfc=col, mec="white", mew=0.35, zorder=3)
    ax.set_xticks([-60, -30, 0, 30, 60])
    ax.set_xticklabels(["−60", "−30", "0", "30", "60"])
    ax.tick_params(axis="x", pad=1.5)
    ax.yaxis.set_major_locator(FixedLocator(YT))
    ax.yaxis.set_minor_locator(NullLocator())
    if i == 0:
        ax.yaxis.set_major_formatter(FixedFormatter(YTL))
        ax.tick_params(axis="y", pad=1.5)
        ax.set_ylabel("RR (95% CI)", labelpad=2)
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.set_title(DISEASE_LABEL[d], fontsize=7, pad=3)
    if i == 0:  # sign convention is shared by all five panels; label it once
        tr = ax.get_xaxis_transform()
        ax.text(-34, -0.165, "Sleep→MVPA", transform=tr, ha="center", va="top", fontsize=5.5, color=GREY)
        ax.text(34, -0.165, "MVPA→sleep", transform=tr, ha="center", va="top", fontsize=5.5, color=GREY)
fig.text((L + R) / 2 / FW, 1 - (C_TOP + C_H + 11.0) / FH, "Change in sleep (min/day)",
         ha="center", va="center", fontsize=7)
letter(axes_c[0], "c", 1.0, C_TOP - 6.0, c_lefts[0], C_TOP, wc, C_H)

save(fig, "FigS4_five_diseases")

# Flatten the raster exports to RGB (production expects RGB TIFF); style.save() is shared and left as is.
from PIL import Image

for ext, kw in [("tiff", {"compression": "tiff_lzw", "dpi": (600, 600)}), ("png", {"dpi": (300, 300)})]:
    p = FIGDIR / f"FigS4_five_diseases.{ext}"
    with Image.open(p) as im:
        rgb = im.convert("RGB")
    rgb.save(p, **kw)
print("saved")
