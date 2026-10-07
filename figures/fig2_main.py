"""Fig. 2 | Reallocating sleep and risk of gallstone disease, and disease specificity.

a  Gallstone disease: dose-response for reallocating 15-60 min/day between sleep and MVPA, light
   activity or inactivity (full cohort, 500 bootstrap replicates), log y-axis.
b  Gallstone disease: RD per 1,000 over 9.8 years for sleep -> MVPA, all diagnoses vs the
   cholecystectomy endpoint (200 replicates).
c  30-min sleep -> MVPA and sleep -> light RR for gallstone disease (primary outcome), the four
   secondary outcomes and the two negative control outcomes; right column: within-replicate ratio
   RR(sleep -> MVPA) / RR(sleep -> light).

Data: figures/data/fig2_full_cohort.csv (figures/data/fig2_data.R); figures/data/fig2m_surgery.csv and
figures/data/fig2m_compare.csv (figures/data/fig2_main_data.R).
"""
from style import *  # noqa: F403

import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
from matplotlib import transforms

full = pd.read_csv(DATADIR / "fig2_full_cohort.csv")
surg = pd.read_csv(DATADIR / "fig2m_surgery.csv")
cmp = pd.read_csv(DATADIR / "fig2m_compare.csv")
MINUS = "−"
REF = dict(color="#8C8C8C", lw=0.6, ls=(0, (3, 2)), zorder=0)
TXT = "#1A1A1A"
LINE_COL = {"mvpa": BEHAVIOUR_COLOR["mvpa"], "light": BEHAVIOUR_COLOR["light"], "inactivity": "#7A7A7A"}


def tint(c, f=0.55):
    r, g, b = mcolors.to_rgb(c)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def m(x):
    return f"{x:g}".replace("-", MINUS)


FW, FH = 183.0, 152.0
fig = plt.figure(figsize=(FW * MM, FH * MM))


def ax_mm(left, top, w, h):
    return fig.add_axes([left / FW, 1 - (top + h) / FH, w / FW, h / FH])


# ------------------------------------------------------------------ a: gallstone dose-response
A_L, A_T, A_W, A_H = 15.0, 9.0, 82.0, 52.0
axa = ax_mm(A_L, A_T, A_W, A_H)
axa.set_yscale("log")
axa.axhline(1, **REF)
axa.axvline(0, color="#BDBDBD", lw=0.5, zorder=0)
ends = {}
for b in ["inactivity", "light", "mvpa"]:
    d = full[(full.disease == "gallstone") & (full["from"] == b)].sort_values("duration")
    xs = np.r_[d.duration[d.duration < 0], 0, d.duration[d.duration > 0]]
    rr = np.r_[d.rr[d.duration < 0], 1.0, d.rr[d.duration > 0]]
    lo = np.r_[d.lower[d.duration < 0], 1.0, d.lower[d.duration > 0]]
    hi = np.r_[d.upper[d.duration < 0], 1.0, d.upper[d.duration > 0]]
    col = LINE_COL[b]
    axa.fill_between(xs, lo, hi, color=col, alpha=0.18 if b == "mvpa" else 0.14, lw=0, zorder=1)
    axa.plot(xs, rr, "-", color=col, lw=1.3 if b == "mvpa" else 1.0, zorder=3)
    axa.plot(xs, rr, "o", ms=3.0 if b == "mvpa" else 2.5, mfc=col, mec="white", mew=0.35, zorder=4)
    ends[b] = (rr[0], rr[-1])
axa.set_ylim(0.56, 1.95)
yt = [0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75]
axa.yaxis.set_major_locator(FixedLocator(yt))
axa.yaxis.set_major_formatter(FixedFormatter([m(v) for v in yt]))
axa.yaxis.set_minor_locator(NullLocator())
axa.set_xlim(-66, 66)
axa.set_xticks([-60, -45, -30, -15, 0, 15, 30, 45, 60])
axa.set_xticklabels([m(v) for v in [-60, -45, -30, -15, 0, 15, 30, 45, 60]])
axa.set_xlabel("Change in sleep (min/day)")
axa.set_ylabel("RR for gallstone disease (95% CI)", labelpad=2)
tr = axa.get_xaxis_transform()
axa.text(-33, 1.02, "Less sleep, more of the other behaviour", transform=tr, ha="center", va="bottom",
         fontsize=5.5, color=GREY)
axa.text(33, 1.02, "More sleep, less of the other behaviour", transform=tr, ha="center", va="bottom",
         fontsize=5.5, color=GREY)
# legend in the empty lower-right quadrant
lab = {"mvpa": "MVPA", "light": "Light activity", "inactivity": "Inactivity"}
axa.legend(handles=[Line2D([], [], color=LINE_COL[k], lw=1.3 if k == "mvpa" else 1.0, marker="o", ms=3,
                           mfc=LINE_COL[k], mec="white", mew=0.35, label=lab[k]) for k in ["mvpa", "light", "inactivity"]],
           title="Behaviour exchanged with sleep", title_fontsize=6, loc="lower right", fontsize=6,
           handlelength=2.0, borderaxespad=0.6, labelspacing=0.35)
axa.set_title("Gallstone disease", fontsize=7, pad=12, fontweight="bold")
panel_label(axa, "a", x=-0.155, y=1.07)

# ------------------------------------------------------------------ b: RD, diagnosis vs cholecystectomy
B_L, B_T, B_W, B_H = 123.0, 9.0, 57.0, 52.0
axb = ax_mm(B_L, B_T, B_W, B_H)
dg = full[(full.disease == "gallstone") & (full.direction == "sleep_to_mvpa")].set_index("minutes").sort_index()
sg = surg[surg.direction == "sleep_to_mvpa"].set_index("minutes").sort_index()
mins = [15, 30, 45, 60]
xw, off = 0.34, 0.19
col = DISEASE_COLOR["gallstone"]
for i, mn in enumerate(mins):
    for src, dx, fc in [(dg, -off, col), (sg, off, "white")]:
        r = src.loc[mn]
        axb.bar(i + dx, r.rd_per_1000, width=xw, color=fc, edgecolor=col, lw=0.8, zorder=2)
        axb.errorbar(i + dx, r.rd_per_1000, yerr=[[r.rd_per_1000 - r.rd_lower], [r.rd_upper - r.rd_per_1000]],
                     fmt="none", ecolor=TXT, elinewidth=0.6, capsize=1.3, capthick=0.6, zorder=3)
axb.axhline(0, color="#4D4D4D", lw=0.6, zorder=1)
axb.set_xticks(range(4))
axb.set_xticklabels([str(v) for v in mins])
axb.tick_params(axis="x", length=0, pad=3)
axb.spines["bottom"].set_visible(False)
axb.set_ylim(-13.6, 0.8)  # room below the bars for the legend
axb.set_yticks([-12, -10, -8, -6, -4, -2, 0])
axb.set_yticklabels([m(v) for v in [-12, -10, -8, -6, -4, -2, 0]])
axb.set_xlabel("Sleep reallocated to MVPA (min/day)")
axb.set_ylabel("RD per 1,000 over 9.8 years (95% CI)", labelpad=2)
ref_d, ref_s = dg.loc[30].ref_per_1000, sg.loc[30].ref_per_1000
axb.legend(handles=[Patch(fc=col, ec=col, lw=0.8, label=f"All diagnoses (reference risk {ref_d:.1f} per 1,000)"),
                    Patch(fc="white", ec=col, lw=0.8,
                          label=f"Cholecystectomy endpoint (reference risk {ref_s:.1f} per 1,000)")],
           loc="lower left", fontsize=5.5, handlelength=1.1, handleheight=0.8, handletextpad=0.5,
           borderaxespad=0.2, labelspacing=0.3)
axb.set_title("Gallstone disease: risk difference", fontsize=7, pad=12, fontweight="bold")
panel_label(axb, "b", x=-0.30, y=1.07)

# ------------------------------------------------------------------ c: disease specificity forest
ROWS = [("gallstone", "Gallstone disease (primary outcome)"), ("peptic_ulcer", "Peptic ulcer"),
        ("gerd", "GERD"), ("diverticular", "Diverticular disease"), ("ibs", "IBS"), None,
        ("nco_cyst", "Follicular cyst (L72)"), ("nco_earwax", "External ear disorders (H61)")]
C_T, ROW = 84.0, 5.6
ypos, y = {}, 0.0
for r in ROWS:
    if r is None:
        y += 0.6
        continue
    if r[0] == "peptic_ulcer":  # room for the "Secondary outcomes" header
        y += 0.6
    ypos[r[0]] = y
    y += 1
nrow = y
C_H = nrow * ROW
X_LAB = 2.0
P1_L, P_W, P_G = 60.0, 38.0, 6.0
P2_L = P1_L + P_W + P_G
X_RAT = P2_L + P_W + 5.0
axes_c = []
for k, (left, frm, title) in enumerate([(P1_L, "mvpa", "Sleep→MVPA, 30 min/day"),
                                         (P2_L, "light", "Sleep→light activity, 30 min/day")]):
    ax = ax_mm(left, C_T, P_W, C_H)
    axes_c.append(ax)
    ax.set_xscale("log")
    ax.set_xlim(0.73, 1.18)
    xt = [0.8, 0.9, 1.0, 1.1]
    ax.xaxis.set_major_locator(FixedLocator(xt))
    ax.xaxis.set_major_formatter(FixedFormatter([m(v) for v in xt]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_ylim(nrow - 0.4, -0.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.patch.set_visible(False)  # let the primary-outcome band run across both forests
    ax.axvline(1, **REF)
    ax.set_title(title, fontsize=6.5, pad=3)
    ax.set_xlabel("RR (95% CI)")
    for o, yy in ypos.items():
        r = cmp[(cmp.outcome == o) & (cmp["from"] == frm)].iloc[0]
        c = DISEASE_COLOR[o]
        prim = o == "gallstone"
        ax.plot([r.lower, r.upper], [yy, yy], color=c, lw=1.4 if prim else 0.9, solid_capstyle="butt", zorder=2)
        ax.plot(r.rr, yy, marker="s", ms=4.2 if prim else 3.4, color=c,
                mec="#7F7F7F" if o == "nco_earwax" else c, mew=0.4 if o == "nco_earwax" else 0, zorder=3)
trf = transforms.blended_transform_factory(fig.transFigure, axes_c[0].transData)
# highlight band for the primary outcome and labels
axes_c[0].add_patch(Rectangle((X_LAB / FW - 0.005, ypos["gallstone"] - 0.5), (FW - 2 * X_LAB) / FW + 0.005, 1.0,
                              transform=trf, facecolor=tint(DISEASE_COLOR["gallstone"], 0.88), edgecolor="none",
                              zorder=-1, clip_on=False))
for r in ROWS:
    if r is None:
        continue
    o, name = r
    axes_c[0].text(X_LAB / FW, ypos[o], name, transform=trf, ha="left", va="center", fontsize=6.5,
                   fontweight="bold" if o == "gallstone" else "normal", color=TXT, clip_on=False)
    q = cmp[(cmp.outcome == o)].iloc[0]
    axes_c[0].text(X_RAT / FW, ypos[o], f"{q.ratio:.2f} ({q.ratio_lower:.2f}–{q.ratio_upper:.2f})", transform=trf,
                   ha="left", va="center", fontsize=6.5, color=TXT, clip_on=False)
axes_c[0].text(X_LAB / FW, ypos["peptic_ulcer"] - 0.8, "Secondary outcomes", transform=trf, ha="left",
               va="center", fontsize=6, style="italic", color=GREY, clip_on=False)
axes_c[0].text(X_LAB / FW, ypos["nco_cyst"] - 0.8, "Negative control outcomes", transform=trf, ha="left",
               va="center", fontsize=6, style="italic", color=GREY, clip_on=False)
fig.text(X_RAT / FW, 1 - (C_T - 3.2) / FH, "RR ratio, MVPA vs light", ha="left", va="bottom", fontsize=6.5,
         fontweight="bold", color=TXT)
fig.text(X_RAT / FW, 1 - (C_T - 0.2) / FH, "(95% CI)", ha="left", va="bottom", fontsize=6.5, color=TXT)
# how to read the ratio, under the column
fig.text(X_RAT / FW, 1 - (C_T + C_H + 3.0) / FH, "<1: lower risk when sleep is\nreplaced by MVPA than by\nlight activity",
         ha="left", va="top", fontsize=5.5, color=GREY, linespacing=1.25)
fig.text(1.0 / FW, 1 - (C_T - 5.0) / FH, "c", fontsize=8, fontweight="bold", va="bottom", ha="left")

save(fig, "Fig2_gallstone_main")
from PIL import Image

for ext, kw in [("tiff", {"compression": "tiff_lzw", "dpi": (600, 600)}), ("png", {"dpi": (300, 300)})]:
    p = FIGDIR / f"Fig2_gallstone_main.{ext}"
    with Image.open(p) as im:
        im.load()
        rgba = im.convert("RGBA")
        rgb = Image.new("RGB", im.size, (255, 255, 255))
        rgb.paste(rgba, mask=rgba.split()[-1])
    rgb.save(p, **kw)
print("saved")
