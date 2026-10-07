"""Fig. 4 | Reallocating sleep and abdominal fat (MRI subsample), and fat depots and subsequent disease.

a  Percent difference in VAT, ASAT and liver PDFF when 15-60 min/day is reallocated between sleep and MVPA,
   light activity or inactivity (primary model; lines and 95% CI bands, drawn as in Fig. 2a), and sleep <-> MVPA
   further adjusted for BMI at imaging (dashed line, open markers, 95% CI bars).
   VAT and ASAT: mean difference as a percentage of the reference mean (the mean under the observed
   compositions: 3.58 L and 6.61 L). PDFF is analysed on the log scale, so its difference is the percent
   difference in the geometric mean. results/mri_primary and results/mri_bmi (mean_differences.csv), 500
   bootstrap replicates.
b  HR per SD of VAT and of log PDFF for each disease after imaging: primary Cox model (filled) and further
   adjusted for BMI (open) (figures/data/fig4_fatcox.csv). The VAT model mutually adjusted for ASAT is reported in
   the text and in Additional file 1: Table S6.
"""
from style import *  # noqa: F403

import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
from matplotlib.transforms import blended_transform_factory

cox = pd.read_csv(DATADIR / "fig4_fatcox.csv")
MD = {mod: pd.read_csv(PROJ / f"results/mri_{mod}/mean_differences.csv") for mod in ("primary", "bmi")}

FW, FH = 183.0, 151.0  # mm
fig = plt.figure(figsize=(FW * MM, FH * MM))


def ax_mm(left, top, w, h):
    return fig.add_axes([left / FW, 1 - (top + h) / FH, w / FW, h / FH])


MINUS = "−"
TEXT = "#1A1A1A"
REF = dict(color="#8C8C8C", lw=0.6, ls=(0, (3, 2)), zorder=0)
BCOL = {"mvpa": BEHAVIOUR_COLOR["mvpa"], "light": BEHAVIOUR_COLOR["light"], "inactivity": "#7A7A7A"}  # as Fig. 2a
BLAB = {"mvpa": "MVPA", "light": "Light activity", "inactivity": "Inactivity"}
BMI_LS = (0, (2.5, 1.5))


def m(v):
    return f"{v:g}".replace("-", MINUS)


def series(model, outcome, frm):
    """Durations -60..60 with the natural course (0) inserted; estimate and 95% CI in percent."""
    d = MD[model]
    d = d[(d.cohort == "full") & (d.outcome == outcome) & (d["from"] == frm)].sort_values("duration")
    assert len(d) == 8 and d.reported.all() and (d.duration.values[:4] < 0).all(), (model, outcome, frm)
    if outcome == "log_pdff":
        pct = lambda v: 100 * (np.exp(v.values) - 1)  # noqa: E731
    else:
        pct = lambda v: 100 * v.values / d.ref_mean.values  # noqa: E731
    ins = lambda a, v: np.insert(np.asarray(a, dtype=float), 4, v)  # noqa: E731
    return ins(d.duration.values, 0), ins(pct(d.md), 0), ins(pct(d.lower), 0), ins(pct(d.upper), 0)


# values quoted in the text (30 min/day sleep -> MVPA; index 2 = -30)
for (mod, o), v in {("primary", "vat"): -8.0, ("primary", "asat"): -6.4, ("primary", "log_pdff"): -6.9,
                    ("bmi", "vat"): -2.8, ("bmi", "log_pdff"): -2.2}.items():
    assert round(series(mod, o, "mvpa")[1][2], 1) == v, (mod, o)

# ----------------------------------------------------------------------------------------------
# a: dose-response for the three fat depots
# ----------------------------------------------------------------------------------------------
A_TOP, A_H, A_W = 12.0, 40.0, 49.0
A_LEFTS = [17.0, 73.5, 130.0]
SPECS = [("vat", "Visceral fat (VAT)"), ("asat", "Abdominal subcutaneous fat (ASAT)"),
         ("log_pdff", "Liver fat (PDFF)")]
axes_a = []
for k, (outc, title) in enumerate(SPECS):
    ax = ax_mm(A_LEFTS[k], A_TOP, A_W, A_H)
    axes_a.append(ax)
    ax.axhline(0, **REF)
    ax.axvline(0, color="#BDBDBD", lw=0.5, zorder=0)
    for b in ["inactivity", "light", "mvpa"]:
        xs, est, lo, hi = series("primary", outc, b)
        col = BCOL[b]
        ax.fill_between(xs, lo, hi, color=col, alpha=0.18 if b == "mvpa" else 0.14, lw=0, zorder=1)
        ax.plot(xs, est, "-", color=col, lw=1.3 if b == "mvpa" else 1.0, zorder=3)
        ax.plot(xs, est, "o", ms=3.0 if b == "mvpa" else 2.5, mfc=col, mec="white", mew=0.35, zorder=4)
    xs, est, lo, hi = series("bmi", outc, "mvpa")
    col, nz, dx = BCOL["mvpa"], xs != 0, 2.5  # BMI-adjusted markers nudged right so they do not hide others
    ax.plot(xs + np.where(nz, dx, 0), est, color=col, lw=0.9, ls=BMI_LS, zorder=3)
    ax.errorbar(xs[nz] + dx, est[nz], yerr=[est[nz] - lo[nz], hi[nz] - est[nz]], fmt="none", ecolor=col,
                elinewidth=0.6, capsize=0, zorder=4)
    ax.plot(xs[nz] + dx, est[nz], "o", ms=2.8, mfc="white", mec=col, mew=0.7, zorder=5)
    ax.set_xlim(-67, 67)
    ax.set_xticks([-60, -30, 0, 30, 60])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: m(v)))
    ax.set_xlabel("Change in sleep (min/day)")
    ax.set_ylim(-17.5, 27.5)
    ax.set_yticks([-10, 0, 10, 20])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: m(v)))
    if k == 0:
        ax.set_ylabel("Difference from reference (%)", labelpad=2)
    else:
        ax.tick_params(labelleft=False)
    ax.set_title(title, fontsize=6.5, pad=3)
handles = [Line2D([], [], color=BCOL[b], lw=1.3 if b == "mvpa" else 1.0, marker="o", ms=3, mfc=BCOL[b], mec="white",
                  mew=0.35, label=BLAB[b]) for b in ["mvpa", "light", "inactivity"]]
handles.append(Line2D([], [], color=BCOL["mvpa"], lw=0.9, ls=BMI_LS, marker="o", ms=2.8, mfc="white",
                      mec=BCOL["mvpa"], mew=0.7, label="MVPA, further adjusted for BMI at imaging"))
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1 - (A_TOP + A_H + 9.5) / FH), ncol=4,
           fontsize=6, frameon=False, handlelength=2.4, columnspacing=1.8, handletextpad=0.5,
           title="Behaviour exchanged with sleep", title_fontsize=6)
fig.text(1.5 / FW, 1 - 4.6 / FH, "a", fontsize=8, fontweight="bold", va="bottom", ha="left")
fig.text(6.0 / FW, 1 - 4.6 / FH, "Reallocating sleep and abdominal fat (MRI subsample)", fontsize=7,
         fontweight="bold", va="bottom", ha="left")

# ----------------------------------------------------------------------------------------------
# b: fat depots and disease after imaging (Cox), primary and BMI-adjusted
# ----------------------------------------------------------------------------------------------
B_HEAD = 76.0  # panel heading baseline (mm from top)
DEPOTS = [("VAT", "VAT + BMI", "VAT"), ("log PDFF", "PDFF + BMI", "Liver PDFF")]
rows, y = [], 0.0
for j, o in enumerate(DISEASES):
    y += 0.45 if j else 0.0
    for lab_p, lab_b, name in DEPOTS:
        rp = cox[(cox.disease == o) & (cox.label == lab_p)]
        rb = cox[(cox.disease == o) & (cox.label == lab_b)]
        assert len(rp) == 1 and len(rb) == 1, (o, name)
        rows.append((o, name, rp.iloc[0], rb.iloc[0], y))
        y += 1
nrow = y
ROW = 5.0
B_TOP = B_HEAD + 9.0  # first row band starts here
X_DIS, X_DEP, X_N = 2.5, 36.0, 72.0
AX_L, AX_W = 78.0, 56.0
X_HR1, X_HR2 = AX_L + AX_W + 4.0, AX_L + AX_W + 25.5
axb = ax_mm(AX_L, B_TOP, AX_W, nrow * ROW)
axb.set_ylim(nrow - 0.5, -0.5)
axb.set_xscale("log")
axb.set_xlim(0.6, 2.4)
axb.xaxis.set_major_locator(FixedLocator([0.7, 1.0, 1.5, 2.0]))
axb.xaxis.set_minor_locator(NullLocator())
axb.xaxis.set_major_formatter(FuncFormatter(lambda v, p: "1" if v == 1 else f"{v:g}"))
axb.set_yticks([])
axb.spines["left"].set_visible(False)
axb.axvline(1, **REF)
axb.set_xlabel("HR per SD (95% CI)")
trb = blended_transform_factory(fig.transFigure, axb.transData)


def fx(mm):
    return mm / FW


# background bands: primary outcome tinted, then alternate grey
for j, o in enumerate(DISEASES):
    ys = [r[4] for r in rows if r[0] == o]
    if o == "gallstone":
        fc = tuple(c + (1 - c) * 0.88 for c in mcolors.to_rgb(DISEASE_COLOR["gallstone"]))
    elif j % 2 == 0:
        fc = "#F2F2F2"
    else:
        continue
    axb.add_artist(Rectangle((fx(1.0), min(ys) - 0.5), fx(FW - 2.0), max(ys) - min(ys) + 1.0, transform=trb,
                             facecolor=fc, edgecolor="none", zorder=-5, clip_on=False))
for o, name, rp, rb, yy in rows:
    c = DISEASE_COLOR[o]
    for r, dy, face in [(rp, -0.16, c), (rb, 0.16, "white")]:
        axb.plot([r.lower, r.upper], [yy + dy] * 2, color=c, lw=0.9, solid_capstyle="butt", zorder=2)
        axb.scatter([r.hr], [yy + dy], s=13, marker="s", facecolor=face, edgecolor=c, linewidth=0.8, zorder=3)
    first = name == "VAT"
    if first:
        axb.text(fx(X_DIS), yy, DISEASE_LABEL[o], transform=trb, ha="left", va="center", fontsize=6.5,
                 fontweight="bold", color=TEXT)
    elif o == "gallstone":
        axb.text(fx(X_DIS), yy, "Primary outcome", transform=trb, ha="left", va="center", fontsize=6,
                 style="italic", color=GREY)
    axb.text(fx(X_DEP), yy, name, transform=trb, ha="left", va="center", fontsize=6.5, color=TEXT)
    axb.text(fx(X_N), yy, f"{int(rp.n):,} / {int(rp.events):,}", transform=trb, ha="right", va="center",
             fontsize=6.5, color=TEXT)
    for xmm, r in [(X_HR1, rp), (X_HR2, rb)]:
        axb.text(fx(xmm), yy, f"{r.hr:.2f} ({r.lower:.2f}–{r.upper:.2f})", transform=trb, ha="left", va="center",
                 fontsize=6.5, color=TEXT)
hy = -1.25
for xmm, lab, ha in [(X_DIS, "Outcome", "left"), (X_DEP, "Fat depot", "left"), (X_N, "n / events", "right"),
                     (X_HR1, "Primary model", "left"), (X_HR2, "Adjusted for BMI", "left")]:
    axb.text(fx(xmm), hy, lab, transform=trb, ha=ha, va="center", fontsize=6.5, fontweight="bold", color=TEXT)
fig.text(1.5 / FW, 1 - B_HEAD / FH, "b", fontsize=8, fontweight="bold", va="bottom", ha="left")
fig.text(6.0 / FW, 1 - B_HEAD / FH, "Fat depots and disease risk after imaging", fontsize=7, fontweight="bold",
         va="bottom", ha="left")
# marker keys in the HR column headers (filled = primary model, open = adjusted for BMI)
for xmm, face in [(X_HR1, "#4D4D4D"), (X_HR2, "white")]:
    axb.plot([fx(xmm - 2.0)], [hy], marker="s", ms=3.4, mfc=face, mec="#4D4D4D", mew=0.8, transform=trb,
             clip_on=False, ls="none")

save(fig, "Fig4_mri_mechanism")
from PIL import Image

for ext, kw in [("tiff", {"compression": "tiff_lzw", "dpi": (600, 600)}), ("png", {"dpi": (300, 300)})]:
    p = FIGDIR / f"Fig4_mri_mechanism.{ext}"
    with Image.open(p) as im:
        im.load()
        rgba = im.convert("RGBA")
        rgb = Image.new("RGB", im.size, (255, 255, 255))
        rgb.paste(rgba, mask=rgba.split()[-1])
    rgb.save(p, **kw)
print(f"saved; panel b rows end at {B_TOP + nrow * ROW:.1f} of {FH} mm")
