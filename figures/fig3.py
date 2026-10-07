"""Fig. 3 | Reallocating 30 min/day from sleep to MVPA: sensitivity analyses, disease-specific checks
and negative control outcomes.

Forest plot of the full-cohort RR (95% percentile bootstrap CI) for sleep -> MVPA, 30 min/day.
Rows are grouped by disease (style.DISEASES order): primary analysis, the three core sensitivity
analyses (lag3, s2, healthy-user adjustment), then the disease-specific checks; the last block holds
the two negative control outcomes. Right-hand columns give RR (95% CI) and n / events.

Input : figures/data/fig3_robustness.csv (written by figures/data/fig3_data.R)
Output: Fig3_robustness.{svg,pdf,tiff,png} and Fig3_robustness_preview.png via style.save()
"""
from style import *  # noqa: F401,F403

import pandas as pd
from matplotlib import transforms
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

df = pd.read_csv(DATADIR / "fig3_robustness.csv").sort_values("order")

# ---- labels (English terms from the manuscript glossary, block E and D10) -------------------------
GROUP_TITLE = {g: DISEASE_LABEL[g] for g in DISEASES}
GROUP_TITLE["gallstone"] = "Gallstone disease (primary outcome)"
GROUP_TITLE["nco"] = "Negative control outcomes"
ROW_LABEL = {
    "primary": "Primary analysis",
    "lag3": "Excluding cases diagnosed in the first 3 years (lag3)",
    "s2": "+ BMI and comorbidities (s2)",
    "healthuser": "Healthy-user adjustment",
    "cholecystectomy": "Cholecystectomy endpoint",
    "s3": "+ gallstone-specific factors (s3)",
    "ml": "Machine-learning exposure",
    "knots4": "Four-knot splines",
    "reduced": "Main-effects model",
    "mrisub": "MRI subsample",
    "noacid": "Excluding acid-suppressant and antacid users",
    "nsaid": "+ aspirin and NSAID use",
    "nco_cyst": "Follicular cyst (L72)",
    "nco_earwax": "External ear disorders (H61)",
}
GROUP_ORDER = DISEASES + ["nco"]

# ---- vertical layout (data units: one row = 1) ------------------------------------------------------
GAP_SPECIFIC = 0.45   # whitespace between core sensitivity analyses and disease-specific checks
GAP_GROUP = 0.6       # whitespace between disease blocks
rows, headers, blocks, sections = [], [], [], []
y = 0.0
for g in GROUP_ORDER:
    sub = df[df["group"] == g]
    if g == "peptic_ulcer":  # section header for the four secondary outcomes
        sections.append(("Secondary outcomes", y))
        y += 1
    y_top = y
    headers.append((g, y))
    y += 1
    first_specific = True
    for _, r in sub.iterrows():
        if r["tier"] == "specific" and first_specific:
            y += GAP_SPECIFIC
            first_specific = False
        rows.append((y, r))
        y += 1
    blocks.append((g, y_top, y - 1))
    y += GAP_GROUP
y_last = rows[-1][0]

# ---- figure geometry ---------------------------------------------------------------------------------
ROW_MM = 3.45
TOP_MM, BOTTOM_MM = 9.0, 11.0
span = (y_last + 0.6) - (-0.6)
H = (span * ROW_MM + TOP_MM + BOTTOM_MM) * MM
# style.save() uses bbox_inches="tight" (0.1-in pad each side); content spans 0..1 of the figure,
# so a canvas 0.2 in narrower than W2 exports at exactly 183 mm wide
FIG_W = W2 - 0.2
fig = plt.figure(figsize=(FIG_W, H))
AX_L, AX_W = 0.318, 0.372
ax = fig.add_axes([AX_L, BOTTOM_MM * MM / H, AX_W, span * ROW_MM * MM / H])
ax.set_ylim(y_last + 0.6, -0.6)
ax.patch.set_visible(False)

X_HEAD, X_LAB = 0.004, 0.018           # disease header / indented analysis label (figure fraction)
X_RR = 0.712                           # RR (95% CI) text, left aligned
X_N, X_SLASH, X_EV = 0.912, 0.9205, 0.929  # n right-aligned / slash / events left-aligned
trans = transforms.blended_transform_factory(fig.transFigure, ax.transData)
FS = 6.5

# pale bands behind alternate disease blocks (span the text columns too, so rows read across)
import matplotlib.colors as mcolors
_gc = mcolors.to_rgb(DISEASE_COLOR["gallstone"])
PRIMARY_BAND = tuple(c + (1 - c) * 0.88 for c in _gc)
for i, (g, y0, y1) in enumerate(blocks):
    if g == "gallstone":
        ax.add_patch(Rectangle((0.0, y0 - 0.5), 1.0, (y1 - y0) + 1.0, transform=trans,
                               facecolor=PRIMARY_BAND, edgecolor="none", zorder=0, clip_on=False))
    elif i % 2 == 0 and g != "nco":
        ax.add_patch(Rectangle((0.0, y0 - 0.5), 1.0, (y1 - y0) + 1.0, transform=trans,
                               facecolor="#F2F2F2", edgecolor="none", zorder=0, clip_on=False))
for title, yy in sections:
    ax.text(X_HEAD, yy, title, transform=trans, ha="left", va="center", fontsize=6.5, style="italic",
            color=GREY, clip_on=False)

# reference line at RR = 1
ax.axvline(1.0, color="#8C8C8C", lw=0.6, ls=(0, (3, 2)), zorder=1)

# group headers
for g, yy in headers:
    ax.text(X_HEAD, yy, GROUP_TITLE[g], transform=trans, ha="left", va="center",
            fontsize=7, fontweight="bold", clip_on=False)

# estimate rows: one marker size for every analysis (the legend describes no size encoding);
# negative-control rows get a slightly heavier CI line, and the pale H61 grey (style.py) a darker
# marker edge, so they stay legible in print and grayscale
for yy, r in rows:
    key = r["key"]
    is_nco = r["group"] == "nco"
    colour = DISEASE_COLOR[key] if is_nco else DISEASE_COLOR[r["group"]]
    ax.plot([r["lower"], r["upper"]], [yy, yy], color=colour, lw=1.1 if is_nco else 0.9,
            solid_capstyle="butt", zorder=3)
    edge = "#7F7F7F" if key == "nco_earwax" else colour
    ax.plot(r["rr"], yy, marker="s", ms=3.6, color=colour, mec=edge,
            mew=0.4 if key == "nco_earwax" else 0, zorder=4)
    ax.text(X_LAB, yy, ROW_LABEL[key], transform=trans, ha="left", va="center", fontsize=FS,
            clip_on=False)
    ax.text(X_RR, yy, f"{r['rr']:.2f} ({r['lower']:.2f}–{r['upper']:.2f})", transform=trans,
            ha="left", va="center", fontsize=FS, clip_on=False)
    ax.text(X_N, yy, f"{int(r['N']):,}", transform=trans, ha="right", va="center", fontsize=FS,
            clip_on=False)
    ax.text(X_SLASH, yy, "/", transform=trans, ha="center", va="center", fontsize=FS, clip_on=False)
    ax.text(X_EV, yy, f"{int(r['events']):,}", transform=trans, ha="left", va="center",
            fontsize=FS, clip_on=False)

# column headers, placed just above the first block
fig_y_head = (BOTTOM_MM + span * ROW_MM + 4.2) * MM / H
head_kw = dict(transform=fig.transFigure, va="center", fontsize=7, fontweight="bold")
fig.text(X_HEAD, fig_y_head, "Disease and analysis", ha="left", **head_kw)
fig.text(AX_L + AX_W / 2, fig_y_head, "Sleep→MVPA, 30 min/day", ha="center", **head_kw)
fig.text(X_RR, fig_y_head, "RR (95% CI)", ha="left", **head_kw)
fig.text(X_SLASH, fig_y_head, "n / events", ha="center", **head_kw)
fig_y_rule = (BOTTOM_MM + span * ROW_MM + 1.3) * MM / H
fig.add_artist(plt.Line2D([0.0, 1.0], [fig_y_rule, fig_y_rule], transform=fig.transFigure,
                          color="#4D4D4D", lw=0.5))

# x axis: shared log scale
ax.set_xscale("log")
ax.set_xlim(0.64, 1.12)
ticks = [0.7, 0.8, 0.9, 1.0, 1.1]
ax.xaxis.set_major_locator(FixedLocator(ticks))
ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
ax.xaxis.set_minor_locator(NullLocator())
ax.set_xlabel("RR (95% CI)", labelpad=2)
ax.spines["left"].set_visible(False)
ax.set_yticks([])

save(fig, "Fig3_robustness")

# Flatten the raster exports to RGB on white (production expects RGB TIFF); style.save() is shared and left as is.
from PIL import Image

for ext, kw in [("tiff", {"compression": "tiff_lzw", "dpi": (600, 600)}), ("png", {"dpi": (300, 300)})]:
    p = FIGDIR / f"Fig3_robustness.{ext}"
    with Image.open(p) as im:
        im.load()
        rgba = im.convert("RGBA")
        rgb = Image.new("RGB", im.size, (255, 255, 255))
        rgb.paste(rgba, mask=rgba.split()[-1])
    rgb.save(p, **kw)
print(f"canvas: {FIG_W / MM:.1f} x {H / MM:.1f} mm (+5.1 mm tight-bbox pad); rows = {len(rows)}")
