"""Figure 1 | Study design and participant flow (conventional two-panel layout).

a  Study design diagram in the format of Schneeweiss et al. (Ann Intern Med 2019;170:398-406):
   day 0 (time zero), exclusion, covariate and exposure assessment windows, follow-up and the abdominal
   MRI substudy, on a time axis in years relative to day 0.
b  STROBE participant flow with the reason and number for every exclusion, branching into the primary
   outcome cohort, the secondary outcome cohorts and the negative control outcome cohorts.
No cross-references to other figures or tables appear inside the figure.

Sources (all checked at run time where possible):
  flow counts         data/accel_base_flow.csv, data/flow_<outcome>.csv and data/analysis_<outcome>.csv
                      (fig1_common.build_data())
  device-QC reasons   sequential breakdown of the combined QC filter in scripts/01_accelerometry.py
                      (p90015 wear time, p90016/p90017 calibration, p90012 >= 100 mg, p90002 data problem)
  timing              recruitment to day 0 median 5.7 y (IQR 4.9-6.5); MRI median 3.68 y (IQR 1.91-4.99)
                      after day 0 (MRI extracted from the gallstone cohort, scripts/07_extract_mri.R)
  "at or before recruitment" includes one participant without linked outcome records (absent from every
  outcome file).
"""
from style import *  # noqa: F401,F403
from fig1_common import build_data, T, box, arrow, line, tint, fmt, EDGE, LIGHT, TXT, FS, FS_SMALL, FS_HEAD

import pandas as pd
from matplotlib.patches import Rectangle, Polygon

W_MM, FIG_H = 183.0, 169.5
GS = DISEASE_COLOR["gallstone"]
WIN = "#C9CED6"      # assessment windows
WIN_DARK = "#5B6779"
GREY_TXT = "#5A5A5A"
EX_EDGE = "#9A9A9A"

# device-QC exclusions, applied sequentially in this order (sum = 9,117)
QC = [("Insufficient wear time", 6820), ("Calibration problems", 4 + 159),
      ("Data problem or mean acceleration ≥100 mg", 2127 + 7)]
MRI_N = 19752


def letter(ax, x, y, s):
    ax.text(x, y, s, fontsize=8, fontweight="bold", ha="left", va="top")


# ---------------------------------------------------------------------------------------- panel a
def panel_a(ax):
    letter(ax, 0.0, 0.3, "a")
    t_lo, t_hi = -7.6, 10.4
    x_lo, x_hi = 60.0, 179.0
    X = lambda t: x_lo + (t - t_lo) * (x_hi - x_lo) / (t_hi - t_lo)  # noqa: E731
    x0 = X(0)
    dy, y1, bh = 7.6, 10.6, 4.4
    rows = [("Exclusion assessment", "Outcome or outcome-specific criteria"),
            ("Covariate assessment", "Baseline visit, 2006–2010; age at day 0"),
            ("Exposure assessment", "24-h composition from wrist accelerometry"),
            ("Follow-up", "To outcome, death, loss or end of data"),
            ("Abdominal MRI substudy", f"n = {fmt(MRI_N)} from the gallstone disease cohort")]
    ys = [y1 + i * dy for i in range(len(rows))]
    for (lab, desc), y in zip(rows, ys):
        T(ax, 3.0, y - 1.3, lab, ha="left", weight="bold")
        T(ax, 3.0, y + 1.6, desc, ha="left", size=FS_SMALL, color=GREY_TXT)
    # day 0
    line(ax, [x0, x0], [5.0, ys[-1] + bh / 2 + 1.8], color="#222222", lw=0.9)
    T(ax, x0, 3.0, "Day 0 (time zero): end of 7-day accelerometer wear, 2013–2015", weight="bold")

    def bar(y, a, b, fc=WIN, z=2):
        ax.add_patch(Rectangle((a, y - bh / 2), b - a, bh, fc=fc, ec="none", zorder=z))

    # exclusion window: all records before day 0
    y = ys[0]
    bar(y, x_lo, x0)
    ax.add_patch(Polygon([[x_lo - 2.2, y], [x_lo, y - bh / 2], [x_lo, y + bh / 2]], closed=True, fc=WIN, ec="none"))
    T(ax, (x_lo + x0) / 2, y, "Any record before day 0", size=FS_SMALL)
    # covariate window: baseline visit, median 5.7 (IQR 4.9-6.5) years before day 0
    y = ys[1]
    bar(y, X(-6.5), X(-4.9))
    ax.plot([X(-5.7)], [y], "o", ms=2.8, color=WIN_DARK, zorder=3)
    T(ax, X(-4.9) + 1.5, y, "Median 5.7 years before day 0", ha="left", size=FS_SMALL)
    # exposure window: days -7 to 0 (widened for visibility)
    y = ys[2]
    bar(y, x0 - 1.3, x0, fc=WIN_DARK)
    T(ax, x0 - 2.3, y, "Days −7 to 0", ha="right", size=FS_SMALL)
    # follow-up
    y = ys[3]
    bar(y, x0, X(9.81))
    ax.add_patch(Polygon([[X(9.81), y - bh / 2], [X(9.81) + 2.2, y], [X(9.81), y + bh / 2]], closed=True, fc=WIN,
                         ec="none"))
    T(ax, (x0 + X(9.81)) / 2, y, "Up to 9.8 years; end of data 2022–2023", size=FS_SMALL)
    # MRI: median 3.68 (IQR 1.91-4.99) years after day 0
    y = ys[4]
    bar(y, X(1.91), X(4.99))
    ax.plot([X(3.68)], [y], "o", ms=2.8, color=WIN_DARK, zorder=3)
    T(ax, X(4.99) + 1.5, y, "Median 3.7 years after day 0", ha="left", size=FS_SMALL)
    # time axis
    ya = ys[-1] + bh / 2 + 3.0
    line(ax, [x_lo, x_hi], [ya, ya], color="#4D4D4D", lw=0.6)
    for t in range(-6, 11, 2):
        line(ax, [X(t), X(t)], [ya, ya + 0.9], color="#4D4D4D", lw=0.6)
        T(ax, X(t), ya + 2.6, f"{t}".replace("-", "−"), size=FS_SMALL)
    T(ax, (x_lo + x_hi) / 2, ya + 5.5, "Years relative to day 0", size=FS_SMALL)
    # key for dots and bars
    kx, ky = 3.0, ya + 1.6
    ax.add_patch(Rectangle((kx, ky - 1.3), 6.0, 2.6, fc=WIN, ec="none"))
    ax.plot([kx + 3.0], [ky], "o", ms=2.8, color=WIN_DARK)
    T(ax, kx + 7.5, ky, "Median (dot) and interquartile range (bar)", ha="left", size=FS_SMALL, color=GREY_TXT)
    return ya + 6.5


# ---------------------------------------------------------------------------------------- panel b
def excl_box(ax, x, y, w, head, items, lh=2.7):
    h = 2.4 + lh * len(items) + 2.2
    box(ax, x, y, w, h, fc="white", ec=EX_EDGE, lw=0.5)
    T(ax, x + 2.0, y + 2.5, head, ha="left", weight="bold", size=FS_SMALL)
    for i, s in enumerate(items):
        T(ax, x + 2.0, y + 2.5 + lh * (i + 1), s, ha="left", size=FS_SMALL)
    return h


def cohort_table(ax, x, y, w, rows, base_n):
    """Mini table: outcome | excluded | analysed | cases."""
    cols = [w - 19.0, w - 9.5, w - 2.0]  # right edges of the three number columns
    rh = 4.6
    h = 4.6 + rh * len(rows) + 0.6
    box(ax, x, y, w, h, fc=LIGHT, ec=EDGE, lw=0.6)
    for c, s in zip(cols, ["Excluded", "Analysed", "Cases"]):
        T(ax, x + c, y + 2.8, s, ha="right", size=FS_SMALL, color=GREY_TXT)
    line(ax, [x + 1.5, x + w - 1.5], [y + 4.4, y + 4.4], color="#BDBDBD", lw=0.4)
    for i, (name, n, e) in enumerate(rows):
        yy = y + 4.6 + rh * (i + 0.6)
        T(ax, x + 2.0, yy, name, ha="left", weight="bold", size=FS_SMALL)
        for c, v in zip(cols, [base_n - n, n, e]):
            T(ax, x + c, yy, fmt(v), ha="right", size=FS_SMALL)
    return h


def panel_b(ax, flow, cohorts, excl, y0):
    letter(ax, 0.0, y0, "b")
    gs_flow = pd.read_csv(PROJ / "data/flow_gallstone.csv", index_col=0)["n"].tolist()
    assert gs_flow == [88924, 87734, 87731, 86596, 86596], gs_flow
    ex_gs = [flow[3] - gs_flow[0], gs_flow[0] - gs_flow[1], gs_flow[1] - gs_flow[2], gs_flow[2] - gs_flow[3]]
    assert ex_gs[3] == excl["gallstone"] and sum(n for _, n in QC) == flow[0] - flow[1]

    # main chain with exclusion boxes on the right
    cx, bw, bh = 91.5, 62.0, 9.0
    ex_x, ex_w = 129.0, 52.0
    b1 = y0 + 3.0
    e1 = b1 + bh + 1.5
    h1 = excl_box(ax, ex_x, e1, ex_w, f"Excluded, device quality control (n = {fmt(flow[0] - flow[1])})",
                  [f"{s} (n = {fmt(n)})" for s, n in QC])
    b2 = e1 + h1 + 1.5
    e2 = b2 + bh + 1.5
    h2 = excl_box(ax, ex_x, e2, ex_w, f"Excluded, 24-h composition (n = {fmt(flow[1] - flow[3])})",
                  [f"Incomplete or implausible 24-h total (n = {fmt(flow[1] - flow[2])})",
                   f"Any behaviour <0.1th or >99.9th percentile (n = {fmt(flow[2] - flow[3])})"])
    b3 = e2 + h2 + 1.5
    chain = [(b1, "UK Biobank participants with accelerometer data", flow[0]),
             (b2, "Valid accelerometer data", flow[1]),
             (b3, "Valid 24-h composition", flow[3])]
    for y, txt, n in chain:
        box(ax, cx - bw / 2, y, bw, bh, fc=LIGHT)
        T(ax, cx, y + 3.0, txt)
        T(ax, cx, y + 6.4, f"n = {fmt(n)}", weight="bold")
    for (ya, _, _), (yb_, _, _) in zip(chain[:-1], chain[1:]):
        arrow(ax, cx, ya + bh, cx, yb_ - 0.3)
    for e, h in [(e1, h1), (e2, h2)]:
        arrow(ax, cx, e + h / 2, ex_x - 0.3, e + h / 2)

    # branches
    cols = {"primary": (2.0, 57.0), "secondary": (63.0, 57.0), "nco": (124.0, 57.0)}
    centres = [x + w / 2 for x, w in cols.values()]
    yb = b3 + bh + 4.5
    line(ax, [cx, cx], [b3 + bh, yb])
    line(ax, [min(centres), max(centres)], [yb, yb])
    heads = {"primary": "Primary outcome", "secondary": "Secondary outcomes", "nco": "Negative control outcomes"}
    yh = yb + 4.2
    for k, (x, w) in cols.items():
        arrow(ax, x + w / 2, yb, x + w / 2, yh - 2.3)
        T(ax, x + w / 2, yh, heads[k], weight="bold", size=FS_HEAD, color=GS if k == "primary" else TXT)
    yt = yh + 3.0

    # primary outcome: exclusions -> cohort -> MRI substudy
    x, w = cols["primary"]
    h = excl_box(ax, x, yt, w, f"Excluded (n = {fmt(sum(ex_gs))})",
                 [f"Gallstone disease at or before recruitment (n = {fmt(ex_gs[0])})",
                  f"Gallstone disease between recruitment and day 0 (n = {fmt(ex_gs[1])})",
                  f"Died or lost to follow-up before day 0 (n = {fmt(ex_gs[2])})",
                  f"Previous cholecystectomy (n = {fmt(ex_gs[3])})"])
    cy, ch = yt + h + 3.2, 12.4
    arrow(ax, x + w / 2, yt + h, x + w / 2, cy - 0.3)
    box(ax, x, cy, w, ch, fc=tint(GS, 0.12), ec=GS, lw=0.9)
    n, e = cohorts["gallstone"]
    T(ax, x + w / 2, cy + 2.9, "Gallstone disease cohort", weight="bold")
    T(ax, x + w / 2, cy + 6.2, f"n = {fmt(n)}; {fmt(e)} incident cases")
    T(ax, x + w / 2, cy + 9.4, f"{fmt(cohorts['gallstone_surgery'][1])} treated with cholecystectomy",
      size=FS_SMALL, color=GREY_TXT)
    my, mh = cy + ch + 3.2, 7.6
    arrow(ax, x + w / 2, cy + ch, x + w / 2, my - 0.3)
    box(ax, x, my, w, mh, fc="white", ec=EDGE, lw=0.6, ls=(0, (2.2, 1.2)))
    T(ax, x + w / 2, my + 2.6, "Abdominal MRI substudy", weight="bold", size=FS_SMALL)
    T(ax, x + w / 2, my + 5.3, f"n = {fmt(MRI_N)}", size=FS_SMALL)
    y_end = my + mh

    # secondary outcomes
    base_n = flow[3]
    x, w = cols["secondary"]
    rows = [("Peptic ulcer", *cohorts["peptic_ulcer"]), ("GERD", *cohorts["gerd"]),
            ("Diverticular disease", *cohorts["diverticular"]), ("IBS", *cohorts["ibs"])]
    h = cohort_table(ax, x, yt, w, rows, base_n)
    note = "Excluded: disease, death or loss to follow-up before day 0"
    T(ax, x + 1.0, yt + h + 2.6, note, ha="left", size=FS_SMALL, color=GREY_TXT)
    T(ax, x + 1.0, yt + h + 5.2, f"For IBS, also IBD, malabsorption or colorectal cancer (n = {fmt(excl['ibs'])})",
      ha="left", size=FS_SMALL, color=GREY_TXT)
    y_end = max(y_end, yt + h + 6.5)

    # negative control outcomes
    x, w = cols["nco"]
    rows = [("Follicular cyst (L72)", *cohorts["nco_cyst"]),
            ("External ear disorders (H61)", *cohorts["nco_earwax"])]
    h = cohort_table(ax, x, yt, w, rows, base_n)
    T(ax, x + 1.0, yt + h + 2.6, note, ha="left", size=FS_SMALL, color=GREY_TXT)
    T(ax, x + 1.0, yt + h + 5.2, f"For H61, also hearing-aid users (n = {fmt(excl['nco_earwax'])})", ha="left",
      size=FS_SMALL, color=GREY_TXT)
    return max(y_end, yt + h + 6.5)


def main():
    comp, flow, cohorts, excl, n_union = build_data()
    assert excl["ibs"] == 2575 and excl["nco_earwax"] == 2311
    fig = plt.figure(figsize=(W2, FIG_H * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(FIG_H, 0)
    ax.axis("off")
    ya = panel_a(ax)
    y_end = panel_b(ax, flow, cohorts, excl, ya + 3.5)
    assert y_end <= FIG_H - 0.5, y_end
    plt.rcParams["savefig.pad_inches"] = 0
    save(fig, "Fig1_study_design")
    from PIL import Image
    for ext, kw in [("tiff", {"compression": "tiff_lzw", "dpi": (600, 600)}), ("png", {"dpi": (300, 300)})]:
        p = FIGDIR / f"Fig1_study_design.{ext}"
        with Image.open(p) as im:
            im.load()
            rgba = im.convert("RGBA")
            rgb = Image.new("RGB", im.size, (255, 255, 255))
            rgb.paste(rgba, mask=rgba.split()[-1])
        rgb.save(p, **kw)
    print(f"saved; content ends at {y_end:.1f} of {FIG_H} mm")


if __name__ == "__main__":
    main()
