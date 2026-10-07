"""Shared code for Fig. 1 (fig1_design.py) and Fig. S1 (figS1_dag.py): cohort counts recomputed from the analysis
datasets, and drawing helpers.

build_data() recomputes the compositional means of the Table 1 population (the union of the five disease cohorts) and
the participant flow from data/accel_base_flow.csv, data/flow_<outcome>.csv and data/analysis_<outcome>.csv, checks
them against the published numbers and writes figures/data/fig1_composition.csv and fig1_flow.csv.
"""
from style import *  # noqa: F401,F403

import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

DATA = PROJ / "data"
OUTCOMES5 = ["gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs"]
PARTS = ["sleep", "inactivity", "light", "mvpa"]


# ----------------------------------------------------------------------------- data
def build_data():
    cols = ["eid", "sleep", "inactivity_thr", "light_thr", "mvpa_thr"]
    files = {o: pd.read_csv(DATA / f"analysis_{o}.csv", usecols=cols + ["gs"]) for o in
             OUTCOMES5 + ["gallstone_surgery", "nco_cyst", "nco_earwax"]}
    union = pd.concat([files[o][cols] for o in OUTCOMES5]).drop_duplicates("eid")
    x = union[cols[1:]].to_numpy()
    assert np.allclose(x.sum(1), 24)
    grp = np.where(union["sleep"] < 7, "short", np.where(union["sleep"] <= 9, "normal", "long"))
    rows = []
    for g, m in [("all", np.ones(len(union), bool)), ("short", grp == "short"),
                 ("normal", grp == "normal"), ("long", grp == "long")]:
        gm = np.exp(np.log(x[m]).mean(0))
        comp = 24 * gm / gm.sum()
        rows.append(dict(group=g, n=int(m.sum()), **{p: c for p, c in zip(PARTS, comp)}))
    comp = pd.DataFrame(rows)
    # check against the published compositional means (Table 1, h/day)
    expect = {"all": (91712, 8.95, 9.09, 4.28, 1.68), "short": (3531, 6.66, 10.66, 4.86, 1.82),
              "normal": (51074, 8.31, 9.40, 4.49, 1.80), "long": (37107, 10.11, 8.47, 3.92, 1.50)}
    for _, r in comp.iterrows():
        e = expect[r.group]
        assert r.n == e[0], (r.group, r.n)
        assert np.allclose([round(r[p], 2) for p in PARTS], e[1:]), (r.group, [r[p] for p in PARTS])
    comp.to_csv(DATADIR / "fig1_composition.csv", index=False)

    flow = pd.read_csv(DATA / "accel_base_flow.csv", index_col=0)["n"].tolist()
    assert flow == [101393, 92276, 92266, 91721]
    cohorts = {o: (len(files[o]), int(files[o]["gs"].sum())) for o in files}
    expect_c = {"gallstone": (86596, 1961), "peptic_ulcer": (89540, 983), "gerd": (82983, 5577),
                "diverticular": (86008, 7167), "ibs": (83646, 1031), "gallstone_surgery": (86596, 1042),
                "nco_cyst": (89481, 646), "nco_earwax": (84632, 577)}
    assert cohorts == expect_c, cohorts
    # outcome-specific exclusions (row d of the flow files)
    excl = {}
    for o in ["gallstone", "ibs", "nco_earwax"]:
        f = pd.read_csv(DATA / f"flow_{o}.csv", index_col=0)["n"]
        excl[o] = int(f.iloc[2] - f.iloc[3])
    assert excl == {"gallstone": 1135, "ibs": 2575, "nco_earwax": 2311}, excl
    fl = pd.DataFrame([dict(item=k, n=v[0], events=v[1]) for k, v in cohorts.items()] +
                      [dict(item=s, n=n, events=np.nan) for s, n in
                       zip(["returned", "device_qc", "complete", "exposure_base"], flow)] +
                      [dict(item=f"excluded_{k}", n=v, events=np.nan) for k, v in excl.items()] +
                      [dict(item="union_table1", n=len(union), events=np.nan)])
    fl.to_csv(DATADIR / "fig1_flow.csv", index=False)
    return comp.set_index("group"), flow, cohorts, excl, len(union)


# ----------------------------------------------------------------------------- drawing helpers
FIG_H = 226.8  # mm
TXT = "#222222"
EDGE = "#6E6E6E"
LIGHT = "#F2F2F2"
# neutral structural colour (not a behaviour colour): timeline track, exposure base, effect measures, MRI arm
STRUCT = "#5B6779"
# three type sizes: body 5.5 pt (5.0 pt only for footnotes and interval/exclusion numbers),
# bold headers 6.5 pt, panel letters 8 pt (style.panel_label)
FS, FS_SMALL, FS_HEAD = 5.5, 5.0, 6.5


def fmt(n):
    return f"{n:,}"


def T(ax, x, y, s, size=FS, weight="normal", color=TXT, ha="center", va="center_baseline", **kw):
    return ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va,
                   linespacing=1.15, **kw)


def box(ax, x, y, w, h, fc="white", ec=EDGE, lw=0.6, r=0.8, ls="-", z=1):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw,
                       ls=ls, zorder=z)
    ax.add_patch(p)
    return p


def arrow(ax, x0, y0, x1, y1, color=EDGE, lw=0.7, style="-|>", ms=5, z=2, ls="-"):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=ms, color=color, lw=lw,
                        shrinkA=0, shrinkB=0, zorder=z, linestyle=ls)
    ax.add_patch(a)
    return a


def line(ax, xs, ys, color=EDGE, lw=0.6, z=2, ls="-"):
    ax.plot(xs, ys, color=color, lw=lw, zorder=z, ls=ls, solid_capstyle="butt")


def letter(ax, x_mm, y_mm, s):
    panel_label(ax, s, x=x_mm / W_MM, y=1 - y_mm / FIG_H)


def tint(hex_color, a):
    """Blend a colour with white (a = weight of the colour)."""
    h = hex_color.lstrip("#")
    rgb = np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)]) / 255
    return tuple(1 - a * (1 - rgb))


W_MM = 183.0
