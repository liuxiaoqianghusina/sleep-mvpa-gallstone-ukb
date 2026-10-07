"""Supplementary Fig. S3 | Exploratory search for the lowest-risk and highest-risk 24-hour compositions.

Method (scripts/14_composition_search.R): synthetic compositions in 15-min steps within the 2.5th-97.5th percentile
of each behaviour, implausible ones (KDE log-density below the 5th percentile of the observed) removed. The
lowest-risk and highest-risk compositions were chosen in a random half of the cohort (training half); the typical
composition is the one with the highest KDE density. All three were then evaluated in the other half (testing
half) with the g-formula refitted there, 200 bootstrap replicates.

a  The three compositions for each disease (h/day), 9.8-year risk per 1,000 in the testing half (95% CI) and RR vs
   the typical composition (95% CI).
b  Cumulative incidence over follow-up in the testing half under each composition (95% CI bands).

Data: results/ideal_composition (gallstone disease) and results/outcomes/<outcome>/ideal_composition:
reference_compositions.csv, cumulative_incidence.csv, contrasts_vs_typical_evalues.csv.
"""
from style import *  # noqa: F403

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch
from matplotlib.ticker import FuncFormatter

COMPS = ["Lowest risk", "Typical", "Highest risk"]
CCOL = {"Lowest risk": BEHAVIOUR_COLOR["mvpa"], "Typical": "#7A7A7A", "Highest risk": "#B2182B"}
BEH = [("sleep", "Sleep"), ("inactivity", "Inactivity"), ("light", "Light activity"), ("mvpa", "MVPA")]
SEG = {"sleep": BEHAVIOUR_COLOR["sleep"], "inactivity": "#CFCECE", "light": BEHAVIOUR_COLOR["light"],
       "mvpa": BEHAVIOUR_COLOR["mvpa"]}
TXT = "#1A1A1A"
MINUS = "−"


def rdir(o):
    return PROJ / "results" / ("ideal_composition" if o == "gallstone" else f"outcomes/{o}/ideal_composition")


def load(o):
    d = rdir(o)
    ref = pd.read_csv(d / "reference_compositions.csv").set_index("composition").loc[COMPS]
    assert np.allclose(ref[[b for b, _ in BEH]].sum(axis=1), 24), o
    ci = pd.read_csv(d / "cumulative_incidence.csv")
    con = pd.read_csv(d / "contrasts_vs_typical_evalues.csv").set_index("contrast")
    end = ci[ci.timegroup == ci.timegroup.max()].set_index("composition").loc[COMPS]
    assert abs(end.years.iloc[0] - 9.81) < 0.01, o
    return ref, ci, con, end


DATA = {o: load(o) for o in DISEASES}
# printed check against the numbers quoted in the text
for o, (ref, ci, con, end) in DATA.items():
    for c in COMPS:
        r = end.loc[c]
        rr = "1 (reference)" if c == "Typical" else \
            "{:.2f} ({:.2f}-{:.2f})".format(*con.loc[c.split()[0] + " vs typical", ["rr", "lower", "upper"]])
        print(f"{o:13s} {c:13s} {'/'.join(f'{v:g}' for v in ref.loc[c, [b for b, _ in BEH]])}  "
              f"risk/1000 {1000 * r.risk:.1f} ({1000 * r.lower:.1f}-{1000 * r.upper:.1f})  RR {rr}")

FW, FH = 183.0, 166.0
fig = plt.figure(figsize=(FW * MM, FH * MM))


def ax_mm(left, top, w, h):
    return fig.add_axes([left / FW, 1 - (top + h) / FH, w / FW, h / FH])


# ----------------------------------------------------------------------------- a: compositions and risks
A_TOP, ROW, GAP = 13.0, 4.3, 1.6
X_LAB, X_BAR0, X_BAR1 = 2.0, 33.0, 95.0
X_H, X_RISK, X_RR = 99.0, 154.0, 181.0          # hours (left), risk (right-aligned), RR (right-aligned)
rows, y = [], 0.0
for o in DISEASES:
    rows.append(("head", o, None, y))
    y += 1
    for c in COMPS:
        rows.append(("row", o, c, y))
        y += 1
    y += GAP / ROW
A_H = y * ROW
axa = ax_mm(0, A_TOP, FW, A_H)
axa.set_xlim(0, FW)
axa.set_ylim(y - 0.2, -0.6)
axa.axis("off")
mm_per_h = (X_BAR1 - X_BAR0) / 24
for kind, o, c, yy in rows:
    ref, ci, con, end = DATA[o]
    if kind == "head":
        axa.text(X_LAB, yy, DISEASE_LABEL[o] + (" (primary outcome)" if o == "gallstone" else ""), ha="left",
                 va="center", fontsize=6.5, fontweight="bold", color=TXT)
        continue
    axa.text(X_LAB + 2.5, yy, c, ha="left", va="center", fontsize=6.2, color=CCOL[c],
             fontweight="bold" if c != "Typical" else "normal")
    x = X_BAR0
    for b, _ in BEH:
        w = ref.loc[c, b] * mm_per_h
        axa.add_patch(Rectangle((x, yy - 0.36), w, 0.72, fc=SEG[b], ec="white", lw=0.4))
        x += w
    hrs = " / ".join(f"{ref.loc[c, b]:.2f}" for b, _ in BEH)
    axa.text(X_H, yy, hrs, ha="left", va="center", fontsize=6.2, color=TXT)
    r = end.loc[c]
    axa.text(X_RISK, yy, f"{1000 * r.risk:.1f} ({1000 * r.lower:.1f}–{1000 * r.upper:.1f})", ha="right", va="center",
             fontsize=6.2, color=TXT)
    if c == "Typical":
        s = "1 (reference)"
    else:
        q = con.loc[c.split()[0] + " vs typical"]
        s = f"{q.rr:.2f} ({q.lower:.2f}–{q.upper:.2f})"
    axa.text(X_RR, yy, s, ha="right", va="center", fontsize=6.2, color=TXT,
             fontweight="bold" if c != "Typical" and (q.lower > 1 or q.upper < 1) else "normal")
# hour scale under the bars
yb = y - 0.55
axa.plot([X_BAR0, X_BAR1], [yb, yb], color="#4D4D4D", lw=0.5)
for h in [0, 6, 12, 18, 24]:
    xx = X_BAR0 + h * mm_per_h
    axa.plot([xx, xx], [yb, yb + 0.25], color="#4D4D4D", lw=0.5)
    axa.text(xx, yb + 0.85, str(h), ha="center", va="center", fontsize=6, color=TXT)
# column headers
hy = 1 - (A_TOP - 2.6) / FH
fig.text((X_BAR0 + X_BAR1) / 2 / FW, hy, "24-hour composition (h/day)", ha="center", va="bottom", fontsize=6.5,
         fontweight="bold")
fig.text(X_H / FW, hy + 2.8 / FH, "Sleep / inactivity /", ha="left", va="bottom", fontsize=6.5, fontweight="bold")
fig.text(X_H / FW, hy, "light / MVPA (h/day)", ha="left", va="bottom", fontsize=6.5, fontweight="bold")
fig.text(X_RISK / FW, hy + 2.8 / FH, "9.8-year risk", ha="right", va="bottom", fontsize=6.5, fontweight="bold")
fig.text(X_RISK / FW, hy, "per 1,000 (95% CI)", ha="right", va="bottom", fontsize=6.5, fontweight="bold")
fig.text(X_RR / FW, hy + 2.8 / FH, "RR vs typical", ha="right", va="bottom", fontsize=6.5, fontweight="bold")
fig.text(X_RR / FW, hy, "(95% CI)", ha="right", va="bottom", fontsize=6.5, fontweight="bold")
fig.legend(handles=[Patch(fc=SEG[b], ec="white", label=lab) for b, lab in BEH], loc="upper left",
           bbox_to_anchor=(X_BAR0 / FW - 0.005, 1 - (A_TOP + A_H + 3.0) / FH), ncol=4, fontsize=6, frameon=False,
           handlelength=1.2, handleheight=0.9, columnspacing=1.2, handletextpad=0.4)
fig.text(1.0 / FW, 1 - 2.0 / FH, "a", fontsize=8, fontweight="bold", va="top", ha="left")

# ----------------------------------------------------------------------------- b: cumulative incidence curves
B_TOP = A_TOP + A_H + 17.0
B_H, B_GAP = 30.0, 8.5
B_L, B_R = 14.0, 181.0
bw = (B_R - B_L - 4 * B_GAP) / 5
for k, o in enumerate(DISEASES):
    ax = ax_mm(B_L + k * (bw + B_GAP), B_TOP, bw, B_H)
    ref, ci, con, end = DATA[o]
    for c in COMPS:
        d = ci[ci.composition == c].sort_values("timegroup")
        xs = np.r_[0, d.years.values]
        ax.fill_between(xs, 1000 * np.r_[0, d.lower.values], 1000 * np.r_[0, d.upper.values], color=CCOL[c],
                        alpha=0.15, lw=0)
        ax.plot(xs, 1000 * np.r_[0, d.risk.values], color=CCOL[c], lw=1.1 if c != "Typical" else 0.9)
    ax.set_xlim(0, 10)
    ax.set_xticks([0, 5, 10])
    ax.set_ylim(0, None)
    ax.set_title(DISEASE_LABEL[o], fontsize=6.5, pad=3)
    ax.set_xlabel("Years", fontsize=6.5, labelpad=1.5)
    if k == 0:
        ax.set_ylabel("Cumulative incidence\nper 1,000 (95% CI)", fontsize=6.5, labelpad=2)
fig.legend(handles=[Line2D([], [], color=CCOL[c], lw=1.1, label=f"{c} composition") for c in COMPS],
           loc="upper center", bbox_to_anchor=(0.5, 1 - (B_TOP + B_H + 8.0) / FH), ncol=3, fontsize=6,
           frameon=False, handlelength=2.0, columnspacing=1.8)
fig.text(1.0 / FW, 1 - (B_TOP - 6.5) / FH, "b", fontsize=8, fontweight="bold", va="top", ha="left")

plt.rcParams["savefig.pad_inches"] = 0.02
save(fig, "FigS3_compositions")
print(f"saved; b ends near {B_TOP + B_H + 12:.0f} of {FH} mm")
