"""Supplementary Fig. S2 | Natural-course check of the parametric g-formula.

For each outcome, cumulative incidence predicted by the g-formula under the observed compositions (natural course;
line) is compared with the nonparametric Aalen-Johansen estimate with death as a competing event (points with 95%
CI) at the end of each of the 10 follow-up intervals. Data: results/natural_course/<outcome>.csv
(scripts/15_natural_course.R).
"""
from style import *  # noqa: F403

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

PANELS = [("gallstone", "Gallstone disease", DISEASE_COLOR["gallstone"]),
          ("gallstone_surgery", "Cholecystectomy endpoint", DISEASE_COLOR["gallstone"]),
          ("peptic_ulcer", "Peptic ulcer", DISEASE_COLOR["peptic_ulcer"]),
          ("gerd", "GERD", DISEASE_COLOR["gerd"]),
          ("diverticular", "Diverticular disease", DISEASE_COLOR["diverticular"]),
          ("ibs", "IBS", DISEASE_COLOR["ibs"]),
          ("nco_cyst", "Follicular cyst (L72)", "#6E6E6E"),
          ("nco_earwax", "External ear disorders (H61)", "#6E6E6E")]
OBS = "#1A1A1A"

FW, FH = 183.0, 92.0
fig = plt.figure(figsize=(FW * MM, FH * MM))
L, R, TOP, H, GX, GY = 14.0, 181.0, 7.0, 29.0, 9.5, 14.0
W = (R - L - 3 * GX) / 4


def ax_mm(left, top, w, h):
    return fig.add_axes([left / FW, 1 - (top + h) / FH, w / FW, h / FH])


summary = []
for k, (o, title, col) in enumerate(PANELS):
    d = pd.read_csv(PROJ / "results/natural_course" / f"{o}.csv").sort_values("timegroup")
    r, c = divmod(k, 4)
    ax = ax_mm(L + c * (W + GX), TOP + r * (H + GY), W, H)
    x = np.r_[0, d.years]
    ax.plot(x, 1000 * np.r_[0, d.pred_risk], color=col, lw=1.4, zorder=2)
    ax.errorbar(d.years, 1000 * d.obs_risk, yerr=[1000 * (d.obs_risk - d.obs_lower), 1000 * (d.obs_upper - d.obs_risk)],
                fmt="o", ms=2.6, color=OBS, mfc="white", mew=0.7, elinewidth=0.6, capsize=0, zorder=3)
    ax.set_xlim(0, 10.3)
    ax.set_xticks([0, 2, 4, 6, 8, 10])
    ax.set_ylim(0, None)
    ax.set_title(title, fontsize=6.5, pad=3)
    if r == 1:
        ax.set_xlabel("Years since day 0")
    if c == 0:
        ax.set_ylabel("Cumulative incidence\nper 1,000", labelpad=2)
    last = d.iloc[-1]
    ax.text(0.03, 0.97, f"9.8 y: predicted {1000 * last.pred_risk:.1f}\nobserved {1000 * last.obs_risk:.1f} "
            f"({1000 * last.obs_lower:.1f}–{1000 * last.obs_upper:.1f})", transform=ax.transAxes, ha="left", va="top",
            fontsize=5.3, color=GREY, linespacing=1.2)
    summary.append((o, 1000 * last.pred_risk, 1000 * last.obs_risk, 1000 * last.obs_lower, 1000 * last.obs_upper,
                    1000 * (d.pred_risk - d.obs_risk).abs().max(),
                    int(((d.pred_risk >= d.obs_lower) & (d.pred_risk <= d.obs_upper)).sum())))
fig.legend(handles=[Line2D([], [], color="#4D4D4D", lw=1.4, label="Predicted by the g-formula (natural course)"),
                    Line2D([], [], color=OBS, marker="o", ms=2.6, mfc="white", mew=0.7, ls="none",
                           label="Observed (Aalen–Johansen, 95% CI)")],
           loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=6, frameon=False, handlelength=2.2,
           columnspacing=2.0)
save(fig, "FigS2_natural_course")
for s in summary:
    print("%-18s pred %.1f  obs %.1f (%.1f-%.1f)  max|diff| %.2f  pred within obs CI at %d/10 intervals" % s)
