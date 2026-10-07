"""Supplementary Fig. S1 | Directed acyclic graph of the assumed causal structure.

Block-level DAG: measured confounders (adjusted in all models), BMI / chronic disease / disability at recruitment
(adjusted only in the s2 sensitivity analysis, because it can be both a confounder of the current composition and a
mediator of long-term activity and sleep habits), health-seeking behaviour (healthy-user sensitivity analysis; it
can also raise the chance that disease is detected), the 24-hour composition at day 0 (exposure), visceral and
liver fat (mediators, not adjusted) and incident disease (outcome, with death as a competing event).
Covariate lists follow the Methods of the paper.
"""
from style import *  # noqa: F401,F403
from fig1_common import T, box, arrow, line, tint, EDGE

from matplotlib.patches import FancyBboxPatch

W_MM, FIG_H = 183.0, 121.0
GREY_TXT = "#5A5A5A"
FILL_ADJ = "#E4E4E4"           # adjusted in all models
FILL_MED = "#FFF3CC"           # mediator, not adjusted
EXP = BEHAVIOUR_COLOR["mvpa"]
OUT = DISEASE_COLOR["gallstone"]
S_HEAD, S_TXT = 6.5, 5.5


def node(ax, x, y, w, h, lines, fc="white", ec=EDGE, lw=0.7, ls="-", n_head=1):
    """Box with n_head bold title lines followed by smaller detail lines, centred."""
    box(ax, x, y, w, h, fc=fc, ec=ec, lw=lw, ls=ls)
    n = len(lines)
    lh = 2.75
    y0 = y + h / 2 - (n - 1) * lh / 2
    for i, s in enumerate(lines):
        head = i < n_head
        T(ax, x + w / 2, y0 + i * lh, s, size=S_HEAD if head else S_TXT, weight="bold" if head else "normal",
          color="#222222" if head else GREY_TXT)


def main():
    fig = plt.figure(figsize=(W2, FIG_H * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(FIG_H, 0)
    ax.axis("off")

    # nodes
    node(ax, 2, 3, 135, 17,
         ["Measured confounders (recruitment; adjusted in all models)",
          "Age at day 0, sex, ethnicity, Townsend deprivation index, education, household income, retirement,",
          "smoking, alcohol, healthy diet score, shift work and psychotropic medication"], fc=FILL_ADJ)
    node(ax, 14, 29, 32, 14, ["Long-term activity", "and sleep habits", "(unmeasured)"], ls=(0, (2.5, 1.5)), n_head=2)
    node(ax, 56, 29, 48, 14, ["BMI, chronic disease and", "disability at recruitment", "(sensitivity analysis)"], n_head=2)
    node(ax, 2, 55, 60, 17, ["Exposure: 24-hour composition at day 0", "sleep, inactivity, light activity and MVPA",
                             "(7-day wrist accelerometry)"], fc=tint(EXP, 0.12), ec=EXP, lw=0.9)
    node(ax, 134, 55, 47, 17, ["Outcome: incident disease", "gallstone disease (primary) and",
                               "secondary outcomes; death as a", "competing event"], fc=tint(OUT, 0.12), ec=OUT, lw=0.9)
    node(ax, 74, 80, 50, 13, ["Visceral and liver fat (MRI)", "on the causal pathway; not adjusted"], fc=FILL_MED)
    node(ax, 2, 96, 60, 16, ["Health-seeking behaviour", "self-rated health, long-standing illness, number of",
                             "medications, cancer screening (sensitivity analysis)"])

    # edges
    a = dict(lw=0.8, ms=6)
    arrow(ax, 7, 20, 7, 55 - 0.3, **a)            # confounders -> exposure
    arrow(ax, 80, 20, 80, 29 - 0.3, **a)          # confounders -> BMI
    arrow(ax, 128, 20, 152, 55 - 0.3, **a)        # confounders -> outcome
    arrow(ax, 46, 36, 56 - 0.3, 36, **a)          # habits -> BMI
    arrow(ax, 30, 43, 30, 55 - 0.3, **a)          # habits -> exposure
    arrow(ax, 64, 43, 54, 55 - 0.3, **a)          # BMI -> exposure
    arrow(ax, 104, 38, 134 - 0.3, 58.5, **a)      # BMI -> outcome
    arrow(ax, 62, 63.5, 134 - 0.3, 63.5, color=EXP, lw=1.2, ms=7)   # exposure -> outcome
    arrow(ax, 50, 72, 74 - 0.3, 85, **a)          # exposure -> fat
    arrow(ax, 124, 86, 147, 72 + 0.3, **a)        # fat -> outcome
    arrow(ax, 24, 96, 24, 72 + 0.3, **a)          # health seeking -> exposure
    line(ax, [62, 162], [104, 104], lw=0.8)       # health seeking -> detection of the outcome
    arrow(ax, 162, 104, 162, 72 + 0.3, **a)
    T(ax, 112, 101.6, "Detection: diagnostic testing and recording of the outcome", size=S_TXT, color=GREY_TXT)

    # legend
    ly = 116.5
    items = [(dict(fc=FILL_ADJ, ec=EDGE, ls="-"), "Adjusted in all models"),
             (dict(fc="white", ec=EDGE, ls="-"), "Adjusted in sensitivity analyses"),
             (dict(fc="white", ec=EDGE, ls=(0, (2.5, 1.5))), "Unmeasured"),
             (dict(fc=FILL_MED, ec=EDGE, ls="-"), "Mediator, not adjusted")]
    x = 2.0
    for kw, lab in items:
        ax.add_patch(FancyBboxPatch((x, ly - 1.6), 6.0, 3.2, boxstyle="round,pad=0,rounding_size=0.5", lw=0.6,
                                    **kw))
        T(ax, x + 7.5, ly, lab, ha="left", size=S_TXT)
        x += 7.5 + len(lab) * 0.95 + 6.0

    plt.rcParams["savefig.pad_inches"] = 0
    save(fig, "FigS1_dag")
    print("saved")


if __name__ == "__main__":
    main()
