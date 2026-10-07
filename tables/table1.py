"""Table 1: characteristics of participants by accelerometer-measured sleep duration.

The population is the union of the five disease cohorts (each participant counted once). Values are n (%) or median
(interquartile range); percentages use participants with non-missing values as the denominator. Age is at day 0 (end
of accelerometer wear); other characteristics were measured at recruitment. Compositional means are the geometric
means of each part, closed to 24 hours.

Inputs: data/analysis_<outcome>.csv (06_build_analysis_datasets.py)
Output: tables/output/table1.csv
Usage:  python3 tables/table1.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parents[1]
DISEASES = ["gallstone", "peptic_ulcer", "gerd", "diverticular", "ibs"]
GROUPS = ["Sleep <7 h/day", "Sleep 7–9 h/day", "Sleep >9 h/day", "All"]
MINUS = "−"


def num(v, d=1):
    return f"{v:.{d}f}".replace("-", MINUS)


def median_iqr(s, d=1):
    s = s.dropna()
    lo, med, hi = s.quantile(0.25), s.median(), s.quantile(0.75)
    sep = " to " if lo < 0 or hi < 0 else "–"
    return f"{num(med, d)} ({num(lo, d)}{sep}{num(hi, d)})"


def n_pct(mask, valid):
    n, k = int(mask[valid].sum()), int(valid.sum())
    return f"{n:,} ({100 * n / k:.1f})"


def main():
    u = pd.concat([pd.read_csv(PROJ / "data" / f"analysis_{o}.csv") for o in DISEASES]).drop_duplicates("eid")
    u = u.rename(columns={"Sex": "sex", "Ethnic": "ethnic", "Townsend.deprivation.index.at.recruitment": "tdi",
                          "Education": "edu", "Body.mass.index..BMI....Instance.0": "bmi"})
    grp = np.select([u.sleep < 7, u.sleep <= 9], GROUPS[:2], GROUPS[2])
    parts = ["sleep", "inactivity_thr", "light_thr", "mvpa_thr"]
    rows = {}
    for g in GROUPS:
        d = u if g == "All" else u[grp == g]
        col = {"Participants, n (%)": f"{len(d):,}" if g == "All" else f"{len(d):,} ({100 * len(d) / len(u):.1f})",
               "Age, years": median_iqr(d.age_accel),
               "Women": n_pct(d.sex == 0, d.sex.notna()),
               "White": n_pct(d.ethnic == 1, d.ethnic.notna()),
               "Townsend deprivation index": median_iqr(d.tdi),
               "College or university degree": n_pct(d.edu == 1, d.edu.notna()),
               "Household income, £/year": ""}
        for k, lab in enumerate(["<18,000", "18,000–30,999", "31,000–51,999", "≥52,000"]):
            col[f"  {lab}"] = n_pct(d.income == k, d.income.notna())
        col["Retired"] = n_pct(d.retired == 1, d.retired.notna())
        col["Smoking"] = ""
        col["  Previous"] = n_pct(d.smoke == 1, d.smoke.notna())
        col["  Current"] = n_pct(d.smoke == 2, d.smoke.notna())
        col["Current drinker"] = n_pct(d.alcohol == 2, d.alcohol.notna())
        col["Healthy diet score (0–7)"] = median_iqr(d.diet_score)
        col["BMI, kg/m²"] = median_iqr(d.bmi)
        for c, lab in [("shift", "Shift work"), ("psych_meds", "Psychotropic medication"),
                       ("CVD", "Cardiovascular disease"), ("diabetes", "Diabetes"), ("hypertension", "Hypertension"),
                       ("Dyslipidemia", "Dyslipidaemia"), ("Cancer", "Cancer"), ("depression", "Depression"),
                       ("anxiety", "Anxiety"), ("sleep_disorder", "Sleep disorder")]:
            col[lab] = n_pct(d[c] == 1, d[c].notna())
        col["24-hour behaviours, h/day"] = ""
        for c, lab in zip(parts, ["Sleep", "Inactivity", "Light activity", "MVPA"]):
            col[f"  {lab}"] = median_iqr(d[c])
        gm = np.exp(np.log(d[parts]).mean())
        col["Compositional mean (sleep / inactivity / light activity / MVPA)"] = " / ".join(
            f"{v:.2f}" for v in 24 * gm / gm.sum())
        rows[g] = col
    tab = pd.DataFrame(rows)
    tab.index.name = "Characteristic"
    out = PROJ / "tables" / "output"
    out.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out / "table1.csv")
    print(f"union of the five cohorts: {len(u):,}; household income missing {int(u.income.isna().sum()):,} "
          f"({100 * u.income.isna().mean():.1f}%)")
    print(tab.to_string())


if __name__ == "__main__":
    main()
