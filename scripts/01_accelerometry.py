"""Accelerometer quality control and the 24-hour composition for every participant who returned accelerometer data.

- Device quality control: good wear time (field 90015), good calibration (90016), calibrated on own data (90017),
  mean acceleration <100 mg (90012) and no data-problem flag (90002).
- Oxford machine-learning behaviours (fields 40034-40041): weekday and weekend hourly averages, weighted 5:2;
  compositions must be complete, non-negative and sum to 23-25 h, and are closed to 24 h; compositions with any
  behaviour below the 0.1th or above the 99.9th percentile are removed.
- Exposure used in the analyses: machine-learning sleep plus acceleration thresholds for waking time (fields 90113 and
  90127, fractions of wear time at <=30 mg and <=100 mg): MVPA >100 mg, light activity >30 to <=100 mg, inactivity
  <=30 mg minus sleep.
- Day 0 (index date) = end of accelerometer wear (field 90011); age at day 0 from year and month of birth (fields 34,
  52; day set to the 15th); date lost to follow-up (field 191).

Inputs (config.sh): UKB_ACCELEROMETER, UKB_POPULATION
Outputs: data/accel_base.csv, data/accel_base_flow.csv
Usage:   python3 scripts/01_accelerometry.py
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[1] / "data"
BEHAVIOURS = ["sleep", "sedentary", "light", "mvpa"]


def hourly_sum(s):
    """24 comma- or bar-separated hourly fractions -> hours per day."""
    return s.map(lambda x: sum(map(float, x.split("|" if "|" in x else ","))) if isinstance(x, str) else np.nan)


def main():
    DATA.mkdir(exist_ok=True)
    wd_f = dict(zip(BEHAVIOURS, ["p40034_i0", "p40035_i0", "p40036_i0", "p40037_i0"]))
    we_f = dict(zip(BEHAVIOURS, ["p40038_i0", "p40039_i0", "p40040_i0", "p40041_i0"]))
    qc_f = ["p90002", "p90010", "p90011", "p90012", "p90015", "p90016", "p90017", "p90051", "p90113", "p90127"]
    acc = pd.read_csv(os.environ["UKB_ACCELEROMETER"], usecols=["eid", *wd_f.values(), *we_f.values(), *qc_f],
                      low_memory=False)
    acc = acc[acc.p90011.notna()].copy()
    flow = {"accelerometer data returned": len(acc)}

    for b in BEHAVIOURS:
        acc[b] = (5 * hourly_sum(acc[wd_f[b]]) + 2 * hourly_sum(acc[we_f[b]])) / 7
    acc = acc[(acc.p90015 == 1) & (acc.p90016 == 1) & (acc.p90017 == 1) & (acc.p90012 < 100) & acc.p90002.isna()]
    flow["passed device QC"] = len(acc)
    acc = acc.dropna(subset=BEHAVIOURS)
    acc = acc[(acc[BEHAVIOURS] >= 0).all(axis=1)]
    acc = acc[acc[BEHAVIOURS].sum(axis=1).between(23, 25)]
    acc[BEHAVIOURS] = acc[BEHAVIOURS].div(acc[BEHAVIOURS].sum(axis=1), axis=0) * 24
    flow["complete 4-part composition"] = len(acc)
    lo, hi = acc[BEHAVIOURS].quantile(0.001), acc[BEHAVIOURS].quantile(0.999)
    acc = acc[((acc[BEHAVIOURS] >= lo) & (acc[BEHAVIOURS] <= hi)).all(axis=1)].copy()
    flow["after trimming extreme compositions"] = len(acc)

    acc["mvpa_thr"] = (1 - acc.p90127) * 24
    acc["light_thr"] = (acc.p90127 - acc.p90113) * 24
    acc["inactivity_thr"] = acc.p90113 * 24 - acc.sleep
    acc["index_date"] = pd.to_datetime(acc.p90011).dt.normalize()

    pop = pd.read_csv(os.environ["UKB_POPULATION"], usecols=["eid", "p34", "p52", "p191"])
    acc = acc.merge(pop, on="eid", how="left")
    birth = pd.to_datetime(dict(year=acc.p34, month=acc.p52, day=15), errors="coerce")
    acc["age_accel"] = (acc.index_date - birth).dt.days / 365.25
    acc["lost_date"] = pd.to_datetime(acc.p191, errors="coerce")
    keep = ["eid", "age_accel", *BEHAVIOURS, "inactivity_thr", "light_thr", "mvpa_thr", "index_date", "lost_date"]
    acc[keep].to_csv(DATA / "accel_base.csv", index=False)
    pd.Series(flow).to_csv(DATA / "accel_base_flow.csv", header=["n"])
    print(pd.Series(flow).to_string())


if __name__ == "__main__":
    main()
