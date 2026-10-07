"""Analysis dataset for each outcome: link the outcome to the accelerometer cohort, apply the exclusions, derive the
follow-up and add the covariates. Every outcome uses the same exposure and covariates.

Exclusions, in order: disease at recruitment (removed by FreeUKB in 05_define_outcomes.R); disease between recruitment
and day 0; death or loss to follow-up before day 0; outcome-specific exclusions using records before day 0:
  gallstone            previous cholecystectomy (OPCS-4 J18, or self-reported operation code 1455 or 1528)
  ibs                  inflammatory bowel disease (K50, K51), intestinal malabsorption (K90) or colorectal cancer
                       (C18-C20, cancer registry)
  nco_earwax           hearing-aid use at recruitment (field 3393)
Follow-up runs from day 0 to the earliest of onset, death, loss to follow-up and administrative censoring. In the
gallstone cohorts, participants who undergo cholecystectomy without a gallstone record are censored at surgery.
gallstone_surgery (cholecystectomy endpoint) uses the gallstone cohort; an event is a first cholecystectomy during
follow-up with gallstone disease recorded no later than 30 days after the operation.
peptic_ulcer is the first of gastric ulcer (K25) and duodenal ulcer (K26).

The outcome event and follow-up are stored in the columns gs / time_to_gs for every outcome, so that the modelling
scripts work unchanged.

Inputs:  data/accel_base.csv (01), data/raw_assessment_fields.csv and data/raw_health_outcome_fields.csv (02),
         data/baseline_covariates.csv (03), data/hearing_baseline.csv (04), data/outcomes/<name>/out.csv (05),
         code_lists/*_codes.csv
Outputs: data/analysis_<outcome>.csv, data/flow_<outcome>.csv
Usage:   python3 scripts/06_build_analysis_datasets.py [outcome ...]
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parents[1]
DATA = PROJ / "data"
CODES = PROJ / "code_lists"

OUTCOMES = {
    "gallstone": {"files": ["gallstone"], "exclude": ["cholecystectomy"]},
    # sensitivity outcome: same cohort as gallstone, event = cholecystectomy for gallstone disease
    "gallstone_surgery": {"files": ["gallstone"], "exclude": ["cholecystectomy"], "event": "cholecystectomy"},
    "gerd": {"files": ["gerd"], "exclude": []},
    "diverticular": {"files": ["diverticular"], "exclude": []},
    "ibs": {"files": ["ibs"], "exclude": ["ibd", "malabsorption", "colorectal_cancer"]},
    "peptic_ulcer": {"files": ["gastric_ulcer", "duodenal_ulcer"], "exclude": []},
    "nco_cyst": {"files": ["cyst"], "exclude": []},
    "nco_earwax": {"files": ["earwax"], "exclude": ["hearing_aid"]},
}
BEHAVIOURS = ["sleep", "sedentary", "light", "mvpa"]


# ------------------------------------------------------------------ covariates
def tokens(x):
    if not isinstance(x, str) or not x.strip():
        return []
    return re.findall(r"-?[A-Z]?\d+", x)


def has_any(series, codes):
    codes = set(codes)
    return series.map(lambda x: bool(codes & set(tokens(x))))


def derive_covariates(ac, female):
    """income (0-3), retired, sick_disabled, shift, psych_meds, aspirin, nsaid, acid_rx, antacid, chol_med, hrt,
    ocp, parity."""
    cov = pd.DataFrame({"eid": ac.eid})
    # household income: <18k, 18-31k, 31-52k, >=52k (the top two bands merged); do not know / prefer not to answer -> missing
    cov["income"] = pd.to_numeric(ac.p738_i0, errors="coerce").map({1: 0, 2: 1, 3: 2, 4: 3, 5: 3})
    emp = ac.p6142_i0
    cov["retired"] = has_any(emp, ["2"]).astype(int)
    cov["sick_disabled"] = has_any(emp, ["4"]).astype(int)
    cov.loc[emp.isna() | has_any(emp, ["-3"]), ["retired", "sick_disabled"]] = np.nan
    # shift work was asked only of those in paid employment; not employed -> no shift work
    shift = pd.to_numeric(ac.p826_i0, errors="coerce")
    cov["shift"] = np.select([shift.isin([2, 3, 4]), shift == 1], [1, 0], default=np.nan)
    cov.loc[shift.isna() & ~has_any(emp, ["1"]) & emp.notna(), "shift"] = 0
    med_cols = [c for c in ac.columns if c.startswith("p20003_i0_")]
    psych = set(pd.read_csv(CODES / "psych_med_codes.csv", dtype=str).code)
    cov["psych_meds"] = ac[med_cols].isin(psych).any(axis=1).astype(int)
    # aspirin / oral NSAID use (peptic ulcer adjustment) and acid-suppressant / antacid use (exclusion analysis)
    for name in ["aspirin", "nsaid", "acid_rx", "antacid"]:
        codes = set(pd.read_csv(CODES / f"{name}_codes.csv", dtype=str).code)
        cov[name] = ac[med_cols].isin(codes).any(axis=1).astype(int)
    # cholesterol-lowering medication (women answer field 6153, men field 6177)
    meds = ac.p6153_i0.fillna(ac.p6177_i0)
    cov["chol_med"] = has_any(meds, ["1"]).astype(int)
    cov.loc[meds.isna() | has_any(meds, ["-1", "-3"]), "chol_med"] = np.nan
    # female reproductive factors, coded 0 for men
    yes_no = lambda s: pd.to_numeric(s, errors="coerce").map({1: 1, 0: 0})  # noqa: E731
    fem = female.reindex(ac.eid).values
    cov["hrt"] = np.where(fem, yes_no(ac.p2814_i0), 0)
    cov["ocp"] = np.where(fem, yes_no(ac.p2784_i0), 0)
    births = pd.to_numeric(ac.p2734_i0, errors="coerce").where(lambda x: x >= 0)
    cov["parity"] = np.where(fem, births.clip(upper=4), 0)  # 4 = four or more
    return cov


# ------------------------------------------------------------------ outcomes
def first_date_from_array(ho, code_prefix):
    """Earliest OPCS-4 date for codes starting with code_prefix (field 41272 aligned with 41282 arrays)."""
    date_cols = sorted([c for c in ho.columns if c.startswith("p41282_a")], key=lambda c: int(c.split("_a")[1]))

    def one(row):
        ops = tokens(row["p41272"])
        hits = [pd.to_datetime(row[c], errors="coerce") for op, c in zip(ops, date_cols[: len(ops)])
                if op.startswith(code_prefix)]
        hits = [h for h in hits if pd.notna(h)]
        return min(hits) if hits else pd.NaT

    return ho.apply(one, axis=1)


def load_outcome(names):
    merged = None
    for n in names:
        f = DATA / "outcomes" / n / "out.csv"
        h = pd.read_csv(f, nrows=0).columns
        ev = [c for c in h if c.endswith("_event_date")][0]
        adm = [c for c in h if c.endswith("_admin_censor_date")][0]
        o = pd.read_csv(f, usecols=["Participant.ID", ev, adm, "death_date"]).rename(
            columns={"Participant.ID": "eid", ev: "ev", adm: "adm"})
        for c in ["ev", "adm", "death_date"]:
            o[c] = pd.to_datetime(o[c], errors="coerce")
        if merged is None:
            merged = o
        else:  # composite outcome: first of the component diagnoses
            merged = merged.merge(o, on="eid", suffixes=("", "_2"))  # inner: drops baseline cases of either
            merged["ev"] = merged[["ev", "ev_2"]].min(axis=1)
            merged["adm"] = merged[["adm", "adm_2"]].min(axis=1)
            merged["death_date"] = merged.death_date.fillna(merged.death_date_2)
            merged = merged[["eid", "ev", "adm", "death_date"]]
    return merged.rename(columns={"ev": "gs_date", "adm": "admin_censor"})


def build(names):
    base = pd.read_csv(DATA / "accel_base.csv", parse_dates=["index_date", "lost_date"])
    base_cov = pd.read_csv(DATA / "baseline_covariates.csv")
    ac = pd.read_csv(DATA / "raw_assessment_fields.csv", dtype=str)
    ac["eid"] = ac.eid.astype(int)
    ho = pd.read_csv(DATA / "raw_health_outcome_fields.csv", dtype=str)
    ho["eid"] = ho.eid.astype(int)
    female = base_cov.set_index("eid").Sex.eq(0)
    extra = derive_covariates(ac, female)

    dates = pd.DataFrame({"eid": ho.eid})
    for col, name in [("p131626", "k50"), ("p131628", "k51"), ("p131688", "k90")]:
        dates[name] = pd.to_datetime(ho[col], errors="coerce")
    dates["j18"] = first_date_from_array(ho, "J18")
    crc = []  # first colorectal cancer (C18-C20) in the cancer registry
    for i in range(len([c for c in ho.columns if c.startswith("p40006_i")])):
        is_crc = ho[f"p40006_i{i}"].astype(str).str.match(r"^C(18|19|20)")
        crc.append(pd.to_datetime(ho[f"p40005_i{i}"], errors="coerce").where(is_crc))
    dates["crc"] = pd.concat(crc, axis=1).min(axis=1)
    self_chole = ac[[c for c in ac.columns if c.startswith("p20004_i0_")]].isin({"1455", "1528"}).any(axis=1)
    dates = dates.merge(pd.DataFrame({"eid": ac.eid, "self_chole": self_chole}), on="eid", how="left")

    summary = []
    for name in names:
        cfg = OUTCOMES[name]
        flow = {}
        d = base.merge(load_outcome(cfg["files"]), on="eid", how="inner")
        flow["linked to outcome data (no disease at recruitment)"] = len(d)
        d = d[~(d.gs_date <= d.index_date)]
        flow["excluding disease between recruitment and accelerometry"] = len(d)
        d = d[~(d.lost_date <= d.index_date) & ~(d.death_date <= d.index_date)]
        flow["excluding death/loss before accelerometry"] = len(d)
        d = d.merge(dates, on="eid", how="left")
        before = lambda c: d[c] <= d.index_date  # noqa: E731
        if "ibd" in cfg["exclude"]:
            d = d[~(before("k50") | before("k51"))]
        if "malabsorption" in cfg["exclude"]:
            d = d[~before("k90")]
        if "colorectal_cancer" in cfg["exclude"]:
            d = d[~before("crc")]
        if "cholecystectomy" in cfg["exclude"]:
            d = d[~(d.self_chole.fillna(False) | before("j18"))]
        if "hearing_aid" in cfg["exclude"]:  # hearing aids cause ear-wax impaction
            hear = pd.read_csv(DATA / "hearing_baseline.csv")
            d = d[~d.eid.isin(hear.eid[hear.hearing_aid == 1])]
        if cfg["exclude"]:
            flow["outcome-specific exclusions (" + ", ".join(cfg["exclude"]) + ")"] = len(d)
        if cfg.get("event") == "cholecystectomy":
            # cholecystectomy with gallstone disease recorded by 30 days after the operation; gallstones that are
            # never operated on do not count as events
            j = d.j18.where(d.j18 > d.index_date)
            d["gs_date"] = j.where(d.gs_date.notna() & (d.gs_date <= j + pd.Timedelta(days=30)))

        end = d[["admin_censor", "death_date", "lost_date"]].min(axis=1)
        d["gs"] = (d.gs_date.notna() & (d.gs_date <= end)).astype(int)
        d["death"] = ((d.gs == 0) & d.death_date.notna() & (d.death_date <= end)).astype(int)
        d["end_date"] = pd.to_datetime(np.where(d.gs == 1, d.gs_date, end))
        d["time_to_gs"] = (d.end_date - d.index_date).dt.days
        if "cholecystectomy" in cfg["exclude"]:  # no longer at risk after surgery
            j = (d.j18 - d.index_date).dt.days
            cens = (d.gs == 0) & j.notna() & (j > 0) & (j < d.time_to_gs)
            d.loc[cens, "time_to_gs"] = j[cens]
            d.loc[cens, "death"] = 0
        d["time_to_death"] = np.where(d.death == 1, d.time_to_gs, np.nan)
        d = d[d.time_to_gs > 0]
        flow["analytic sample"] = len(d)

        d = d.merge(base_cov, on="eid", how="left").merge(extra, on="eid", how="left")
        keep = ["eid", "age_accel", *BEHAVIOURS, "inactivity_thr", "light_thr", "mvpa_thr", "index_date", "end_date",
                "gs", "death", "time_to_gs", "time_to_death",
                *[c for c in base_cov.columns if c != "eid"], *[c for c in extra.columns if c != "eid"]]
        d = d[keep].rename(columns={"sedentary": "inactivity"})
        out = DATA / f"analysis_{name}.csv"
        tmp = out.with_suffix(".tmp")
        d.to_csv(tmp, index=False)
        tmp.replace(out)  # atomic, so a running analysis never reads a half-written file
        pd.Series(flow).to_csv(DATA / f"flow_{name}.csv", header=["n"])
        g = np.select([d.sleep < 7, d.sleep <= 9], ["<7", "7-9"], ">9")
        summary.append({"outcome": name, "N": len(d), "cases": int(d.gs.sum()), "deaths": int(d.death.sum()),
                        **{f"cases {k}": int(d.gs[g == k].sum()) for k in ["<7", "7-9", ">9"]},
                        "median FU": round(d.time_to_gs.median() / 365.25, 2)})
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    build(sys.argv[1:] or list(OUTCOMES))
