"""Baseline covariates at the recruitment visit for participants with valid accelerometer data.

Definitions (Additional file 1, Table M1); column names are those used by the analysis scripts:
  Age.at.recruitment                          field 21022
  Sex                                         field 31: 0 female, 1 male
  Ethnic                                      field 21000: 1 White (codes 1, 1001-1003), 2 other; do not know or
                                              prefer not to answer -> missing
  Townsend.deprivation.index.at.recruitment   field 22189
  Education                                   field 6138: 1 college or university degree, 2 other; prefer not to
                                              answer -> missing
  Body.mass.index..BMI....Instance.0          field 21001
  smoke, alcohol                              fields 20116, 20117: 0 never, 1 previous, 2 current; prefer not to
                                              answer -> missing
  diet_score, healthy_diet                    seven components from fields 1289-1468 (see diet_score());
                                              healthy_diet = diet_score >= 4
  CVD, diabetes, hypertension, Dyslipidemia,  first-occurrence date before the recruitment visit (field 53):
  anxiety, depression, sleep_disorder         I20-I25, I60-I64 or I69; E10-E14; I10; E78; F40-F41; F32-F33; G47
  Cancer                                      cancer registry (fields 40005, 40006): C00-C97 except C44, diagnosed
                                              before the recruitment visit

Inputs:  data/raw_population_fields.csv, data/raw_assessment_fields.csv, data/raw_health_outcome_fields.csv
         (02_extract_fields.R)
Output:  data/baseline_covariates.csv
Usage:   python3 scripts/03_baseline_covariates.py
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[1] / "data"
CHRONIC = {  # first-occurrence date fields
    "CVD": [131296, 131298, 131300, 131302, 131304, 131306, 131360, 131362, 131364, 131366, 131368, 131378],
    "diabetes": [130706, 130708, 130710, 130712, 130714],
    "hypertension": [131286],
    "Dyslipidemia": [130814],
    "anxiety": [130904, 130906],
    "depression": [130894, 130896],
    "sleep_disorder": [131060],
}


def num(s):
    return pd.to_numeric(s, errors="coerce")


def tokens(x):
    """Codes from a multi-value cell such as '1|2', '[1,2]' or '1'."""
    if not isinstance(x, str) or not x.strip():
        return []
    return re.findall(r"-?\d+", x)


def counts(s):
    """Servings per day or week: 'less than one' (-10) -> 0; do not know / prefer not to answer -> missing."""
    return num(s).replace({-10: 0, -1: np.nan, -3: np.nan}).astype(float)


def weekly(s):
    """Frequency coding 100377 -> times per week."""
    return num(s).replace({-1: np.nan, -3: np.nan}).map({0: 0.0, 1: 0.5, 2: 1.0, 3: 3.0, 4: 5.5, 5: 7.0})


def component(value, healthy):
    return np.where(value.isna(), np.nan, healthy.astype(int))


def diet_score(ac):
    """Healthy diet score: one point for each of seven components (0-7).

    fruit (fresh + dried) >= 3 servings/day; vegetables (cooked + raw) >= 3 servings/day; fish (oily + non-oily)
    >= 2 times/week; processed meat <= 1 time/week; red meat (beef + lamb + pork) <= 1.5 times/week; whole grains
    >= 3 servings/day; refined grains <= 1.5 servings/day. Bread and cereal intakes are per week; whole grains are
    wholemeal bread (bread type 3) and bran, biscuit, oat or muesli cereals (cereal types 1-4); refined grains are
    white, brown or other bread (types 1, 2, 4) and other cereals (type 5). The score is missing only when all seven
    components are missing.
    """
    fruit = counts(ac.p1309_i0) + counts(ac.p1319_i0)
    veg = counts(ac.p1289_i0) + counts(ac.p1299_i0)
    fish = weekly(ac.p1329_i0) + weekly(ac.p1339_i0)
    processed = weekly(ac.p1349_i0)
    red = weekly(ac.p1369_i0) + weekly(ac.p1379_i0) + weekly(ac.p1389_i0)
    bread_day, cereal_day = counts(ac.p1438_i0) / 7.0, counts(ac.p1458_i0) / 7.0
    bread_type = num(ac.p1448_i0).replace({-1: np.nan, -3: np.nan})
    cereal_type = num(ac.p1468_i0).replace({-1: np.nan, -3: np.nan})
    whole = (np.where(bread_type == 3, bread_day, np.where(bread_type.isna(), np.nan, 0.0))
             + np.where(cereal_type.isin([1, 2, 3, 4]), cereal_day, np.where(cereal_type.isna(), np.nan, 0.0)))
    refined = (np.where(bread_type.isin([1, 2, 4]), bread_day, np.where(bread_type.isna(), np.nan, 0.0))
               + np.where(cereal_type == 5, cereal_day, np.where(cereal_type.isna(), np.nan, 0.0)))
    whole, refined = pd.Series(whole, index=ac.index), pd.Series(refined, index=ac.index)
    comp = pd.DataFrame({
        "fruit": component(fruit, fruit >= 3), "vegetables": component(veg, veg >= 3),
        "fish": component(fish, fish >= 2), "processed_meat": component(processed, processed <= 1),
        "red_meat": component(red, red <= 1.5), "whole_grains": component(whole, whole >= 3),
        "refined_grains": component(refined, refined <= 1.5)})
    return comp.sum(axis=1, min_count=1)


def main():
    po = pd.read_csv(DATA / "raw_population_fields.csv", dtype=str)
    ac = pd.read_csv(DATA / "raw_assessment_fields.csv", dtype=str)
    ho = pd.read_csv(DATA / "raw_health_outcome_fields.csv", dtype=str)
    d = po[["eid"]].merge(ac[["eid"]], on="eid").merge(ho[["eid"]], on="eid")
    po, ac, ho = (x.set_index("eid").loc[d.eid] for x in (po, ac, ho))
    out = pd.DataFrame({"eid": d.eid.astype(int).values})
    out["Age.at.recruitment"] = num(po.p21022).values
    out["Sex"] = num(po.p31).values

    eth = num(ac.p21000_i0)
    out["Ethnic"] = np.where(eth.isin([1, 1001, 1002, 1003]), 1.0, np.where(eth > 0, 2.0, np.nan))
    out["Townsend.deprivation.index.at.recruitment"] = num(po.p22189).values

    edu = ac.p6138_i0.map(tokens)
    out["Education"] = np.where(edu.map(lambda t: "1" in t), 1.0,
                                np.where(edu.map(lambda t: len(t) == 0 or t == ["-3"]), np.nan, 2.0))
    out["Body.mass.index..BMI....Instance.0"] = num(ac.p21001_i0).values
    out["smoke"] = num(ac.p20116_i0).where(lambda x: x >= 0).values
    out["alcohol"] = num(ac.p20117_i0).where(lambda x: x >= 0).values
    score = diet_score(ac)
    out["diet_score"] = score.values
    out["healthy_diet"] = np.where(score.isna(), np.nan, (score >= 4).astype(int))

    recruit = pd.to_datetime(ac.p53_i0, errors="coerce")
    for name, fields in CHRONIC.items():
        first = pd.concat([pd.to_datetime(ho[f"p{f}"], errors="coerce") for f in fields], axis=1).min(axis=1)
        out[name] = (first < recruit).astype(int).values

    n_inst = len([c for c in ho.columns if c.startswith("p40006_i")])
    dates = []
    for i in range(n_inst):
        icd = ho[f"p40006_i{i}"].fillna("").str.strip().str.upper()
        code = pd.to_numeric(icd.str.extract(r"^C(\d{2})")[0], errors="coerce")
        target = code.between(0, 97) & (code != 44)
        dates.append(pd.to_datetime(ho[f"p40005_i{i}"], errors="coerce").where(target))
    first_cancer = pd.concat(dates, axis=1).min(axis=1)
    out["Cancer"] = (first_cancer < recruit).astype(int).values
    out = out[["eid", "Age.at.recruitment", "Sex", "Ethnic", "Townsend.deprivation.index.at.recruitment",
               "Education", "Body.mass.index..BMI....Instance.0", "smoke", "alcohol", "diet_score", "healthy_diet",
               "CVD", "diabetes", "hypertension", "Dyslipidemia", "Cancer", "anxiety", "depression",
               "sleep_disorder"]]
    out.to_csv(DATA / "baseline_covariates.csv", index=False)
    print(f"{len(out)} participants; missing values:")
    print(out.isna().sum()[lambda s: s > 0].to_string())


if __name__ == "__main__":
    main()
