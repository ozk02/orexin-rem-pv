"""Step 6 (FAERS only): REM-related reports with onset after the end of therapy (Table 5).

Onset after the end of therapy = event_dt (report level) later than end_dt of the suspect drug of that class;
only reports in which both dates are complete (YYYYMMDD) are classified. Also counts REM-related reports
co-reported with withdrawal/rebound PTs.

Usage: python src/faers_post_discontinuation.py   -> output/faers_post_discontinuation.csv / .txt
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from build_data import load_faers
from config import FAERS_PT, OUT_DIR

CLASSES = {"DORA": "DORA", "BZZ": "BZ/Z", "RAM": "Ramelteon", "SSRI_SNRI": "SSRI/SNRI",
           "TCA": "Tricyclic antidepressants", "MIR_TRZ": "Mirtazapine/trazodone"}
HALL = FAERS_PT["other_halluc"] + FAERS_PT["hypnagogic"]
WD = "withdrawal|discontinuation syndrome|rebound effect|rebound insomnia|rebound anxiety"

demo, reac, drug, ther = load_faers()
sus = drug[drug["role_cod"].isin(["PS", "SS"])]
df = demo[["event_dt"]].copy()
df["REM"] = df.index.isin(reac.loc[reac["pt"].isin(FAERS_PT["REM"]), "primaryid"]).astype(int)
df["halluc"] = df.index.isin(reac.loc[reac["pt"].isin(HALL), "primaryid"]).astype(int)
ev = df[(df.REM == 1) | (df.halluc == 1)].copy()
ev["onset"] = pd.to_datetime(ev["event_dt"].where(ev["event_dt"].str.len() == 8), format="%Y%m%d", errors="coerce")
th = ther.copy()
th["end"] = pd.to_datetime(th["end_dt"].where(th["end_dt"].str.len() == 8), format="%Y%m%d", errors="coerce")
wd_ids = set(reac.loc[reac["pt"].str.contains(WD, na=False), "primaryid"])

rows = []
for k, name in CLASSES.items():
    seqs = sus.loc[sus["cls"] == k, ["primaryid", "drug_seq"]].drop_duplicates()
    e = ev[ev.index.isin(seqs["primaryid"])]
    end = th.merge(seqs, on=["primaryid", "drug_seq"]).groupby("primaryid")["end"].max()
    e = e.join(end, how="left")
    known = e.dropna(subset=["onset", "end"])
    after = known[known["onset"] > known["end"]]
    gap = (after["onset"] - after["end"]).dt.days
    rows.append(dict(drug_class=name, REM_reports=int(e.REM.sum()), REM_complete_dates=int(known.REM.sum()),
                     REM_after_end=int(after.REM.sum()), halluc_reports=int(e.halluc.sum()),
                     halluc_complete_dates=int(known.halluc.sum()), halluc_after_end=int(after.halluc.sum()),
                     median_days_after_end=float(gap.median()) if len(gap) else np.nan,
                     REM_with_withdrawal_PT=int(e[e.REM == 1].index.isin(wd_ids).sum())))
t = pd.DataFrame(rows)
t.to_csv(OUT_DIR / "faers_post_discontinuation.csv", index=False)

lines = [t.to_string(index=False), ""]
r = t.set_index("drug_class")
for comp in ["SSRI/SNRI", "BZ/Z"]:
    for a_col, n_col, label in [("REM_after_end", "REM_complete_dates", "onset after end of therapy"),
                                ("REM_with_withdrawal_PT", "REM_reports", "co-reported withdrawal/rebound PT")]:
        a, n1 = r.loc[comp, a_col], r.loc[comp, n_col]
        b, n2 = r.loc["DORA", a_col], r.loc["DORA", n_col]
        o, p = fisher_exact([[a, n1 - a], [b, n2 - b]])
        lines.append(f"{comp} vs DORA, {label}: {a}/{n1} vs {b}/{n2}, OR {o:.1f}, Fisher p={p:.2g}")
print("\n".join(lines))
(OUT_DIR / "faers_post_discontinuation.txt").write_text("\n".join(lines), encoding="utf-8")
