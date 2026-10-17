"""Step 4: sensitivity analyses of Model B (Table 4) and the name-matching check (Supplementary Table S1).

(1) excluding consumer reports / healthcare professionals only
(2) by reporting period (notoriety bias)
(3) name-matching check: most frequent names captured by each drug pattern
(4) by country of occurrence (FAERS)

Usage: python src/sensitivity.py jader|faers   -> output/{db}_sensitivity.txt / .csv
"""
import sys

import pandas as pd

from build_data import JADER_CLASS, JADER_CONC, build
from config import LADDER, OUT_DIR, REF_LEVELS
from firth import fit_ror, gradient_test
from main_models import covariates, model_b_subset

db = sys.argv[1]
df, drug = build(db)
lines = []


def P(s=""):
    print(s, flush=True)
    lines.append(s)


_, full = covariates(df)
base = model_b_subset(df)
yr = base["year"].astype(int)
subsets = [("primary", base.index == base.index),
           ("(1) excluding consumer reports", base["reporter"] != "consumer"),
           ("(1) healthcare professionals only", base["reporter"] == "HCP"),
           ("(2) 2014-2017", yr <= 2017),
           ("(2) 2018-2021", (yr >= 2018) & (yr <= 2021)),
           ("(2) 2022-2026", yr >= 2022)]
if db == "faers":
    subsets += [("(4) USA only", base["country"] == "US"), ("(4) Japan only", base["country"] == "JP"),
                ("(4) outside USA", base["country"] != "US")]

P(f"=== Sensitivity analyses, {db.upper()}: DORA-only vs BZ/Z-only (Firth, full adjustment) ===")
P("subset | DORA reports | " + " | ".join(f"{o}: events DORA/BZZ, ROR (95% CI)" for o in LADDER) + " | p_trend")
rows = []
for name, m in subsets:
    h = base[m]
    cells, fits = [], []
    for o in LADDER:
        a, c = int(h.loc[h.DORA == 1, o].sum()), int(h.loc[h.DORA == 0, o].sum())
        r = fit_ror(h, o, ["DORA"], full, REF_LEVELS)["DORA"]
        fits.append(r)
        cells.append(f"{a}/{c} {r[0]:.1f} ({r[1]:.1f}-{r[2]:.1f})")
    g = gradient_test([x[3] for x in fits], [x[4] for x in fits])
    P(f"{name} | {int(h.DORA.sum()):,} | " + " | ".join(cells) + f" | {g['p_trend']:.2g}")
    rows.append(dict(subset=name, dora_reports=int(h.DORA.sum()), **dict(zip(LADDER, cells)), p_trend=g["p_trend"]))

# (3) Name-matching check. Counts are drug records in any role.
P("\n(3) Name-matching check (most frequent matched names; counts = drug records)")
if db == "jader":
    gn = drug["医薬品（一般名）"].fillna("")
    for k, p in {**JADER_CLASS, **JADER_CONC}.items():
        v = gn[gn.str.contains(p, regex=True)].value_counts()
        P(f"  {k} ({len(v)} names): " + "; ".join(f"{n} {c}" for n, c in v.head(40).items()))
else:
    d = drug.copy()
    d["name"] = d["prod_ai"].fillna(d["drugname"]).str.upper().str.strip()
    for k, g in d.groupby("cls"):
        v = g["name"].value_counts()
        P(f"  {k} ({len(v)} names): " + "; ".join(f"{n} {c}" for n, c in v.head(25).items()))

pd.DataFrame(rows).to_csv(OUT_DIR / f"{db}_sensitivity.csv", index=False)
(OUT_DIR / f"{db}_sensitivity.txt").write_text("\n".join(lines), encoding="utf-8")
