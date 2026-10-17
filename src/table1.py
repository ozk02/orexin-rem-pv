"""Step 5: characteristics of DORA-only and BZ/Z-only reports (Table 1) and report counts (flow).

Usage: python src/table1.py jader|faers   -> output/{db}_table1.txt
"""
import glob
import sys

import pandas as pd

from build_data import build
from config import FAERS_DIR, OUT_DIR

db = sys.argv[1]
df, _ = build(db)
lines = []


def P(s=""):
    print(s, flush=True)
    lines.append(s)


if db == "faers":
    raw = sum(len(pd.read_parquet(f, columns=["primaryid"]))
              for f in glob.glob(str(FAERS_DIR / "quarters" / "*_demo.parquet")))
    P(f"reports before deduplication {raw:,} -> after {len(df):,}")
else:
    P(f"cases from fiscal year 2014 onwards {len(df):,}")
P(f"suspect DORA {int(df.DORA.sum()):,} / BZ/Z {int(df.BZZ.sum()):,} / ramelteon {int(df.RAM.sum()):,} / "
  f"DORA and BZ/Z both {int((df.DORA & df.BZZ).sum()):,}")

groups = {"DORA-only": df[(df.DORA == 1) & (df.BZZ == 0) & (df.RAM == 0)],
          "BZ/Z-only": df[(df.DORA == 0) & (df.BZZ == 1) & (df.RAM == 0)],
          "all reports": df}
for col in ["age", "sex", "reporter"] + (["country"] if db == "faers" else []):
    P(f"\n[{col}]")
    for name, g in groups.items():
        v = g[col].value_counts()
        if col == "country":
            v = pd.concat([v[v.index.isin(["US", "JP"])], pd.Series({"other": v[~v.index.isin(["US", "JP"])].sum()})])
        P(f"  {name} (n={len(g):,}): " + ", ".join(f"{k} {c:,} ({c / len(g) * 100:.1f}%)" for k, c in v.items()))
P("\n[drug indicators]")
for c in [c for c in df.columns if c.startswith("conc_")]:
    P(f"  {c[5:]}: " + ", ".join(f"{n} {int(g[c].sum()):,} ({g[c].mean() * 100:.1f}%)" for n, g in groups.items()))
(OUT_DIR / f"{db}_table1.txt").write_text("\n".join(lines), encoding="utf-8")
