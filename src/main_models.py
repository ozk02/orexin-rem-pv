"""Step 3: primary analyses (Tables 2 and 3 of the manuscript).

Model A: all other reports as reference; DORA, BZ/Z and ramelteon entered together.
Model B (primary): active comparator, DORA-only vs BZ/Z-only (reports with both, or with ramelteon, excluded).
Firth logistic regression; basic adjustment = age, sex; full adjustment = + reporter type and drug indicators.
Trend test across the event tiers REM_core -> hypnagogic -> other_halluc.

Usage: python src/main_models.py jader|faers   -> output/{db}_main.txt
"""
import sys

from build_data import build
from config import EXPOSURES, LADDER, OUT_DIR, OUTCOMES, REF_LEVELS
from firth import fit_ror, gradient_test, ratio_test


def fmt(t):
    return f"{t[0]:.1f} ({t[1]:.1f}-{t[2]:.1f})"


def covariates(df):
    conc = [c for c in df.columns if c.startswith("conc_")]
    return ["age", "sex"], ["age", "sex", "reporter"] + conc


def model_b_subset(df):
    return df[((df["DORA"] + df["BZZ"]) == 1) & (df["RAM"] == 0)]


def run(df, label):
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)

    basic, full = covariates(df)
    P(f"=== {label} ===  reports {len(df):,}")
    P("Drug indicators (reports): " + ", ".join(f"{c[5:]} {int(df[c].sum()):,}" for c in full[3:]))
    for o in OUTCOMES:
        P(f"  {o}: total {int(df[o].sum()):,} | " + " ".join(f"{e} {int((df[e] & df[o]).sum())}" for e in EXPOSURES))

    P("\n[A] All other reports as reference (three groups entered together), Firth. basic = age, sex; "
      "full = + reporter, drug indicators")
    for o in OUTCOMES:
        rb = fit_ror(df, o, EXPOSURES, basic, REF_LEVELS)
        rf = fit_ror(df, o, EXPOSURES, full, REF_LEVELS)
        for e in EXPOSURES:
            P(f"  {o:13s} {e:5s} basic {fmt(rb[e]):24s} full {fmt(rf[e])}")

    P("\n[B] Active comparator: DORA-only vs BZ/Z-only (Firth, full adjustment)")
    h = model_b_subset(df)
    res = {}
    for o in OUTCOMES:
        a, c = int(h.loc[h.DORA == 1, o].sum()), int(h.loc[h.DORA == 0, o].sum())
        res[o] = fit_ror(h, o, ["DORA"], full, REF_LEVELS)["DORA"]
        P(f"  {o:13s} DORA {a}/{int(h.DORA.sum()):,}  BZ/Z {c}/{int((h.DORA == 0).sum()):,}  ROR {fmt(res[o])}")

    P("\n[Trend] REM_core -> hypnagogic -> other_halluc (Model B)")
    g = gradient_test([res[o][3] for o in LADDER], [res[o][4] for o in LADDER])
    P(f"  ROR ratio per step {g['ratio_per_step']:.2f} (slope {g['slope']:.2f}, SE {g['se']:.2f})  "
      f"p_trend={g['p_trend']:.2g}  Cochran Q={g['Q']:.1f} p={g['p_het']:.2g}")
    for i, j in [(0, 2), (1, 2), (0, 1)]:
        a, b = LADDER[i], LADDER[j]
        rt = ratio_test(res[a][3], res[a][4], res[b][3], res[b][4])
        P(f"  {a} / {b}: {rt[0]:.1f} ({rt[1]:.1f}-{rt[2]:.1f})  p={rt[3]:.2g}")

    P("\n[Trend, Model A] DORA ROR against all other reports (full adjustment)")
    ra = {o: fit_ror(df, o, EXPOSURES, full, REF_LEVELS)["DORA"] for o in LADDER}
    for o in LADDER:
        P(f"  {o:13s} {fmt(ra[o])}")
    g2 = gradient_test([ra[o][3] for o in LADDER], [ra[o][4] for o in LADDER])
    P(f"  ratio per step {g2['ratio_per_step']:.2f}  p_trend={g2['p_trend']:.2g}  p_het={g2['p_het']:.2g}")
    return lines


if __name__ == "__main__":
    db = sys.argv[1]
    df, _ = build(db)
    label = "JADER (fiscal year 2014 onwards)" if db == "jader" else "FAERS (2014Q3-2026Q2)"
    (OUT_DIR / f"{db}_main.txt").write_text("\n".join(run(df, label)), encoding="utf-8")
