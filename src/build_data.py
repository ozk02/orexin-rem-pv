"""Step 2: build one analysis table per database (one row per deduplicated report).

Columns: DORA, BZZ, RAM (1 = suspect drug), conc_* (covariate drug indicators), age, sex, reporter,
year, country, and the event indicators REM, REM_core, hypnagogic, other_halluc.
"""
import glob

import pandas as pd

from config import (AGE_BANDS, AGE_LABELS, FAERS_DIR, FAERS_PT, JADER_DIR, JADER_FIRST_FISCAL_YEAR, JADER_PT)

# ---------------------------------------------------------------- JADER (Japanese generic names)
JADER_CLASS = {
    "DORA": r"スボレキサント|レンボレキサント|ダリドレキサント|ボルノレキサント",
    # midazolam/remimazolam (procedural sedation) and two misspellings found in the name check are excluded
    "BZZ": r"ゾルピデム|ゾピクロン|(?<!ミダ)(?<!ミタ)(?<!レミマ)(?<!シタロ)ゾラム|ゼパム|リルマザホン",
    "RAM": r"ラメルテオン",
}
JADER_CONC = {
    "antipsychotic": r"アリピプラゾール|ブレクスピプラゾール|クエチアピン|オランザピン|リスペリドン|パリペリドン|ハロペリドール|"
                     r"クロザピン|ルラシドン|アセナピン|ブロナンセリン|ペロスピロン|クロルプロマジン|レボメプロマジン|スルピリド|"
                     r"チアプリド|ゾテピン",
    "antidepressant": r"パロキセチン|セルトラリン|フルボキサミン|エスシタロプラム|デュロキセチン|ミルナシプラン|ベンラファキシン|"
                      r"アミトリプチリン|イミプラミン|クロミプラミン|ノルトリプチリン|アモキサピン|ミルタザピン|トラゾドン|ボルチオキセチン",
    "dopaminergic": r"レボドパ|カルビドパ|プラミペキソール|ロピニロール|ロチゴチン|アマンタジン|セレギリン|ラサギリン|"
                    r"サフィナミド|エンタカポン|オピカポン|アポモルヒネ",
    "beta_blocker": r"プロプラノロール|メトプロロール|ピンドロール",
    "varenicline": r"バレニクリン",
    "montelukast": r"モンテルカスト",
    "opioid": r"オキシコドン|(?<!アポ)モルヒネ|フェンタニル|ヒドロモルフォン|トラマドール|タペンタドール|メサドン|"
              r"ブプレノルフィン|コデイン",
    "gabapentinoid": r"プレガバリン|ガバペンチン|ミロガバリン",
}


def _jader_age(a):
    a = str(a)
    if a.endswith("歳代") and a[:-2].isdigit():  # "60歳代" = 60-69 years
        n = int(a[:-2])
        return "<20" if n < 20 else ("20-39" if n < 40 else ("40-59" if n < 60 else ("60-79" if n < 80 else "80+")))
    if a in ("10歳未満", "新生児", "小児", "乳幼児", "胎児", "青少年") or "トライメスター" in a:
        return "<20"
    return "60-79" if a == "高齢者" else "unknown"  # 高齢者 = "elderly"


def build_jader():
    """Returns (df, drug). Restricted to reports from fiscal year 2014 onwards."""
    reac = pd.read_parquet(JADER_DIR / "reac.parquet")
    drug = pd.read_parquet(JADER_DIR / "drug.parquet")
    demo = pd.read_parquet(JADER_DIR / "demo.parquet").drop_duplicates("識別番号").set_index("識別番号")
    demo = demo[demo["報告年度・四半期"].str[:4].astype(int) >= JADER_FIRST_FISCAL_YEAR]
    gn = drug["医薬品（一般名）"].fillna("")
    is_sus = drug["医薬品の関与"] == "被疑薬"  # suspect drug

    df = pd.DataFrame(index=demo.index)
    for k, p in JADER_CLASS.items():
        m = gn.str.contains(p, regex=True)
        df[k] = df.index.isin(drug.loc[m & is_sus, "識別番号"]).astype(int)
        df[f"conc_{k}_nonsuspect"] = df.index.isin(drug.loc[m & ~is_sus, "識別番号"]).astype(int)
    for k, p in JADER_CONC.items():
        df[f"conc_{k}"] = df.index.isin(drug.loc[gn.str.contains(p, regex=True), "識別番号"]).astype(int)
    for k, v in JADER_PT.items():
        df[k] = df.index.isin(reac.loc[reac["有害事象"].isin(v), "識別番号"]).astype(int)

    df["age"] = demo["年齢"].map(_jader_age)
    df["sex"] = demo["性別"].map({"女性": "female", "男性": "male"}).fillna("unknown")
    r = demo["報告者の資格"].fillna("")
    # Several reporter qualifications can be joined by ";". Any physician/pharmacist/other HCP -> HCP.
    hcp = r.str.contains("医師|薬剤師|医療関係者")
    df["reporter"] = "unknown"
    df.loc[r.str.contains("消費者") & ~hcp, "reporter"] = "consumer"
    df.loc[hcp, "reporter"] = "HCP"
    df["year"] = demo["報告年度・四半期"].str[:4]
    df["country"] = "JP"
    return df, drug


# ---------------------------------------------------------------- FAERS
FAERS_CONC = {"antipsychotic": ["ANTIPSY"], "antidepressant": ["SSRI_SNRI", "TCA", "MIR_TRZ"],
              "dopaminergic": ["DOPA"], "beta_blocker": ["BB"], "varenicline": ["VAREN"], "montelukast": ["MONTE"],
              "opioid": ["OPIOID"], "gabapentinoid": ["GABA"]}


def load_faers():
    """Concatenate the quarterly extracts and deduplicate: latest caseversion per caseid
    (ties broken by FDA receipt date, then primaryid). Returns (demo, reac, drug, ther)."""
    q = FAERS_DIR / "quarters"
    cat = lambda t: pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(q / f"*_{t}.parquet")))],
                              ignore_index=True)
    demo, reac, drug, ther = cat("demo"), cat("reac"), cat("drug"), cat("ther")
    demo["cv"] = pd.to_numeric(demo["caseversion"], errors="coerce").fillna(0)
    demo["pid"] = pd.to_numeric(demo["primaryid"], errors="coerce")
    demo = demo.sort_values(["caseid", "cv", "fda_dt", "pid"]).drop_duplicates("caseid", keep="last")
    demo = demo.set_index("primaryid")
    reac, drug, ther = (x[x["primaryid"].isin(demo.index)] for x in (reac, drug, ther))
    # Name-matching corrections: apomorphine is not an opioid; "CELEXA" misspellings captured celecoxib
    nm = drug["prod_ai"].fillna("").str.upper()
    drug = drug[~((drug["cls"] == "OPIOID") & nm.str.contains("APOMORPHINE"))
                & ~((drug["cls"] == "SSRI_SNRI") & nm.str.contains("CELECOXIB"))]
    return demo, reac, drug, ther


def build_faers(loaded=None):
    """Returns (df, drug). `loaded` can be the output of load_faers() to avoid reading twice."""
    demo, reac, drug, _ = loaded or load_faers()
    is_sus = drug["role_cod"].isin(["PS", "SS"])
    df = pd.DataFrame(index=demo.index)
    for k in ["DORA", "BZZ", "RAM"]:
        m = drug["cls"] == k
        df[k] = df.index.isin(drug.loc[m & is_sus, "primaryid"]).astype(int)
        df[f"conc_{k}_nonsuspect"] = df.index.isin(drug.loc[m & ~is_sus, "primaryid"]).astype(int)
    for k, v in FAERS_CONC.items():
        df[f"conc_{k}"] = df.index.isin(drug.loc[drug["cls"].isin(v), "primaryid"]).astype(int)
    for k, v in FAERS_PT.items():
        df[k] = df.index.isin(reac.loc[reac["pt"].isin(v), "primaryid"]).astype(int)

    years = pd.to_numeric(demo["age"], errors="coerce") * demo["age_cod"].map(
        {"YR": 1, "DEC": 10, "MON": 1 / 12, "WK": 1 / 52, "DY": 1 / 365, "HR": 1 / 8760})
    df["age"] = pd.cut(years, AGE_BANDS, labels=AGE_LABELS).astype(str).replace("nan", "unknown")
    df["sex"] = demo["sex"].map({"F": "female", "M": "male"}).fillna("unknown")
    df["reporter"] = demo["occp_cod"].map({"MD": "HCP", "PH": "HCP", "OT": "HCP", "HP": "HCP",
                                           "CN": "consumer"}).fillna("unknown")
    df["year"] = demo["fda_dt"].str[:4]
    df["country"] = demo["occr_country"].fillna("unknown")
    return df, drug


def build(db):
    return build_jader() if db == "jader" else build_faers()
