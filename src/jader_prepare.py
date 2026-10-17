"""Step 1 (JADER): read the PMDA CSV files, keep the latest report version of each case, save as parquet.

Input : JADER_DIR/demo{release}.csv, drug{release}.csv, reac{release}.csv (Shift_JIS / cp932)
Output: JADER_DIR/demo.parquet, drug.parquet, reac.parquet
Usage : python src/jader_prepare.py
"""
import pandas as pd

from config import JADER_DIR, JADER_RELEASE


def read(name, cols):
    return pd.read_csv(JADER_DIR / f"{name}{JADER_RELEASE}.csv", encoding="cp932", dtype=str, usecols=cols,
                       encoding_errors="replace")


demo = read("demo", ["識別番号", "報告回数", "性別", "年齢", "報告年度・四半期", "報告者の資格"])
drug = read("drug", ["識別番号", "報告回数", "医薬品の関与", "医薬品（一般名）", "投与開始日", "投与終了日", "医薬品の処置"])
reac = read("reac", ["識別番号", "報告回数", "有害事象", "転帰", "有害事象の発現日"])

# Keep only the latest report version (報告回数) of each case (識別番号)
last = demo.groupby("識別番号")["報告回数"].max().reset_index()
demo = demo.merge(last, on=["識別番号", "報告回数"])
reac = reac.merge(last, on=["識別番号", "報告回数"]).drop_duplicates(["識別番号", "有害事象"])
drug = drug.merge(last, on=["識別番号", "報告回数"])

for name, df in [("demo", demo), ("reac", reac), ("drug", drug)]:
    df.to_parquet(JADER_DIR / f"{name}.parquet")
print(f"cases {demo['識別番号'].nunique():,}  reaction rows {len(reac):,}  drug rows {len(drug):,}")
