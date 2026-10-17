"""Step 1 (FAERS): download the FDA quarterly ASCII files one by one, keep only the columns and rows needed,
save them as parquet and delete the zip (about 2.8 GB in total is downloaded; little is kept on disk).
Quarters that are already processed are skipped, so the script can be re-run after an interruption.

Output: FAERS_DIR/quarters/{quarter}_{demo,reac,drug,ther}.parquet
Usage : python src/faers_download.py
"""
import io
import re
import time
import urllib.request
import zipfile

import pandas as pd

from config import FAERS_DIR, FAERS_FIRST_QUARTER, FAERS_LAST_QUARTER, FAERS_PT

OUT = FAERS_DIR / "quarters"
TMP = FAERS_DIR / "tmp"
OUT.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)

# Drug groups. Exposure groups are DORA, BZZ and RAM; SSRI_SNRI, TCA and MIR_TRZ are used in the
# post-discontinuation analysis; the remaining groups are covariates (concomitant drugs).
CLASS = {
    "DORA": r"SUVOREXANT|LEMBOREXANT|DARIDOREXANT|BELSOMRA|DAYVIGO|QUVIVIQ",
    "BZZ": r"ZOLPIDEM|ZOPICLONE|ZALEPLON|TRIAZOLAM|BROTIZOLAM|ETIZOLAM|FLUNITRAZEPAM|NITRAZEPAM|TEMAZEPAM|"
           r"ESTAZOLAM|FLURAZEPAM|QUAZEPAM|LORMETAZEPAM|RILMAZAFONE|ALPRAZOLAM|LORAZEPAM|DIAZEPAM|CLONAZEPAM|"
           r"BROMAZEPAM|CHLORDIAZEPOXIDE|OXAZEPAM|CLORAZEPATE|CLOTIAZEPAM|CLOBAZAM|AMBIEN|LUNESTA|SONATA|HALCION|"
           r"XANAX|ATIVAN|VALIUM|KLONOPIN|RESTORIL",
    "RAM": r"RAMELTEON|ROZEREM",
    "SSRI_SNRI": r"PAROXETINE|SERTRALINE|FLUVOXAMINE|CITALOPRAM|FLUOXETINE|DULOXETINE|MILNACIPRAN|VENLAFAXINE|"
                 r"PAXIL|ZOLOFT|LEXAPRO|CELEXA|PROZAC|CYMBALTA|EFFEXOR|PRISTIQ",
    "TCA": r"AMITRIPTYLINE|IMIPRAMINE|CLOMIPRAMINE|NORTRIPTYLINE|DESIPRAMINE|DOXEPIN|AMOXAPINE|TRIMIPRAMINE|PROTRIPTYLINE",
    "MIR_TRZ": r"MIRTAZAPINE|TRAZODONE|REMERON|DESYREL",
    "ANTIPSY": r"ARIPIPRAZOLE|BREXPIPRAZOLE|CARIPRAZINE|QUETIAPINE|OLANZAPINE|RISPERIDONE|PALIPERIDONE|HALOPERIDOL|"
               r"CLOZAPINE|LURASIDONE|ZIPRASIDONE|ASENAPINE|BLONANSERIN|PEROSPIRONE|CHLORPROMAZINE|LEVOMEPROMAZINE|"
               r"SULPIRIDE|LUMATEPERONE|ILOPERIDONE|PIMAVANSERIN",
    "DOPA": r"LEVODOPA|CARBIDOPA|PRAMIPEXOLE|ROPINIROLE|ROTIGOTINE|AMANTADINE|SELEGILINE|RASAGILINE|SAFINAMIDE|"
            r"ENTACAPONE|OPICAPONE|APOMORPHINE",
    "BB": r"PROPRANOLOL|METOPROLOL|PINDOLOL",
    "VAREN": r"VARENICLINE|CHANTIX|CHAMPIX",
    "MONTE": r"MONTELUKAST|SINGULAIR",
    "OPIOID": r"OXYCODONE|MORPHINE|FENTANYL|HYDROCODONE|TRAMADOL|TAPENTADOL|HYDROMORPHONE|METHADONE|BUPRENORPHINE|CODEINE",
    "GABA": r"PREGABALIN|GABAPENTIN|MIROGABALIN|LYRICA|NEURONTIN",
}
KEEP_PT = set(FAERS_PT["REM"]) | set(FAERS_PT["hypnagogic"]) | set(FAERS_PT["other_halluc"])
WD_RE = r"withdrawal|rebound|discontinuation syndrome"


def quarters():
    y, q = FAERS_FIRST_QUARTER
    while (y, q) <= FAERS_LAST_QUARTER:
        yield y, q
        y, q = (y + 1, 1) if q == 4 else (y, q + 1)


def read(zf, table, cols):
    # e.g. ascii/DEMO18Q1_new.txt (2018Q1) or ASCII/demo14q3.txt
    name = next(n for n in zf.namelist() if re.search(rf"(^|/){table}\d\dq\d(_new)?\.txt$", n, re.I))
    with zf.open(name) as f:
        df = pd.read_csv(io.TextIOWrapper(f, encoding="latin-1"), sep="$", dtype=str, quoting=3,
                         on_bad_lines="skip", low_memory=False)
    df.columns = [c.strip().lower() for c in df.columns]
    return df[[c for c in cols if c in df.columns]]


def download(y, q, zp):
    # The FDA site uses both "2025q3" and "2025Q4" style names, so try both.
    for tag in (f"{y}q{q}", f"{y}Q{q}"):
        url = f"https://fis.fda.gov/content/Exports/faers_ascii_{tag}.zip"
        for _ in range(3):
            try:
                urllib.request.urlretrieve(url, zp)
                if zipfile.is_zipfile(zp):
                    return
            except Exception:
                pass
            time.sleep(5)
    raise RuntimeError(f"download failed for {y}Q{q}")


def process(y, q):
    tag = f"{y}q{q}"
    if (OUT / f"{tag}_demo.parquet").exists():
        return "skip"
    zp = TMP / f"{tag}.zip"
    download(y, q, zp)
    with zipfile.ZipFile(zp) as zf:
        demo = read(zf, "demo", ["primaryid", "caseid", "caseversion", "fda_dt", "event_dt", "age", "age_cod", "sex",
                                 "occp_cod", "occr_country"])
        reac = read(zf, "reac", ["primaryid", "pt"])
        reac["pt"] = reac["pt"].str.strip().str.lower()
        reac = reac[reac["pt"].isin(KEEP_PT) | reac["pt"].str.contains(WD_RE, na=False)]
        drug = read(zf, "drug", ["primaryid", "drug_seq", "role_cod", "drugname", "prod_ai"])
        txt = (drug["prod_ai"].fillna("") + " " + drug["drugname"].fillna("")).str.upper()
        parts = [drug.loc[txt.str.contains(p, regex=True), ["primaryid", "drug_seq", "role_cod", "prod_ai", "drugname"]]
                 .assign(cls=k) for k, p in CLASS.items()]
        drug = pd.concat(parts, ignore_index=True)
        ther = read(zf, "ther", ["primaryid", "dsg_drug_seq", "start_dt", "end_dt"])
        ther = ther.merge(drug.loc[drug["role_cod"].isin(["PS", "SS"]), ["primaryid", "drug_seq"]].drop_duplicates(),
                          left_on=["primaryid", "dsg_drug_seq"], right_on=["primaryid", "drug_seq"])
    reac.to_parquet(OUT / f"{tag}_reac.parquet")
    drug.to_parquet(OUT / f"{tag}_drug.parquet")
    ther.to_parquet(OUT / f"{tag}_ther.parquet")
    demo.to_parquet(OUT / f"{tag}_demo.parquet")  # written last: marks the quarter as complete
    zp.unlink()
    return f"demo {len(demo):,} reac {len(reac):,} drug {len(drug):,} ther {len(ther):,}"


if __name__ == "__main__":
    for y, q in quarters():
        t = time.time()
        try:
            msg = process(y, q)
        except Exception as e:  # one failed quarter should not stop the others; re-run to retry
            msg = f"ERROR {e!r}"
        print(f"{y}Q{q} {msg} ({time.time() - t:.0f}s)", flush=True)
