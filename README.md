# orexin-rem-pv

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23447366.svg)](https://doi.org/10.5281/zenodo.23447366)

Analysis code for:

> **Orexin receptor antagonists and REM sleep–related experiences at the sleep–wake boundary: a disproportionality analysis of individual case safety reports in JADER and FAERS** (preprint, [link to be added])

This project is a disproportionality analysis of two spontaneous reporting databases:
- the Japanese Adverse Drug Event Report database (**JADER**, PMDA);
- the FDA Adverse Event Reporting System (**FAERS**).

Reports with a suspected dual orexin receptor antagonist (DORA: suvorexant, lemborexant, daridorexant) are compared with reports with a suspected benzodiazepine or Z-drug (BZ/Z). The events are grouped into three tiers by their proximity to REM sleep:

1. REM-core: nightmare, abnormal dreams, sleep paralysis
2. Hypnagogic/hypnopompic hallucinations
3. Other hallucinations

Reporting odds ratios are estimated with Firth-penalised logistic regression. A trend across the three tiers is tested, and FAERS reports with onset after the end of therapy are described. The reporting follows the READUS-PV guideline.

**No data are included in this repository.** Both databases are public and must be obtained from their custodians (see below).

## Requirements

- Python 3.12 (tested with 3.12.10)
- Packages: see `requirements.txt` (tested with pandas 3.0.6, NumPy 2.4.4, SciPy 1.18.1, statsmodels 0.15.0, pyarrow)
- About 3 GB of free disk space while FAERS is being downloaded, and 16 GB of RAM for the FAERS steps
- Internet access for the FAERS download (about 2.8 GB in 48 quarterly files)

```bash
pip install -r requirements.txt
```

## Data

### JADER (manual download)
1. Open the PMDA page "医薬品副作用データベース" (https://www.pmda.go.jp/safety/info-services/drugs/adr-info/suspected-adr/0003.html), read and accept the terms of use, and download the full dataset (one zip, about 140 MB). An image CAPTCHA has to be entered by hand.
2. Extract `demoYYYYMM.csv`, `drugYYYYMM.csv` and `reacYYYYMM.csv` into `data/jader/`.
3. If the release is not `202610`, set `JADER_RELEASE` (for example `JADER_RELEASE=202701`).

Under the JADER terms of use, any publication must state that the results were obtained using the PMDA Japanese Adverse Drug Event Report database (JADER).

### FAERS (automatic download)
`src/faers_download.py` downloads the quarterly ASCII files for 2014Q3–2026Q2 from https://fis.fda.gov. To save disk space, it keeps only the needed columns, deletes each zip after use, and skips quarters that have already been processed.

### Paths
The default folders are `data/jader`, `data/faers` and `output`. You can change them with environment variables: `JADER_DIR`, `JADER_RELEASE`, `FAERS_DIR` and `OUT_DIR` (see `src/config.py`).

## How to run

Run all commands from the repository root.

```bash
# 1. Prepare data
python src/jader_prepare.py                 # JADER CSV -> parquet (latest version of each case)
python src/faers_download.py                # FAERS quarterly files -> parquet (about 1-1.5 h)

# 2. Primary analyses (Tables 2 and 3)
python src/main_models.py jader
python src/main_models.py faers

# 3. Sensitivity analyses and name-matching check (Table 4, Supplementary Table S1)
python src/sensitivity.py jader
python src/sensitivity.py faers

# 4. Report characteristics (Table 1)
python src/table1.py jader
python src/table1.py faers

# 5. Onset after the end of therapy, FAERS only (Table 5)
python src/faers_post_discontinuation.py
```

Results are written to `output/` as text and CSV files.

## Files

| File | Purpose |
|---|---|
| `src/config.py` | Paths, study period, age bands, event definitions (MedDRA PTs) |
| `src/jader_prepare.py` | Reads JADER CSVs and keeps the latest report version of each case |
| `src/faers_download.py` | Downloads and extracts FAERS; drug-group patterns for FAERS |
| `src/build_data.py` | Builds the analysis table (one row per report); drug-group patterns for JADER; FAERS deduplication |
| `src/firth.py` | Firth logistic regression on aggregated data, trend test, ratio of RORs |
| `src/main_models.py` | Model A (all other reports as reference) and Model B (active comparator), trend test |
| `src/sensitivity.py` | Sensitivity analyses (reporter type, reporting period, country) and name-matching check |
| `src/table1.py` | Characteristics of DORA-only and BZ/Z-only reports, report counts |
| `src/faers_post_discontinuation.py` | Onset after the end of therapy and co-reported withdrawal/rebound terms |

## Methods in brief

- **Deduplication**
  - JADER: latest report version (報告回数) for each case ID.
  - FAERS: latest `caseversion` for each `caseid`; ties are broken by `fda_dt`, then by `primaryid`.
- **Exposure**: suspect drugs only (JADER 被疑薬; FAERS role PS/SS). Drugs are matched by regular expressions on generic names (JADER) or on `prod_ai` + `drugname` (FAERS). Midazolam and remimazolam are excluded from BZ/Z. Apomorphine is excluded from opioids, and celecoxib (captured through a misspelled brand name) from SSRIs.
- **Covariates**: age band, sex, reporter type, and indicators for antipsychotics, antidepressants, dopaminergic drugs, lipophilic β-blockers, varenicline, montelukast, opioids, gabapentinoids, and non-suspect DORA/BZ/Z/ramelteon.
- **Models**: Firth-penalised logistic regression fitted to aggregated strata. Wald 95% CIs are reported. A signal is defined as a lower confidence limit above 1.
- **Trend test**: inverse-variance weighted regression of the log RORs on the tier rank (0, 1, 2).

Spontaneous reports cannot show causality or incidence; the results are hypothesis-generating.

## Notes

- Comments inside the code are partly in Japanese. JADER field names are in Japanese because the source files are.
- The FAERS download URLs and the JADER file layout may change. If a FAERS quarter fails, re-run the script; finished quarters are skipped.
- Running the analysis with a newer data release will give slightly different numbers.

## License

MIT (see `LICENSE`).

## Citation

Code (all versions): Ozaki T. orexin-rem-pv. Zenodo. https://doi.org/10.5281/zenodo.23447366 (v1.0.0: https://doi.org/10.5281/zenodo.23447367)

Article:

Ozaki T. Orexin receptor antagonists and REM sleep–related experiences at the sleep–wake boundary: a disproportionality analysis of JADER and FAERS. medRxiv [year]. doi: [to be added]
