"""Paths and constants shared by all scripts.

Data are NOT included in this repository. Override the defaults with environment variables:
  JADER_DIR      folder with the PMDA JADER CSV files (demoYYYYMM.csv, drugYYYYMM.csv, reacYYYYMM.csv)
  JADER_RELEASE  release suffix of the JADER files (default 202610)
  FAERS_DIR      working folder for the FAERS quarterly extracts
  OUT_DIR        folder for result files
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
JADER_DIR = Path(os.environ.get("JADER_DIR", ROOT / "data" / "jader"))
JADER_RELEASE = os.environ.get("JADER_RELEASE", "202610")
FAERS_DIR = Path(os.environ.get("FAERS_DIR", ROOT / "data" / "faers"))
OUT_DIR = Path(os.environ.get("OUT_DIR", ROOT / "output"))
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Study period
JADER_FIRST_FISCAL_YEAR = 2014           # DORA era (suvorexant launched in Japan in 2014)
FAERS_FIRST_QUARTER = (2014, 3)
FAERS_LAST_QUARTER = (2026, 2)

# Age bands (years) used for adjustment; reference level is 40-59
AGE_BANDS = [-1, 19.999, 39.999, 59.999, 79.999, 200]
AGE_LABELS = ["<20", "20-39", "40-59", "60-79", "80+"]
REF_LEVELS = {"age": "40-59", "sex": "female", "reporter": "HCP"}

# Event tiers ordered by proximity to REM sleep
FAERS_PT = {
    "REM": ["nightmare", "abnormal dreams", "sleep paralysis", "hypnagogic hallucination"],
    "REM_core": ["nightmare", "abnormal dreams", "sleep paralysis"],
    "hypnagogic": ["hypnagogic hallucination", "hypnopompic hallucination"],
    "other_halluc": ["hallucination", "hallucination, auditory", "hallucination, visual", "hallucination, tactile",
                     "hallucination, olfactory", "hallucination, gustatory", "hallucinations, mixed",
                     "somatic hallucination"],
}
JADER_PT = {
    "REM": ["悪夢", "異常な夢", "睡眠時麻痺", "入眠時幻覚"],
    "REM_core": ["悪夢", "異常な夢", "睡眠時麻痺"],
    "hypnagogic": ["入眠時幻覚", "覚醒時幻覚"],
    "other_halluc": ["幻覚", "幻聴", "幻視", "体感幻覚", "混合性幻覚", "幻嗅", "触覚性幻覚"],
}
LADDER = ["REM_core", "hypnagogic", "other_halluc"]
OUTCOMES = ["REM"] + LADDER
EXPOSURES = ["DORA", "BZZ", "RAM"]
