#!/usr/bin/env python3
"""Build pheno/cov files from two projection sscore files.

Outputs are named by semantic labels (not raw source file prefixes).

Primary sscore (A) outputs:
1) <label_a>.pheno.tsv
2) <label_a>.cov.sex.tsv
3) <label_a>.cov.sex_age_agez.tsv
4) <label_a>.age_na.fid_iid

Secondary sscore (B) outputs:
1) <label_b>.cov.sex.tsv
2) <label_b>.cov.sex_age_agez.tsv
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Optional, Tuple

import pandas as pd


@dataclass
class BuildResult:
    label: str
    sscore_path: str
    pheno_path: Optional[str]
    cov_sex_path: str
    cov_agez_path: str
    n_total: int
    n_age_na: int
    n_age_kept: int


def _normalize_series_as_str(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip()


def _read_sscore(sscore_path: str) -> pd.DataFrame:
    df = pd.read_csv(sscore_path, sep=r"\s+", engine="python")
    if "#FID" not in df.columns and "FID" in df.columns:
        df = df.rename(columns={"FID": "#FID"})

    required = ["#FID", "IID", "PHENO1"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{sscore_path}: missing required columns: {missing}")

    pc_cols = [c for c in df.columns if c.startswith("PC")]
    if not pc_cols:
        raise ValueError(f"{sscore_path}: no PC columns found (expected columns starting with 'PC')")

    return df


def _build_lookup(sample_info_path: str, sample_id_col: str, sex_col: str, age_col: str) -> pd.DataFrame:
    sample_df = pd.read_excel(sample_info_path)
    required = [sample_id_col, sex_col, age_col]
    missing = [c for c in required if c not in sample_df.columns]
    if missing:
        raise ValueError(f"sample_info missing required columns: {missing}")

    sample_df = sample_df[[sample_id_col, sex_col, age_col]].copy()
    sample_df[sample_id_col] = _normalize_series_as_str(sample_df[sample_id_col])
    sample_df = sample_df.drop_duplicates(subset=[sample_id_col], keep="first")
    return sample_df


def _map_sex(raw: pd.Series, female_value: str, male_value: str) -> pd.Series:
    female_norm = str(female_value).strip().lower()
    male_norm = str(male_value).strip().lower()

    raw_norm = raw.astype(str).str.strip().str.lower()
    mapped = pd.Series(pd.NA, index=raw.index, dtype="Int64")
    mapped.loc[raw_norm == male_norm] = 1
    mapped.loc[raw_norm == female_norm] = 2
    return mapped


def _zscore(s: pd.Series) -> pd.Series:
    mean = s.mean()
    std = s.std(ddof=0)
    if pd.isna(std) or std == 0:
        return pd.Series([0.0] * len(s), index=s.index)
    return (s - mean) / std


def _build_one(
    output_label: str,
    sscore_path: str,
    lookup_df: pd.DataFrame,
    sample_id_col: str,
    sex_col: str,
    age_col: str,
    female_value: str,
    male_value: str,
    write_pheno_file: bool,
    write_age_na_file: bool,
) -> Tuple[BuildResult, str | None]:
    df = _read_sscore(sscore_path).copy()
    df["IID"] = _normalize_series_as_str(df["IID"])
    df["#FID"] = _normalize_series_as_str(df["#FID"])

    pc_cols = [c for c in df.columns if c.startswith("PC")]

    merged = df.merge(
        lookup_df,
        left_on="IID",
        right_on=sample_id_col,
        how="left",
        validate="m:1",
    )

    merged["SEX"] = _map_sex(merged[sex_col], female_value=female_value, male_value=male_value)
    if merged["SEX"].isna().any():
        n_missing = int(merged["SEX"].isna().sum())
        raise ValueError(
            f"{sscore_path}: {n_missing} samples cannot map SEX from sample_info "
            f"using female='{female_value}' and male='{male_value}'"
        )

    merged["AGE"] = pd.to_numeric(merged[age_col], errors="coerce")

    pheno_path = None
    if write_pheno_file:
        pheno_path = f"{output_label}.pheno.tsv"
        pheno_df = merged[["#FID", "IID", "PHENO1"]].copy()
        pheno_df.to_csv(pheno_path, sep="\t", index=False)

    cov_sex_path = f"{output_label}.cov.sex.tsv"
    cov_agez_path = f"{output_label}.cov.sex_age_agez.tsv"

    cov_sex_df = merged[["#FID", "IID", "SEX"] + pc_cols].copy()
    cov_sex_df["SEX"] = cov_sex_df["SEX"].astype(int)
    cov_sex_df.to_csv(cov_sex_path, sep="\t", index=False)

    age_na_df = merged.loc[merged["AGE"].isna(), ["#FID", "IID"]].copy()

    cov_age_df = merged.loc[~merged["AGE"].isna(), ["#FID", "IID", "SEX", "AGE"] + pc_cols].copy()
    cov_age_df["SEX"] = cov_age_df["SEX"].astype(int)
    cov_age_df["AGE_Z"] = _zscore(cov_age_df["AGE"].astype(float))

    final_cols = ["#FID", "IID", "SEX", "AGE", "AGE_Z"] + pc_cols
    cov_age_df = cov_age_df[final_cols]
    cov_age_df.to_csv(cov_agez_path, sep="\t", index=False)

    age_na_path = None
    if write_age_na_file:
        age_na_path = f"{output_label}.age_na.fid_iid"
        age_na_df.to_csv(age_na_path, sep="\t", index=False, header=False)

    result = BuildResult(
        label=output_label,
        sscore_path=sscore_path,
        pheno_path=pheno_path,
        cov_sex_path=cov_sex_path,
        cov_agez_path=cov_agez_path,
        n_total=int(len(merged)),
        n_age_na=int(len(age_na_df)),
        n_age_kept=int(len(cov_age_df)),
    )
    return result, age_na_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build cov/pheno files from two sscore files")
    parser.add_argument("--sscore-a", required=True)
    parser.add_argument("--sscore-b", required=True)
    parser.add_argument("--label-a", required=True)
    parser.add_argument("--label-b", required=True)
    parser.add_argument("--sample-info", required=True)
    parser.add_argument("--sample-id-col", required=True)
    parser.add_argument("--sex-col", required=True)
    parser.add_argument("--sex-female-value", required=True)
    parser.add_argument("--sex-male-value", required=True)
    parser.add_argument("--age-col", required=True)
    parser.add_argument("--out-log", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    lookup_df = _build_lookup(
        sample_info_path=args.sample_info,
        sample_id_col=args.sample_id_col,
        sex_col=args.sex_col,
        age_col=args.age_col,
    )

    result_a, age_na_path = _build_one(
        output_label=args.label_a,
        sscore_path=args.sscore_a,
        lookup_df=lookup_df,
        sample_id_col=args.sample_id_col,
        sex_col=args.sex_col,
        age_col=args.age_col,
        female_value=args.sex_female_value,
        male_value=args.sex_male_value,
        write_pheno_file=True,
        write_age_na_file=True,
    )

    result_b, _ = _build_one(
        output_label=args.label_b,
        sscore_path=args.sscore_b,
        lookup_df=lookup_df,
        sample_id_col=args.sample_id_col,
        sex_col=args.sex_col,
        age_col=args.age_col,
        female_value=args.sex_female_value,
        male_value=args.sex_male_value,
        write_pheno_file=False,
        write_age_na_file=False,
    )

    with open(args.out_log, "w", encoding="utf-8") as f:
        f.write("Cov/Pheno build summary\n")
        f.write(f"SCORE_A_LABEL: {result_a.label}\n")
        f.write(f"SCORE_A_PATH: {result_a.sscore_path}\n")
        if result_a.pheno_path is not None:
            f.write(f"  PHENO: {result_a.pheno_path}\n")
        f.write(f"  COV_SEX: {result_a.cov_sex_path}\n")
        f.write(f"  COV_SEX_AGE_AGEZ: {result_a.cov_agez_path}\n")
        f.write(f"  N_TOTAL: {result_a.n_total}\n")
        f.write(f"  N_AGE_NA: {result_a.n_age_na}\n")
        f.write(f"  N_AGE_KEPT: {result_a.n_age_kept}\n")
        f.write("\n")
        f.write(f"SCORE_B_LABEL: {result_b.label}\n")
        f.write(f"SCORE_B_PATH: {result_b.sscore_path}\n")
        f.write("  PHENO: not generated (primary score only)\n")
        f.write(f"  COV_SEX: {result_b.cov_sex_path}\n")
        f.write(f"  COV_SEX_AGE_AGEZ: {result_b.cov_agez_path}\n")
        f.write(f"  N_TOTAL: {result_b.n_total}\n")
        f.write(f"  N_AGE_NA: {result_b.n_age_na}\n")
        f.write(f"  N_AGE_KEPT: {result_b.n_age_kept}\n")
        f.write("\n")
        if age_na_path is not None:
            f.write(f"AGE_NA_FID_IID_FILE: {age_na_path}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
