"""Core metric computation utilities for variant-level QC.

This module runs PLINK2-based metric generation and produces a consolidated
`*.variant_qc_summary.tsv` table used by downstream filtering steps.

Scope:
- Compute whole-cohort and subgroup-level AAF/MAF, VMISS, and HWE statistics.
- Support depth-based groups (30X/15X) using sample metadata.
- Produce high-throughput summary output with chunked multiprocessing.
"""

import subprocess
import os
import tempfile
import uuid
import logging
import time
from typing import Optional
import concurrent.futures
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('default')

LOGGER = logging.getLogger("VariantQC.Metrics")


def _log_info(phase: str, message: str) -> None:
    LOGGER.info(f"[{phase}] {message}")


def _log_warn(phase: str, message: str) -> None:
    LOGGER.warning(f"[{phase}] {message}")


def _log_error(phase: str, message: str) -> None:
    LOGGER.error(f"[{phase}] {message}")


def _build_summary_chunk_worker(
    idx,
    tmpdir,
    all_af_chunk,
    all_vmiss_chunk,
    case_af_chunk,
    ctrl_af_chunk,
    case_vmiss_chunk,
    ctrl_vmiss_chunk,
    case_hwe_chunk,
    ctrl_hwe_chunk,
    dp30_vmiss_chunk,
    dp15_vmiss_chunk,
    vmiss_mode=None,
    thr_vmiss=None,
    thr_case=None,
    thr_ctrl=None,
    thr_dp30=None,
    thr_dp15=None,
):
    """Worker to build one summary chunk and write it to a temp TSV.

    This is used by the streaming multi-process implementation of
    run_plink2_variant_qc. It assumes that all input DataFrames are
    already aligned by variant ID and correspond to the same slice of
    the PLINK outputs.
    """

    base_ids = all_af_chunk["ID"].astype(str)

    aaf_all = all_af_chunk["ALT_FREQS"].astype(float)
    maf_all = aaf_all.where(aaf_all <= 0.5, 1.0 - aaf_all)

    vmiss_overall = all_vmiss_chunk["F_MISS"].astype(float)
    case_vmiss_vals = case_vmiss_chunk["F_MISS"].astype(float)
    ctrl_vmiss_vals = ctrl_vmiss_chunk["F_MISS"].astype(float)

    case_aaf = case_af_chunk["ALT_FREQS"].astype(float)
    ctrl_aaf = ctrl_af_chunk["ALT_FREQS"].astype(float)
    case_maf = case_aaf.where(case_aaf <= 0.5, 1.0 - case_aaf)
    ctrl_maf = ctrl_aaf.where(ctrl_aaf <= 0.5, 1.0 - ctrl_aaf)

    case_hwe_p = case_hwe_chunk["P"].astype(float)
    ctrl_hwe_p = ctrl_hwe_chunk["P"].astype(float)

    if dp30_vmiss_chunk is not None:
        dp30_vmiss_vals = dp30_vmiss_chunk["F_MISS"].astype(float)
    else:
        dp30_vmiss_vals = pd.Series([float("nan")] * len(all_af_chunk))

    if dp15_vmiss_chunk is not None:
        dp15_vmiss_vals = dp15_vmiss_chunk["F_MISS"].astype(float)
    else:
        dp15_vmiss_vals = pd.Series([float("nan")] * len(all_af_chunk))

    summary_df = pd.DataFrame({
        "VARIANT_ID": base_ids,
        "MAF": maf_all,
        "VMISS": vmiss_overall,
        "CASE_VMISS": case_vmiss_vals,
        "CTRL_VMISS": ctrl_vmiss_vals,
        "30X_VMISS": dp30_vmiss_vals,
        "15X_VMISS": dp15_vmiss_vals,
        "CASE_AAF": case_aaf,
        "CASE_MAF": case_maf,
        "CTRL_AAF": ctrl_aaf,
        "CTRL_MAF": ctrl_maf,
        "CASE_HWE": case_hwe_p,
        "CTRL_HWE": ctrl_hwe_p,
    })

    # Optionally encode VMISS pass/fail directly into the summary table.
    # This avoids re-loading a large VMISS pass list for downstream HWE QC.
    if vmiss_mode is not None:
        pass_mask = None
        if vmiss_mode == "mix" and thr_vmiss is not None:
            thr = float(thr_vmiss)
            pass_mask = vmiss_overall.astype(float) < thr
        elif vmiss_mode == "dp" and thr_dp30 is not None and thr_dp15 is not None:
            thr30 = float(thr_dp30)
            thr15 = float(thr_dp15)
            pass_mask = (dp30_vmiss_vals.astype(float) <= thr30) & (dp15_vmiss_vals.astype(float) <= thr15)
        elif vmiss_mode == "case_ctrl" and thr_ctrl is not None and thr_case is not None:
            thr_ctrl_f = float(thr_ctrl)
            thr_case_f = float(thr_case)
            pass_mask = (ctrl_vmiss_vals.astype(float) <= thr_ctrl_f) & (case_vmiss_vals.astype(float) <= thr_case_f)

        if pass_mask is not None:
            summary_df["PASS_VMISS"] = pd.Series(pass_mask, index=summary_df.index).astype(bool)

    tmp_output = os.path.join(tmpdir, f"summary_chunk_{idx}_{uuid.uuid4().hex}.tsv")
    summary_df.to_csv(tmp_output, sep="\t", header=False, index=False, lineterminator="\n")

    return tmp_output, len(summary_df), idx

def run_plink2_variant_qc(
    bed_prefix: str,
    tmpdir: str = "/tmp/variant_qc",
    plink2_path: str = "plink2",
    threads: int = 8,
    output_prefix: str = "cteph_agp3k",
    verbose: bool = True,
    info_path: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx",
    sample_col: str = "ID",
    target_dp_col: str = "Target DP (JHRPv4)",
    pheno_col: str = "STATUS",
    case_value: str = "CASE",
    ctrl_value: str = "CTRL",
    pihat_vertex_cover_tsv_for_hwe: Optional[str] = None,
    summary_workers: int = 4,
    summary_chunksize: int = 100000,
    vmiss_mode: Optional[str] = None,
    vmiss_json_path: Optional[str] = None,
) -> str:
    """Run PLINK2-driven variant QC and write `<output_prefix>.variant_qc_summary.tsv`."""

    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()
    else:
        os.makedirs(tmpdir, exist_ok=True)

    def run_cmd(cmd, desc: Optional[str] = None):
        if verbose and desc:
            _log_info("PHASE-1", f"Running PLINK2 step: {desc}")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            _log_error("PHASE-1", f"Command failed: {' '.join(cmd)}")
            raise e

    out_all = os.path.join(tmpdir, "all_samples")
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--threads", str(threads),
        "--freq", "--missing", "--out", out_all
    ], "All-sample AAF + VMISS")

    _log_info(
        "PHASE-1",
        f"Starting run_plink2_variant_qc: bed_prefix={bed_prefix}, threads={threads}, tmpdir={tmpdir}"
    )
    try:
        info_df = pd.read_excel(info_path)
    except Exception as e:
        raise RuntimeError(f"Failed to read info_path: {info_path}. Reason: {e}")

    if sample_col not in info_df.columns or target_dp_col not in info_df.columns:
        raise ValueError(f"Missing required columns in info_path: {sample_col} or {target_dp_col}")

    sample_missing_mask = info_df[sample_col].isna() | info_df[sample_col].astype(str).str.strip().eq("")
    sample_missing_n = int(sample_missing_mask.sum())
    if sample_missing_n > 0:
        _log_warn("PHASE-1", f"Dropping {sample_missing_n} metadata rows with missing/blank {sample_col}")
        info_df = info_df.loc[~sample_missing_mask].copy()
    _log_info("PHASE-1", f"Metadata rows after {sample_col} cleanup: {len(info_df):,}")

    fam_df = pd.read_csv(f"{bed_prefix}.fam", sep=r"\s+", header=None)
    fam_df.columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
    fam_df["IID_str"] = fam_df["IID"].astype(str)

    info_df[sample_col] = info_df[sample_col].astype(str).str.strip()
    _log_info("PHASE-1", f"Loaded metadata: path={info_path}, rows={len(info_df):,}")
    case_samples = info_df[info_df[pheno_col] == case_value][sample_col]
    ctrl_samples = info_df[info_df[pheno_col] == ctrl_value][sample_col]

    case_iids = fam_df[fam_df["IID_str"].isin(case_samples)][["FID", "IID"]]
    ctrl_iids = fam_df[fam_df["IID_str"].isin(ctrl_samples)][["FID", "IID"]]
    matched_unique = pd.concat([case_iids, ctrl_iids], ignore_index=True).drop_duplicates().shape[0]
    _log_info(
        "PHASE-1",
        "FAM matching summary: "
        f"case_matched={len(case_iids):,}, "
        f"ctrl_matched={len(ctrl_iids):,}, "
        f"unique_matched={matched_unique:,}, "
        f"fam_total={len(fam_df):,}"
    )

    case_iids_vmiss = case_iids.copy()
    ctrl_iids_vmiss = ctrl_iids.copy()
    case_iids_hwe = case_iids.copy()
    ctrl_iids_hwe = ctrl_iids.copy()

    if pihat_vertex_cover_tsv_for_hwe:
        if not os.path.exists(pihat_vertex_cover_tsv_for_hwe):
            _log_warn("PHASE-1", f"PIHAT vertex-cover file not found, skip HWE exclusion: {pihat_vertex_cover_tsv_for_hwe}")
        else:
            try:
                pihat_df = pd.read_csv(pihat_vertex_cover_tsv_for_hwe, sep="\t")
                required_cols = {"IID", "SELECTED_FOR_REMOVAL"}
                if not required_cols.issubset(set(pihat_df.columns)):
                    _log_warn("PHASE-1", f"PIHAT vertex-cover file missing required columns {required_cols}, skip HWE exclusion")
                else:
                    removal_mask = pihat_df["SELECTED_FOR_REMOVAL"].astype(str).str.strip().str.lower().isin(["true", "1", "yes", "y"])
                    removal_iids = set(pihat_df.loc[removal_mask, "IID"].astype(str))

                    if len(removal_iids) > 0:
                        case_before = len(case_iids_hwe)
                        ctrl_before = len(ctrl_iids_hwe)
                        case_iids_hwe = case_iids_hwe[~case_iids_hwe["IID"].astype(str).isin(removal_iids)]
                        ctrl_iids_hwe = ctrl_iids_hwe[~ctrl_iids_hwe["IID"].astype(str).isin(removal_iids)]
                        _log_info(
                            "PHASE-1",
                            "Applied PIHAT-based exclusion for HWE metrics only: "
                            f"remove_iids={len(removal_iids)}, "
                            f"case_hwe_kept={len(case_iids_hwe)}/{case_before}, "
                            f"ctrl_hwe_kept={len(ctrl_iids_hwe)}/{ctrl_before}"
                        )
                    else:
                        _log_info("PHASE-1", "No samples with SELECTED_FOR_REMOVAL=true in PIHAT file; use all samples for HWE")
            except Exception as e:
                _log_warn("PHASE-1", f"Failed to parse PIHAT vertex-cover file, skip HWE exclusion: {e}")

    case_iid_path_vmiss = os.path.join(tmpdir, "case_iids_vmiss.txt")
    ctrl_iid_path_vmiss = os.path.join(tmpdir, "ctrl_iids_vmiss.txt")
    case_iid_path_hwe = os.path.join(tmpdir, "case_iids_hwe.txt")
    ctrl_iid_path_hwe = os.path.join(tmpdir, "ctrl_iids_hwe.txt")
    _log_info(
        "PHASE-1",
        "IID counts for downstream metrics: "
        f"case_vmiss={len(case_iids_vmiss):,}, ctrl_vmiss={len(ctrl_iids_vmiss):,}, "
        f"case_hwe={len(case_iids_hwe):,}, ctrl_hwe={len(ctrl_iids_hwe):,}"
    )
    case_iids_vmiss.to_csv(case_iid_path_vmiss, sep="\t", index=False, header=False)
    ctrl_iids_vmiss.to_csv(ctrl_iid_path_vmiss, sep="\t", index=False, header=False)
    case_iids_hwe.to_csv(case_iid_path_hwe, sep="\t", index=False, header=False)
    ctrl_iids_hwe.to_csv(ctrl_iid_path_hwe, sep="\t", index=False, header=False)

    # Build depth-based IID groups (30X/15X) from metadata.
    norm_dp = info_df[target_dp_col].astype(str).str.strip().str.lower()
    info_df = info_df.assign(_norm_dp=norm_dp)

    fam_cols = ["FID", "IID", "PHENO"]
    fam_sub = fam_df[fam_cols]
    merged = fam_sub.merge(info_df[[sample_col, "_norm_dp"]], left_on="IID", right_on=sample_col, how="inner")

    dp30 = merged[merged["_norm_dp"].isin(["30x", "30"])][["FID", "IID"]]
    dp15 = merged[merged["_norm_dp"].isin(["15x", "15"])][["FID", "IID"]]

    dp30_iid_path = os.path.join(tmpdir, "dp30x_iids.txt")
    dp15_iid_path = os.path.join(tmpdir, "dp15x_iids.txt")
    dp30.to_csv(dp30_iid_path, sep="\t", index=False, header=False)
    dp15.to_csv(dp15_iid_path, sep="\t", index=False, header=False)

    out_dp30 = os.path.join(tmpdir, "dp30x")
    out_dp15 = os.path.join(tmpdir, "dp15x")
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", dp30_iid_path,
        "--threads", str(threads), "--missing", "--out", out_dp30
    ], "30X subgroup VMISS")
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", dp15_iid_path,
        "--threads", str(threads), "--missing", "--out", out_dp15
    ], "15X subgroup VMISS")
    _log_info(
        "PHASE-1",
        f"Depth group sizes from metadata/FAM merge: dp30={len(dp30):,}, dp15={len(dp15):,}"
    )

    out_case_vmiss = os.path.join(tmpdir, "case_vmiss")
    out_ctrl_vmiss = os.path.join(tmpdir, "ctrl_vmiss")
    out_case_hwe = os.path.join(tmpdir, "case_hwe")
    out_ctrl_hwe = os.path.join(tmpdir, "ctrl_hwe")
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", case_iid_path_vmiss,
        "--threads", str(threads), "--freq", "--missing", "--out", out_case_vmiss
    ], "Case AAF + VMISS (full case set)")

    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", ctrl_iid_path_vmiss,
        "--threads", str(threads), "--freq", "--missing", "--out", out_ctrl_vmiss
    ], "Control AAF + VMISS (full control set)")

    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", case_iid_path_hwe,
        "--threads", str(threads), "--hardy", "--out", out_case_hwe
    ], "Case HWE (HWE exclusion applied if configured)")

    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", ctrl_iid_path_hwe,
        "--threads", str(threads), "--hardy", "--out", out_ctrl_hwe
    ], "Control HWE (HWE exclusion applied if configured)")

    # ------------------------------------------------------------------
    # Memory-efficient streaming summary construction
    # ------------------------------------------------------------------
    # Instead of loading >40M variants x multiple dictionaries into
    # memory, we stream all PLINK outputs in lockstep, chunk by chunk,
    # assuming PLINK preserves BIM order across all metrics.

    # Clamp summary parameters to sane minimums.
    if summary_workers < 1:
        summary_workers = 1
    if summary_chunksize < 10000:
        summary_chunksize = 10000

    output_file = output_prefix + ".variant_qc_summary.tsv"

    # Pre-compute VMISS thresholds if configuration is provided. When available,
    # these are pushed down into the summary chunks to materialize a PASS_VMISS
    # boolean column directly in the summary table.
    vmiss_mode_eff: Optional[str] = None
    thr_vmiss = None
    thr_case = None
    thr_ctrl = None
    thr_dp30 = None
    thr_dp15 = None
    if vmiss_mode is not None and vmiss_json_path is not None:
        try:
            import json  # local import to avoid polluting module namespace unnecessarily
            with open(vmiss_json_path, "r") as f:
                vmiss_cfg = json.load(f)
            if vmiss_mode in {"dp", "case_ctrl", "mix"} and vmiss_mode in vmiss_cfg:
                cfg = vmiss_cfg[vmiss_mode]
                vmiss_mode_eff = vmiss_mode
                if vmiss_mode == "mix":
                    thr_vmiss = float(cfg.get("VMISS", 0.05))
                elif vmiss_mode == "dp":
                    thr_dp30 = float(cfg.get("30X_VMISS", 0.05))
                    thr_dp15 = float(cfg.get("15X_VMISS", 0.05))
                elif vmiss_mode == "case_ctrl":
                    thr_ctrl = float(cfg.get("CTRL_VMISS", 0.05))
                    thr_case = float(cfg.get("CASE_VMISS", 0.05))
        except Exception as e:
            _log_warn("PHASE-1", f"Failed to parse VMISS config for PASS_VMISS annotation; skipping PASS_VMISS column: {e}")

    if verbose:
        _log_info(
            "PHASE-1",
            "Starting streaming summary build from PLINK outputs "
            f"(chunksize={summary_chunksize})"
        )

    # Prepare chunked readers for all required PLINK outputs.
    reader_all_af = pd.read_csv(out_all + ".afreq", sep=r"\s+", chunksize=summary_chunksize)
    reader_all_vmiss = pd.read_csv(out_all + ".vmiss", sep=r"\s+", usecols=["ID", "F_MISS"], chunksize=summary_chunksize)

    reader_case_af = pd.read_csv(out_case_vmiss + ".afreq", sep=r"\s+", usecols=["ID", "ALT_FREQS"], chunksize=summary_chunksize)
    reader_ctrl_af = pd.read_csv(out_ctrl_vmiss + ".afreq", sep=r"\s+", usecols=["ID", "ALT_FREQS"], chunksize=summary_chunksize)

    reader_case_vmiss = pd.read_csv(out_case_vmiss + ".vmiss", sep=r"\s+", usecols=["ID", "F_MISS"], chunksize=summary_chunksize)
    reader_ctrl_vmiss = pd.read_csv(out_ctrl_vmiss + ".vmiss", sep=r"\s+", usecols=["ID", "F_MISS"], chunksize=summary_chunksize)

    reader_case_hwe = pd.read_csv(out_case_hwe + ".hardy", sep=r"\s+", usecols=["ID", "P"], chunksize=summary_chunksize)
    reader_ctrl_hwe = pd.read_csv(out_ctrl_hwe + ".hardy", sep=r"\s+", usecols=["ID", "P"], chunksize=summary_chunksize)

    reader_dp30_vmiss = None
    reader_dp15_vmiss = None
    has_dp30_vmiss = os.path.exists(out_dp30 + ".vmiss")
    has_dp15_vmiss = os.path.exists(out_dp15 + ".vmiss")
    if has_dp30_vmiss:
        reader_dp30_vmiss = pd.read_csv(out_dp30 + ".vmiss", sep=r"\s+", usecols=["ID", "F_MISS"], chunksize=summary_chunksize)
    if has_dp15_vmiss:
        reader_dp15_vmiss = pd.read_csv(out_dp15 + ".vmiss", sep=r"\s+", usecols=["ID", "F_MISS"], chunksize=summary_chunksize)

    total_variants = 0
    chunk_files = []

    # For explicit order sanity logging: track first/last IDs
    first_chunk_first_id = None
    first_chunk_last_id = None
    last_chunk_first_id = None
    last_chunk_last_id = None
    last_chunk_index = None

    # Dispatch chunk processing to a process pool, keeping memory bounded
    # by `summary_chunksize` and the number of workers.
    with concurrent.futures.ProcessPoolExecutor(max_workers=summary_workers) as executor:
        futures = []

        for chunk_idx, all_af_chunk in enumerate(reader_all_af):
            try:
                all_vmiss_chunk = next(reader_all_vmiss)
                case_af_chunk = next(reader_case_af)
                ctrl_af_chunk = next(reader_ctrl_af)
                case_vmiss_chunk = next(reader_case_vmiss)
                ctrl_vmiss_chunk = next(reader_ctrl_vmiss)
                case_hwe_chunk = next(reader_case_hwe)
                ctrl_hwe_chunk = next(reader_ctrl_hwe)

                dp30_vmiss_chunk = next(reader_dp30_vmiss) if reader_dp30_vmiss is not None else None
                dp15_vmiss_chunk = next(reader_dp15_vmiss) if reader_dp15_vmiss is not None else None
            except StopIteration:
                raise RuntimeError("Encountered premature end of one of the PLINK metric files while streaming")

            # Sanity check: IDs should be aligned across all chunks.
            base_ids = all_af_chunk["ID"].astype(str)

            # Record boundary IDs for the first and last chunks for
            # downstream sanity checks on ordering.
            if first_chunk_first_id is None:
                if len(base_ids) > 0:
                    first_chunk_first_id = base_ids.iloc[0]
                    first_chunk_last_id = base_ids.iloc[-1]
            if len(base_ids) > 0:
                last_chunk_first_id = base_ids.iloc[0]
                last_chunk_last_id = base_ids.iloc[-1]
                last_chunk_index = chunk_idx

            def _check_ids(df, label: str) -> None:
                if not base_ids.equals(df["ID"].astype(str)):
                    raise RuntimeError(f"ID mismatch between all_samples.afreq and {label} during streaming summary build")

            _check_ids(all_vmiss_chunk, "all_samples.vmiss")
            _check_ids(case_af_chunk, "case_vmiss.afreq")
            _check_ids(ctrl_af_chunk, "ctrl_vmiss.afreq")
            _check_ids(case_vmiss_chunk, "case_vmiss.vmiss")
            _check_ids(ctrl_vmiss_chunk, "ctrl_vmiss.vmiss")
            _check_ids(case_hwe_chunk, "case_hwe.hardy")
            _check_ids(ctrl_hwe_chunk, "ctrl_hwe.hardy")
            if dp30_vmiss_chunk is not None:
                _check_ids(dp30_vmiss_chunk, "dp30x.vmiss")
            if dp15_vmiss_chunk is not None:
                _check_ids(dp15_vmiss_chunk, "dp15x.vmiss")

            future = executor.submit(
                _build_summary_chunk_worker,
                chunk_idx,
                tmpdir,
                all_af_chunk,
                all_vmiss_chunk,
                case_af_chunk,
                ctrl_af_chunk,
                case_vmiss_chunk,
                ctrl_vmiss_chunk,
                case_hwe_chunk,
                ctrl_hwe_chunk,
                dp30_vmiss_chunk,
                dp15_vmiss_chunk,
                vmiss_mode_eff,
                thr_vmiss,
                thr_case,
                thr_ctrl,
                thr_dp30,
                thr_dp15,
            )
            futures.append(future)

        total_chunks = len(futures)
        if verbose:
            _log_info(
                "PHASE-1",
                f"Dispatching {total_chunks} frequency chunks to workers (chunksize={summary_chunksize}, workers={summary_workers})"
            )

        # Start timing the streaming summary build for progress and ETA logs.
        stream_start = time.time()

        for i, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            chunk_path, n_rows, idx = future.result()
            chunk_files.append((idx, chunk_path, n_rows))
            total_variants += n_rows

            if verbose and (i % 10 == 0 or i == total_chunks):
                elapsed = max(time.time() - stream_start, 1e-6)
                # Average throughput based on processed chunks so far.
                avg_variants_per_sec = total_variants / elapsed
                est_total_variants = (total_variants / i) * total_chunks
                remaining_variants = max(est_total_variants - total_variants, 0)
                eta_sec = remaining_variants / avg_variants_per_sec if avg_variants_per_sec > 0 else None

                if eta_sec is not None and np.isfinite(eta_sec):
                    eta_min = eta_sec / 60.0
                    eta_str = f"ETA≈{eta_min:.1f} min"
                else:
                    eta_str = "ETA≈unknown"

                _log_info(
                    "PHASE-1",
                    f"Processed {i}/{total_chunks} chunks (~{total_variants:,} variants, "
                    f"speed≈{avg_variants_per_sec:,.0f} variants/s, {eta_str}) for summary table"
                )

    # Optional explicit boundary logging for additional order sanity.
    if verbose and total_chunks > 0 and first_chunk_first_id is not None and last_chunk_index is not None:
        _log_info(
            "PHASE-1",
            "Streaming summary ID boundaries: "
            f"first_chunk[0]: first_id={first_chunk_first_id}, last_id={first_chunk_last_id}; "
            f"last_chunk[{last_chunk_index}]: first_id={last_chunk_first_id}, last_id={last_chunk_last_id}"
        )

    # Merge all temporary chunk files into the final summary in order.
    if verbose:
        _log_info("PHASE-1", f"Merging {len(chunk_files)} chunk files into {output_file}")

    with open(output_file, "w", newline="") as fout:
        header_cols = [
            "VARIANT_ID",
            "MAF",
            "VMISS",
            "CASE_VMISS",
            "CTRL_VMISS",
            "30X_VMISS",
            "15X_VMISS",
            "CASE_AAF",
            "CASE_MAF",
            "CTRL_AAF",
            "CTRL_MAF",
            "CASE_HWE",
            "CTRL_HWE",
        ]
        # Append PASS_VMISS when thresholds were successfully configured.
        if vmiss_mode_eff is not None:
            header_cols.append("PASS_VMISS")
        fout.write("\t".join(header_cols) + "\n")

        for idx, chunk_path, _n_rows in sorted(chunk_files, key=lambda x: x[0]):
            with open(chunk_path, "r") as fin:
                for line in fin:
                    fout.write(line)

    if verbose:
        _log_info("PHASE-1", f"QC summary written: {output_file} (total_variants={total_variants:,})")

    return output_file



def plot_vmiss_scatter(
    variant_qc_summary: str,
    vmiss_json_path: str,
    mode: str = "dp",  # "dp" or "case_ctrl"
    variant_id_col: str = "VARIANT_ID",
    output_tsv: str = "vmiss_pass_variants.tsv",
    output_prefix: str = "cteph_agp3k",
    plot_style: str = "hex",          # "hex" | "hist2d" | "scatter" | "kde2d" | "density"
    gridsize: int = 75,
    hist_bins: int = 75,
    density_norm: str = "log",
    bw_adjust: float = 1.0,
    density_thresh: float = 0.02,
    overlay_points: bool = False
) -> str:
    """Plot VMISS QC relationships and export pass variants.

    Modes:
    - `dp`: 30X_VMISS vs 15X_VMISS
    - `case_ctrl`: CTRL_VMISS vs CASE_VMISS
    """
    import json
    import gc
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import seaborn as sns
    from matplotlib.colors import LogNorm, LinearSegmentedColormap, to_rgba
    import pandas as pd
    import numpy as np

    if mode not in {"dp", "case_ctrl"}:
        raise ValueError("mode must be 'dp' or 'case_ctrl'")

    with open(vmiss_json_path, "r") as f:
        cfg = json.load(f)
    if mode not in cfg:
        raise KeyError(f"No configuration found for mode '{mode}'")

    if mode == "dp":
        x_col, y_col = "30X_VMISS", "15X_VMISS"
        thr_x = float(cfg[mode].get("30X_VMISS", 0.05))
        thr_y = float(cfg[mode].get("15X_VMISS", 0.05))
    else:
        x_col, y_col = "CTRL_VMISS", "CASE_VMISS"
        thr_x = float(cfg[mode].get("CTRL_VMISS", 0.05))
        thr_y = float(cfg[mode].get("CASE_VMISS", 0.05))

    usecols = [variant_id_col, x_col, y_col]
    
    data_chunks = []
    reader = pd.read_csv(variant_qc_summary, sep="	", usecols=usecols, chunksize=500000)
    for chunk in reader:
        data_chunks.append(chunk)
        del chunk
        gc.collect()
    
    df = pd.concat(data_chunks, ignore_index=True)

    def _make_density_cmap(base_hex: str) -> LinearSegmentedColormap:
        c0 = to_rgba(base_hex, 0.15)
        c1 = to_rgba(base_hex, 0.40)
        c2 = to_rgba(base_hex, 0.80)
        c3 = to_rgba(base_hex, 1.00)
        return LinearSegmentedColormap.from_list("density_" + base_hex, [c0, c1, c2, c3])

    fig = plt.figure(figsize=(8, 8))
    outer_gs = gridspec.GridSpec(1, 1, width_ratios=[1], height_ratios=[1])
    
    x_raw = pd.to_numeric(df[x_col], errors="coerce")
    y_raw = pd.to_numeric(df[y_col], errors="coerce")
    pass_mask = (x_raw <= thr_x) & (y_raw <= thr_y)
    
    q_pass = pass_mask.sum()
    q_fail_x = ((x_raw > thr_x) & (y_raw <= thr_y)).sum()
    q_fail_y = ((x_raw <= thr_x) & (y_raw > thr_y)).sum()
    q_fail_both = ((x_raw > thr_x) & (y_raw > thr_y)).sum()
    
    x = x_raw
    y = y_raw
    valid = x.notna() & y.notna()
    xv = x[valid]
    yv = y[valid]
    
    color = "#1f77b4"
    cmap_local = _make_density_cmap(color)
    norm_obj = LogNorm(vmin=1) if density_norm == "log" else None

    inner_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer_gs[0, 0],
                                                width_ratios=[6, 1.5], height_ratios=[1.5, 6],
                                                wspace=0.05, hspace=0.05)
    ax_main = plt.Subplot(fig, inner_gs[1, 0])
    ax_top = plt.Subplot(fig, inner_gs[0, 0], sharex=ax_main)
    ax_right = plt.Subplot(fig, inner_gs[1, 1], sharey=ax_main)
    fig.add_subplot(ax_main)
    fig.add_subplot(ax_top)
    fig.add_subplot(ax_right)

    if plot_style == "hex":
        hb = ax_main.hexbin(xv, yv, extent=[0, 1, 0, 1], gridsize=gridsize, mincnt=1, linewidths=0, norm=norm_obj, cmap=cmap_local)
    elif plot_style == "hist2d":
        h = ax_main.hist2d(xv, yv, bins=hist_bins, range=[[0, 1], [0, 1]], norm=norm_obj, cmap=cmap_local)
    elif plot_style == "density":
        sns.kdeplot(x=xv, y=yv, ax=ax_main, fill=True, thresh=density_thresh, levels=100, bw_adjust=bw_adjust, cmap=cmap_local, linewidths=0)
        if overlay_points:
            ax_main.scatter(xv, yv, alpha=0.05, c=color, s=2, linewidths=0)
    elif plot_style == "kde2d":
        sns.kdeplot(x=xv, y=yv, ax=ax_main, fill=True, thresh=0, levels=30, bw_adjust=bw_adjust, cmap=cmap_local, linewidths=0.8)
    else:
        ax_main.scatter(xv, yv, alpha=0.1, c=color, s=5)

    ax_main.axvline(x=thr_x, color='red', linestyle='--')
    ax_main.axhline(y=thr_y, color='blue', linestyle='--')
    ax_main.set_xlim(0, 1)
    ax_main.set_ylim(0, 1)
    ax_main.set_xlabel(f"{x_col}")
    ax_main.set_ylabel(f"{y_col}")
    ax_main.grid(True, linestyle=":", linewidth=0.5)

    x_plot = xv
    if x_plot.nunique() >= 2:
        sns.kdeplot(x=x_plot, ax=ax_top, fill=True, color='gray', linewidth=1.5, cut=0)
        ax_top.set_xlim(0, 1)
        ax_top.axvline(x=thr_x, color='red', linestyle='--')
    else:
        ax_top.text(0.5, 0.5, 'KDE skipped', ha='center', va='center', transform=ax_top.transAxes, fontsize=8)
    ax_top.set_xlabel('')
    ax_top.set_ylabel('')
    ax_top.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
    ax_top.grid(False)

    y_plot = yv
    if y_plot.nunique() >= 2:
        sns.kdeplot(y=y_plot, ax=ax_right, fill=True, color='gray', linewidth=1.5, cut=0)
        ax_right.set_ylim(0, 1)
        ax_right.axhline(y=thr_y, color='blue', linestyle='--')
    else:
        ax_right.text(0.5, 0.5, 'KDE skipped', ha='center', va='center', transform=ax_right.transAxes, fontsize=8)
    ax_right.set_xlabel('')
    ax_right.set_ylabel('')
    ax_right.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
    ax_right.grid(False)

    top_ylim = ax_top.get_ylim()
    right_xlim = ax_right.get_xlim()
    inset_ax = inset_axes(ax_main, width="25%", height="25%", loc='upper right',
                          bbox_to_anchor=(0.28, 0.28, 1, 1), bbox_transform=ax_main.transAxes, borderpad=0)
    inset_ax.set_xlim(right_xlim)
    inset_ax.set_ylim(top_ylim)
    for spine in inset_ax.spines.values():
        spine.set_visible(True)
    inset_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    inset_ax.grid(False)
    xmid = (right_xlim[0] + right_xlim[1]) / 2
    ymid = (top_ylim[0] + top_ylim[1]) / 2
    x_left = (right_xlim[0] + xmid) / 2
    x_right = (xmid + right_xlim[1]) / 2
    y_bottom = (top_ylim[0] + ymid) / 2
    y_top = (ymid + top_ylim[1]) / 2
    
    inset_ax.axvline(x=xmid, linestyle='--', color='red', linewidth=1.5)
    inset_ax.axhline(y=ymid, linestyle='--', color='blue', linewidth=1.5)
    inset_ax.text(x_left, y_top, f'{q_fail_y:,}', ha='center', va='center', fontsize=7, fontweight='bold')
    inset_ax.text(x_right, y_top, f'{q_fail_both:,}', ha='center', va='center', fontsize=7, fontweight='bold')
    inset_ax.text(x_left, y_bottom, f'{q_pass:,}', ha='center', va='center', fontsize=7, fontweight='bold')
    inset_ax.text(x_right, y_bottom, f'{q_fail_x:,}', ha='center', va='center', fontsize=7, fontweight='bold')

    title_suffix = "30X_VMISS vs 15X_VMISS" if mode == "dp" else "CTRL_VMISS vs CASE_VMISS"
    fig.suptitle(f"{title_suffix}\n{x_col} Threshold ≤ {thr_x} | {y_col} Threshold ≤ {thr_y}", fontsize=16)

    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
    out_png = f"{output_prefix}.vmiss.{mode}.png"
    plt.savefig(out_png, dpi=600)
    # plt.show()

    if variant_id_col not in df.columns:
        _log_warn("PHASE-2", f"Column not found, skipping pass export: {variant_id_col}")
        return output_tsv

    pass_variants = df.loc[pass_mask, variant_id_col].dropna().unique()
    pd.Series(pass_variants).to_csv(output_tsv, sep="\t", index=False, header=False)
    return output_tsv