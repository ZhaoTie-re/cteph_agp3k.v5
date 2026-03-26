#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Variant Quality Control Pipeline (CLI Version)
==============================================
A robust, professional, and concise command-line tool for large-scale genotype variant QC.

Author: ZHAO TIE
"""

import os
import sys
import json
import time
import logging
import argparse
import shutil
from pathlib import Path


def setup_logger(log_file=None):
    """Setup a formatted, clear logger to stdout and optional file.

    All VariantQC submodules (Metrics/Filters) log via the root logger,
    so we configure handlers on the root and return a namespaced logger
    for convenience.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Avoid duplicate handlers when running in embedded environments.
    if logger.handlers:
        logger.handlers.clear()

    fmt = logging.Formatter(
        "[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    stream_handler = logging.StreamHandler(stream=sys.stdout)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    if log_file is not None:
        try:
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(fmt)
            logger.addHandler(file_handler)
        except OSError:
            # Fallback gracefully if file cannot be created; keep console logs.
            logger.warning(f"[INIT] Failed to create log file: {log_file}")

    return logging.getLogger("VariantQC")

def parse_args():
    """Parse command-line arguments.

    Naming is aligned with other pipeline scripts (run_sample_qc.py,
    build_sample_qc_metrics_table.py, run_pihat_network_qc.py), while
    preserving backward compatibility with the original underscore-style
    options used in this workflow.
    """
    parser = argparse.ArgumentParser(
        description="Comprehensive Variant QC Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Core PLINK IO
    parser.add_argument(
        "--bed-prefix", "--bed_prefix",
        dest="bed_prefix",
        required=True,
        help="Prefix of input PLINK binary files (--bfile prefix)",
    )
    parser.add_argument(
        "--out-prefix", "--output_prefix",
        dest="output_prefix",
        required=True,
        help="Output file prefix for all Variant QC products",
    )

    # Sample metadata and phenotype information
    parser.add_argument(
        "--sample-info-xlsx", "--info_path",
        dest="info_path",
        required=True,
        help="Sample info Excel/CSV path (same file as other QC steps)",
    )
    parser.add_argument(
        "--sample-id-col", "--sample_col",
        dest="sample_col",
        default="ID",
        help="Column name containing sample IDs in sample info",
    )
    parser.add_argument(
        "--target-dp-col", "--target_dp_col",
        dest="target_dp_col",
        default="DEPTH",
        help="Column name containing target depth in sample info",
    )
    parser.add_argument(
        "--phenotype-col", "--pheno_col",
        dest="pheno_col",
        default="STATUS",
        help="Column name containing phenotype labels in sample info",
    )
    parser.add_argument(
        "--case-value", "--case_value",
        dest="case_value",
        default="CASE",
        help="Phenotype value representing cases",
    )
    parser.add_argument(
        "--ctrl-value", "--ctrl_value",
        dest="ctrl_value",
        default="CTRL",
        help="Phenotype value representing controls",
    )

    # Variant-level QC configuration (VMISS / HWE / MAF)
    parser.add_argument(
        "--vmiss-config", "--vmiss_json_path",
        dest="vmiss_json_path",
        default="vqc_config_vmiss.json",
        help="VMISS thresholds JSON config (same format as variant_qc_vmiss_config)",
    )
    parser.add_argument(
        "--vmiss-mode", "--vmiss_mode",
        dest="vmiss_mode",
        default="dp",
        choices=["dp", "case_ctrl", "mix"],
        help="VMISS evaluation mode: dp, case_ctrl, or mix",
    )
    parser.add_argument(
        "--hwe-config", "--hwe_json_path",
        dest="hwe_json_path",
        default="vqc_config_hwe.json",
        help="HWE thresholds JSON config (same format as variant_qc_hwe_config)",
    )
    parser.add_argument(
        "--maf-group", "--maf_group",
        dest="maf_group",
        default="ctrl",
        choices=["ctrl", "case", "all"],
        help="MAF grouping for HWE evaluation (ctrl, case, or all)",
    )
    parser.add_argument(
        "--maf-threshold", "--maf_threshold",
        dest="maf_threshold",
        type=float,
        default=0.01,
        help="MAF threshold separating low vs high MAF bins",
    )

    # Optional final MAF stratification of passed variants
    parser.add_argument(
        "--stratify-by-maf", "--stratify_by_maf",
        dest="stratify_by_maf",
        action="store_true",
        default=True,
        help="Stratify final genotype files by MAF (default: enabled)",
    )
    parser.add_argument(
        "--no-stratify-by-maf", "--no_stratify_by_maf",
        dest="stratify_by_maf",
        action="store_false",
        help="Disable MAF stratification for final genotype files",
    )
    parser.add_argument(
        "--stratify-maf-source", "--stratify_maf_source",
        dest="stratify_maf_source",
        default=None,
        choices=["ctrl", "case", "all"],
        help="MAF source for stratification (default: same as --maf-group)",
    )
    parser.add_argument(
        "--stratify-maf-threshold", "--stratify_maf_threshold",
        dest="stratify_maf_threshold",
        type=float,
        default=0.01,
        help="MAF threshold used only for final stratification",
    )

    # Runtime / environment settings
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of CPU threads to use for PLINK2 and plotting",
    )
    parser.add_argument(
        "--summary-workers",
        dest="summary_workers",
        type=int,
        default=6,
        help="Number of worker processes for building variant_qc_summary.tsv",
    )
    parser.add_argument(
        "--summary-chunksize",
        dest="summary_chunksize",
        type=int,
        default=100000,
        help="Chunk size (number of variants) per worker task when building variant_qc_summary.tsv",
    )
    parser.add_argument(
        "--tmpdir",
        default="./variant_qc_tmp",
        help="Temporary working directory for intermediate PLINK2 outputs",
    )
    parser.add_argument(
        "--delete-tmpdir", "--delete_tmpdir",
        dest="delete_tmpdir",
        action="store_true",
        help="Delete temporary directory after successful completion",
    )
    parser.add_argument(
        "--plink2-path", "--plink2_path",
        dest="plink2_path",
        default="plink2",
        help="Path to plink2 executable",
    )
    parser.add_argument(
        "--script-path", "--script_path",
        dest="script_path",
        default=str(Path(__file__).resolve().parent),
        help="Path to Python modules (variant_qc_metrics.py, variant_qc_filters.py)",
    )
    parser.add_argument(
        "--pihat-vertex-cover-tsv-for-hwe", "--pihat_vertex_cover_tsv_for_hwe",
        dest="pihat_vertex_cover_tsv_for_hwe",
        default=None,
        help=(
            "Optional PIHAT vertex-cover TSV; when provided, samples with "
            "SELECTED_FOR_REMOVAL=true are excluded from HWE calculations only"
        ),
    )

    return parser.parse_args()

def main():
    t0 = time.time()
    args = parse_args()

    # Create a dedicated log file alongside standard stdout logging.
    log_file = f"{args.output_prefix}.log"
    logger = setup_logger(log_file=log_file)

    logger.info("[INIT] Initializing Variant QC pipeline")
    logger.info(f"[INIT] Logging to file: {os.path.abspath(log_file)}")

    logger.info("[INIT] Configuration summary (CLI-style names):")
    logger.info(f"[INIT]   --bed-prefix              = {args.bed_prefix}")
    logger.info(f"[INIT]   --out-prefix              = {args.output_prefix}")
    logger.info(f"[INIT]   --sample-info-xlsx        = {args.info_path}")
    logger.info(f"[INIT]   --sample-id-col           = {args.sample_col}")
    logger.info(f"[INIT]   --target-dp-col           = {args.target_dp_col}")
    logger.info(f"[INIT]   --phenotype-col           = {args.pheno_col} (case={args.case_value}, ctrl={args.ctrl_value})")
    logger.info(f"[INIT]   --vmiss-config            = {args.vmiss_json_path} (mode={args.vmiss_mode})")
    logger.info(f"[INIT]   --hwe-config              = {args.hwe_json_path} (maf_group={args.maf_group}, maf_threshold={args.maf_threshold})")
    logger.info(f"[INIT]   --threads                 = {args.threads}")
    logger.info(f"[INIT]   --summary-workers         = {args.summary_workers}")
    logger.info(f"[INIT]   --summary-chunksize       = {args.summary_chunksize}")
    logger.info(f"[INIT]   --tmpdir                  = {args.tmpdir}")
    if args.pihat_vertex_cover_tsv_for_hwe:
        logger.info(f"[INIT]   PIHAT HWE exclusion: ENABLED (vertex_cover_tsv={args.pihat_vertex_cover_tsv_for_hwe})")
    else:
        logger.info("[INIT]   PIHAT HWE exclusion: DISABLED (no vertex-cover TSV provided)")
    script_dir = os.path.abspath(args.script_path)
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    try:
        from variant_qc_metrics import run_plink2_variant_qc, plot_vmiss_scatter
        from variant_qc_filters import plot_vmiss_distribution, plot_hwe_scatter_by_maf_category, extract_pass_variants_by_intersection
    except ImportError as e:
        logger.error(f"[INIT] Failed to load modules from {script_dir}: {e}")
        sys.exit(1)

    logger.info(f"[PHASE-1] Running PLINK2 variant QC ({args.threads} threads)")
    try:
        variant_qc_summary = run_plink2_variant_qc(
            bed_prefix=args.bed_prefix,
            tmpdir=args.tmpdir,
            plink2_path=args.plink2_path,
            threads=args.threads,
            output_prefix=args.output_prefix,
            verbose=True,
            info_path=args.info_path,
            sample_col=args.sample_col,
            target_dp_col=args.target_dp_col,
            pheno_col=args.pheno_col,
            case_value=args.case_value,
            ctrl_value=args.ctrl_value,
            pihat_vertex_cover_tsv_for_hwe=args.pihat_vertex_cover_tsv_for_hwe,
            summary_workers=args.summary_workers,
            summary_chunksize=args.summary_chunksize,
            vmiss_mode=args.vmiss_mode,
            vmiss_json_path=args.vmiss_json_path,
        )
        logger.info(f"[PHASE-1] Completed. QC summary: {variant_qc_summary}")
    except Exception as e:
        logger.exception(f"[PHASE-1] Failed: {e}")
        sys.exit(1)

    maf_col_map = {"ctrl": "CTRL_MAF", "case": "CASE_MAF", "all": "MAF"}
    maf_col = maf_col_map.get(args.maf_group, "CTRL_MAF")

    pass_vmiss_path = None
    pass_vmiss_tsv = f"{args.output_prefix}.vmiss_pass_variants.tsv"
    logger.info(f"[PHASE-2] Evaluating VMISS (mode={args.vmiss_mode})")
    try:
        if args.vmiss_mode in ("dp", "case_ctrl"):
            pass_vmiss_path = plot_vmiss_scatter(
                variant_qc_summary=variant_qc_summary,
                vmiss_json_path=args.vmiss_json_path,
                mode=args.vmiss_mode,
                output_tsv=pass_vmiss_tsv,
                plot_style="hex",
                density_norm="log",
                output_prefix=args.output_prefix,
            )
        elif args.vmiss_mode == "mix":
            with open(args.vmiss_json_path, "r") as f:
                vmiss_threshold = float(json.load(f)["mix"]["VMISS"])
            pass_vmiss_path = plot_vmiss_distribution(
                variant_qc_summary=variant_qc_summary,
                vmiss_threshold=vmiss_threshold,
                output_tsv=pass_vmiss_tsv,
                output_prefix=args.output_prefix,
            )
        logger.info(f"[PHASE-2] Completed. VMISS pass variants: {pass_vmiss_path}")
    except Exception as e:
        logger.exception(f"[PHASE-2] Failed: {e}")
        sys.exit(1)

    logger.info("[PHASE-3] Evaluating HWE")
    pass_hwe_tsv = f"{args.output_prefix}.hwe_pass_variants.tsv"
    try:
        with open(args.hwe_json_path, "r") as f:
            hwe_thresholds = json.load(f)
        pass_hwe_path = plot_hwe_scatter_by_maf_category(
            variant_qc_summary=variant_qc_summary,
            hwe_thresholds=hwe_thresholds,
            maf_col=maf_col,
            maf_threshold=args.maf_threshold,
            output_tsv=pass_hwe_tsv,
            output_prefix=args.output_prefix,
            pass_vmiss_path=pass_vmiss_path,
        )
        logger.info(f"[PHASE-3] Completed. HWE pass variants: {pass_hwe_path}")
    except Exception as e:
        logger.exception(f"[PHASE-3] Failed: {e}")
        sys.exit(1)

    logger.info("[PHASE-4] Intersecting pass lists and extracting PLINK subset")
    pass_intersection_tsv = f"{args.output_prefix}.pass_variants.tsv"
    try:
        pass_variants_path = extract_pass_variants_by_intersection(
            pass_vmiss_path=pass_vmiss_path,
            pass_hwe_path=pass_hwe_path,
            output_path=pass_intersection_tsv,
            bed_prefix=args.bed_prefix,
            output_prefix=args.output_prefix,
            threads=args.threads,
            plink2=args.plink2_path,
        )
        logger.info(f"[PHASE-4] Completed. Final subset prefix: {args.output_prefix}")
    except Exception as e:
        logger.exception(f"[PHASE-4] Failed: {e}")
        sys.exit(1)

    # PHASE-5: Optional MAF stratification
    if args.stratify_by_maf:
        logger.info("[PHASE-5] MAF stratification enabled")
        
        # Determine MAF source: use stratify_maf_source if specified, otherwise use maf_group
        stratify_maf_source = args.stratify_maf_source if args.stratify_maf_source is not None else args.maf_group
        stratify_maf_col = maf_col_map.get(stratify_maf_source, "CTRL_MAF")
        
        try:
            from variant_qc_filters import stratify_by_maf
            
            below_threshold_list, above_threshold_list = stratify_by_maf(
                variant_qc_summary=variant_qc_summary,
                pass_variants_path=pass_variants_path,
                maf_col=stratify_maf_col,
                maf_threshold=args.stratify_maf_threshold,
                bed_prefix=args.bed_prefix,
                output_prefix=args.output_prefix,
                threads=args.threads,
                plink2=args.plink2_path,
            )
            logger.info(f"[PHASE-5] Completed. Stratified by {stratify_maf_col} at threshold {args.stratify_maf_threshold}")
            logger.info(f"[PHASE-5] Below threshold: {below_threshold_list}")
            logger.info(f"[PHASE-5] Above threshold: {above_threshold_list}")
        except Exception as e:
            logger.exception(f"[PHASE-5] Failed: {e}")
            sys.exit(1)
    else:
        logger.info("[PHASE-5] MAF stratification disabled (--no_stratify_by_maf)")


    if args.delete_tmpdir and os.path.exists(args.tmpdir):
        try:
            shutil.rmtree(args.tmpdir)
            logger.info(f"[CLEANUP] Temporary directory deleted: {args.tmpdir}")
        except Exception as e:
            logger.warning(f"[CLEANUP] Failed to delete temporary directory {args.tmpdir}: {e}")

    logger.info(f"[DONE] Pipeline completed in {(time.time() - t0):.1f} seconds")
    print(pass_variants_path)

if __name__ == "__main__":
    main()
