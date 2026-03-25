#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Professional CLI Entry Point for VMISS Threshold Tuning Pipeline
"""

import os
import sys
import argparse
import logging
import json
import time
import multiprocessing as mp
from pathlib import Path

# Ensure the local scripts directory is in the import path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from tuning_vmiss_tools import calculate_variant_metrics, analyze_vmiss_thresholds

logger = logging.getLogger(__name__)


def setup_logging(log_file: str):
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.handlers.clear()

    formatter = logging.Formatter(
        fmt='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    root.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setFormatter(formatter)
    root.addHandler(file_handler)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Variant Missingness (VMISS) Calculation and Threshold Tuning Pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 1. Main Pipeline Configuration
    pipeline_group = parser.add_argument_group('Pipeline Configuration')
    pipeline_group.add_argument("--out-dir", "-o", default="./vmiss_results", help="Output directory for all results")
    pipeline_group.add_argument("--threads", "-t", type=int, default=16, help="Number of parallel threads/processes")
    pipeline_group.add_argument("--skip-calc", action="store_true", help="Skip computation step and run analysis on existing results")
    pipeline_group.add_argument("--skip-analysis", action="store_true", help="Only run computation, skip threshold analysis")
    pipeline_group.add_argument("--log-file", help="Path to pipeline log file (default: <out-dir>/run_vmiss_tuning.log)")
    
    # 2. Input Data Configuration (Required for calculation step)
    data_group = parser.add_argument_group('Input Data Configuration')
    data_group.add_argument("--bfile", "-b", help="Prefix of the input PLINK binary files (no .bed/.bim/.fam extension)")
    data_group.add_argument("--info", "-i", help="Path to the sample phenotype/grouping metadata file (.xlsx/.csv/.txt)")
    
    # 3. Calculation Specific Parameters
    calc_group = parser.add_argument_group('Calculation Parameters')
    calc_group.add_argument("--dp-col", required=False, help="Column name in the info file indicating sequencing depth (e.g., 'Target DP (JHRPv4)')")
    calc_group.add_argument("--pheno-col", required=False, help="Column name in the info file indicating phenotypes or groups (e.g., 'OUTCOME2')")
    calc_group.add_argument("--case-value", required=False, help="Value in phenotype column representing cases (e.g., 'CTEPH')")
    calc_group.add_argument("--ctrl-value", required=False, help="Value in phenotype column representing controls (e.g., 'AGP3K')")
    calc_group.add_argument("--maf-group", "-m", default="ctrl", choices=["ctrl", "case", "all"], help="Subgroup used to calculate Minor Allele Frequency (MAF)")
    calc_group.add_argument("--metrics-file", help="Path to existing metrics file (used only if --skip-calc is set)")
    calc_group.add_argument("--id-col", default="ID", help="Sample ID column name in metadata file")
    calc_group.add_argument("--plink2-path", default="plink2", help="plink2 executable path or command")

    # 4. Analysis Specific Parameters
    analysis_group = parser.add_argument_group('Analysis Parameters')
    analysis_group.add_argument("--analysis-mode", type=int, choices=[1, 2, 3, 4], default=4, help="Analysis mode (1: All-VMISS, 2: All-15X/30X, 3: MAF-VMISS, 4: MAF-15X/30X)")
    analysis_group.add_argument("--maf-thresholds", type=float, nargs=2, default=[0.01, 0.05], metavar=('RARE', 'COMMON'), help="Thresholds for grouping by MAF")
    analysis_group.add_argument("--cdf-step", type=float, default=None, help="CDF calculation step size for modes 1 and 3 (default: 0.001)")
    analysis_group.add_argument("--cdf-step-15x", type=float, default=None, help="CDF calculation step size for 15X group in modes 2 and 4 (default: 0.001)")
    analysis_group.add_argument("--cdf-step-30x", type=float, default=None, help="CDF calculation step size for 30X group in modes 2 and 4 (default: 0.001)")
    analysis_group.add_argument("--knee-weight-x", type=float, default=1.0, help="Kneedle algorithm X-axis weight (default: 1.0)")
    analysis_group.add_argument("--knee-weight-y", type=float, default=1.0, help="Kneedle algorithm Y-axis weight (default: 1.0)")
    analysis_group.add_argument("--threshold-decimals", type=int, default=2, help="Number of decimal places for threshold rounding (default: 2)")
    
    args = parser.parse_args()
    return args

def print_banner(args):
    print("\n" + "="*70)
    print("   GENOTYPE VARIANT MISSINGNESS (VMISS) THRESHOLD TUNING PIPELINE   ")
    print("   Publication-Quality Analysis & Visualization                     ")
    print("="*70)
    print(f" Output Directory : {os.path.abspath(args.out_dir)}")
    print(f" Processing Threads: {args.threads}")
    print(f" Pipeline Steps   : Calculate={'SKIP' if args.skip_calc else 'RUN'}, Analyze={'SKIP' if args.skip_analysis else 'RUN'}")
    if not args.skip_calc:
        print(f" PLINK Bfile      : {args.bfile}")
        print(f" Metadata Info    : {args.info}")
        print(f" DP Column        : {args.dp_col}")
        print(f" Pheno Column     : {args.pheno_col}")
        print(f" Case Value       : {args.case_value}")
        print(f" Ctrl Value       : {args.ctrl_value}")
        print(f" MAF Group        : {args.maf_group}")
        print(f" ID Column        : {args.id_col}")
        print(f" plink2 Path      : {args.plink2_path}")
    print("="*70 + "\n")


def validate_args(args):
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")

    max_threads = max(1, mp.cpu_count())
    if args.threads > max_threads:
        logger.warning(f"--threads={args.threads} exceeds available CPUs ({max_threads}), clamping to {max_threads}")
        args.threads = max_threads

    maf_low, maf_mid = args.maf_thresholds
    if not (0 <= maf_low <= 1 and 0 <= maf_mid <= 1):
        raise ValueError("--maf-thresholds must be in [0, 1]")
    if maf_low >= maf_mid:
        raise ValueError("--maf-thresholds requires RARE < COMMON")

    if args.skip_calc and args.skip_analysis:
        raise ValueError("Nothing to do: both --skip-calc and --skip-analysis are enabled")

    if not args.skip_calc:
        if not args.bfile or not args.info:
            raise ValueError("--bfile and --info are required unless --skip-calc is specified")
        for ext in ('.bed', '.bim', '.fam'):
            required = Path(f"{args.bfile}{ext}")
            if not required.exists():
                raise FileNotFoundError(f"Missing PLINK input file: {required}")

        info_file = Path(args.info)
        if not info_file.exists():
            raise FileNotFoundError(f"Metadata file not found: {info_file}")

        if not args.dp_col:
            raise ValueError("--dp-col is required for calculation step")
        if not args.pheno_col:
            raise ValueError("--pheno-col is required for calculation step")
        if not args.case_value:
            raise ValueError("--case-value is required for calculation step")
        if not args.ctrl_value:
            raise ValueError("--ctrl-value is required for calculation step")

    if args.skip_calc and not args.skip_analysis and not args.metrics_file:
        default_metrics = Path(args.out_dir) / "variant_metrics_compiled.txt"
        if default_metrics.exists():
            args.metrics_file = str(default_metrics)
        else:
            raise FileNotFoundError(
                "When using --skip-calc, provide --metrics-file or ensure "
                "variant_metrics_compiled.txt exists in --out-dir"
            )

    if args.metrics_file and not Path(args.metrics_file).exists():
        raise FileNotFoundError(f"metrics file not found: {args.metrics_file}")

def main():
    args = parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = Path(args.log_file).resolve() if args.log_file else out_dir / "run_vmiss_tuning.log"
    setup_logging(str(log_path))

    started_at = time.time()
    
    try:
        validate_args(args)
        # Re-assign args.metrics_file to metrics_output so that it correctly passes to analyze_vmiss_thresholds
        metrics_output = args.metrics_file
        print_banner(args)

        if not args.skip_calc:
            logger.info(">>> [STEP 1/2] Initiating Variant Metric Calculations (PLINK2 backend) <<<")
            metrics_output = calculate_variant_metrics(
                bed_prefix=args.bfile,
                info_path=args.info,
                dp_column=args.dp_col,
                pheno_column=args.pheno_col,
                case_value=args.case_value,
                ctrl_value=args.ctrl_value,
                output_dir=str(out_dir),
                id_column=args.id_col,
                maf_group=args.maf_group,
                plink2_path=args.plink2_path,
                n_threads=args.threads
            )
            logger.info(f"✓ Calculation step completed. Metrics saved to: {metrics_output}")
        else:
            logger.info(f">>> [STEP 1/2] SKIPPED. Using existing metrics: {metrics_output} <<<")

        analysis_results = None
        if not args.skip_analysis:
            logger.info("\n>>> [STEP 2/2] Initiating VMISS Cumulative Analysis & Threshold Tuning <<<")
            analysis_results = analyze_vmiss_thresholds(
                variant_metrics_file=metrics_output,
                output_dir=str(out_dir),
                analysis_mode=args.analysis_mode,
                maf_thresholds=tuple(args.maf_thresholds),
                cdf_step=args.cdf_step,
                cdf_step_15x=args.cdf_step_15x,
                cdf_step_30x=args.cdf_step_30x,
                knee_weight_x=args.knee_weight_x,
                knee_weight_y=args.knee_weight_y,
                threshold_decimals=args.threshold_decimals,
                dpi=600,  # Publication-quality resolution (high-quality for journals)
                n_threads=args.threads
            )
            logger.info("✓ Threshold analysis and publication-quality plots completed successfully.")
        else:
            logger.info("\n>>> [STEP 2/2] SKIPPED Threshold Analysis. <<<")

        elapsed = round(time.time() - started_at, 3)
        manifest = {
            'status': 'success',
            'elapsed_seconds': elapsed,
            'out_dir': str(out_dir),
            'log_file': str(log_path),
            'metrics_file': metrics_output,
            'analysis_json': analysis_results.get('json_path') if isinstance(analysis_results, dict) else None,
            'analysis_plot': analysis_results.get('plot_path') if isinstance(analysis_results, dict) else None,
            'args': vars(args)
        }
        manifest_path = out_dir / 'run_manifest.json'
        with open(manifest_path, 'w', encoding='utf-8') as manifest_f:
            json.dump(manifest, manifest_f, indent=2, ensure_ascii=False)

        print("\n" + "="*70)
        print(" ✓ PIPELINE EXECUTION FINISHED SUCCESSFULLY ")
        print(" All outputs meet publication-quality standards")
        print("="*70 + "\n")
        logger.info(f"Run manifest written: {manifest_path}")

    except KeyboardInterrupt:
        logger.error("Interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        elapsed = round(time.time() - started_at, 3)
        fail_manifest = {
            'status': 'failed',
            'elapsed_seconds': elapsed,
            'out_dir': str(out_dir),
            'log_file': str(log_path),
            'error': str(e),
            'args': vars(args)
        }
        manifest_path = out_dir / 'run_manifest.json'
        with open(manifest_path, 'w', encoding='utf-8') as manifest_f:
            json.dump(fail_manifest, manifest_f, indent=2, ensure_ascii=False)
        sys.exit(1)

if __name__ == "__main__":
    main()
