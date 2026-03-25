#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genotype Variant Quality Control Tools
For calculating variant missingness (VMISS) by groups and Minor Allele Frequency (MAF)
"""

import pandas as pd
import numpy as np
import subprocess
import logging
import os
import shutil
import tempfile
import multiprocessing as mp
from typing import Optional, Literal

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def _normalize_thread_count(n_threads: int) -> int:
    """Normalize thread count to a safe range based on local CPU availability."""
    if n_threads < 1:
        logger.warning(f"Invalid n_threads={n_threads}, fallback to 1")
        return 1

    available = max(1, mp.cpu_count())
    if n_threads > available:
        logger.warning(f"Requested n_threads={n_threads} exceeds available CPUs={available}, clamping to {available}")
        return available

    return n_threads


def _resolve_plink2_executable(plink2_path: str) -> str:
    """Resolve plink2 executable path (priority: explicit parameter > PATH > fallback)"""
    fallback = "/home/b/b37974/plink2"

    if os.path.exists(plink2_path):
        return plink2_path

    found = shutil.which(plink2_path)
    if found is not None:
        return found

    if os.path.exists(fallback):
        logger.warning(
            f"plink2 not found in PATH, using fallback: {fallback} "
            f"(Recommend setting export PATH=/home/b/b37974/:$PATH in shell)"
        )
        return fallback

    raise FileNotFoundError(
        f"Cannot find plink2 command or executable: {plink2_path}, and fallback path does not exist: {fallback}"
    )


def calculate_variant_metrics(
    bed_prefix: str,
    info_path: str,
    dp_column: str,
    pheno_column: str,
    case_value: str,
    ctrl_value: str,
    output_dir: str = "./vmiss_results",
    id_column: str = "ID",
    maf_group: Literal["ctrl", "case", "all"] = "ctrl",
    plink2_path: str = "plink2",
    n_threads: int = 4
) -> str:
    """
    Calculate VMISS and MAF for each variant in genotype files
    
    Parameters:
        bed_prefix: Path prefix of PLINK2 genotype files (without .bed/.bim/.fam extension)
        info_path: Path to sample information table
        dp_column: Column name for sequencing depth (in info file)
        pheno_column: Column name for phenotypes (in info file)
        case_value: Value in phenotype column representing cases (e.g., "CTEPH")
        ctrl_value: Value in phenotype column representing controls (e.g., "AGP3K")
        output_dir: Output directory, results saved as output_dir/variant_metrics_compiled.txt
        id_column: Sample ID column name, default "ID"
        maf_group: Sample group for MAF calculation, options: "ctrl"/"case"/"all", default "ctrl"
        plink2_path: plink2 command or executable path
        n_threads: Number of parallel threads
    
    Returns:
        str: Path to tab-delimited output file containing:
            - VARIANT_ID: Variant identifier
            - VMISS: Variant missingness across all samples
            - VMISS_15X: Variant missingness for 15X depth samples
            - VMISS_30X: Variant missingness for 30X depth samples
            - MAF: Minor Allele Frequency, calculated as min(ALT_FREQ, 1-ALT_FREQ)
    
    Memory Optimization Strategy:
        1. Fully streaming: Never load complete plink2 output into memory
        2. Chunked indexing: Build dictionaries in chunks using chunksize=100000
        3. Direct write: Process line-by-line and write directly to final file
        4. Explicit cleanup: Use gc.collect() to proactively release memory
        
        This design can handle millions of variants with memory footprint determined
        primarily by variant ID dictionary size rather than full dataset size.
        Optimized for memory-constrained environments.
    """
    
    import gc
    gc.collect()  # Clean up memory before starting
    n_threads = _normalize_thread_count(n_threads)
    
    logger.info("="*60)
    logger.info("Starting variant quality metrics calculation")
    logger.info(f"Genotype file prefix: {bed_prefix}")
    logger.info(f"Sample info file: {info_path}")
    logger.info(f"MAF calculation group: {maf_group}")
    logger.info(f"Number of threads: {n_threads}")
    logger.info("="*60)
    
    # Resolve plink2 executable (handles conda run and other PATH-isolated scenarios)
    plink2_path = _resolve_plink2_executable(plink2_path)
    logger.info(f"plink2 path: {plink2_path}")

    # Check if input files exist
    _check_input_files(bed_prefix, info_path, plink2_path)
    
    # Read sample information
    logger.info("Reading sample information file...")
    sample_info = _read_sample_info(info_path, id_column, dp_column, pheno_column)
    logger.info(f"Successfully loaded {len(sample_info)} samples' information")
    
    # Read FAM file to get sample list
    logger.info("Reading FAM file...")
    fam_samples = _read_fam_file(bed_prefix)
    logger.info(f"Genotype file contains {len(fam_samples)} samples")
    
    # Merge sample information with FAM file
    logger.info("Matching sample information...")
    merged_samples = _merge_sample_info(fam_samples, sample_info, id_column, dp_column, pheno_column)
    logger.info(f"Successfully matched {len(merged_samples)} samples")
    
    # Group samples by sequencing depth
    dp_groups = _group_samples_by_dp(merged_samples, dp_column)
    logger.info(f"Sequencing depth groups: {dict(dp_groups['count'])}")
    
    # Set output path (consistent with run_vmiss_tuning.py)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "variant_metrics_compiled.txt")
    
    logger.info(f"Output file path: {output_path}")
    
    # Create temporary directory
    with tempfile.TemporaryDirectory(prefix="vmiss_calc_") as tmpdir:
        logger.info(f"Created temporary directory: {tmpdir}")
        
        # Generate sample list files
        sample_lists = _create_sample_list_files(
            tmpdir, merged_samples, dp_groups, maf_group, pheno_column, case_value, ctrl_value
        )
        
        # Calculate overall VMISS
        logger.info("Calculating overall VMISS...")
        vmiss_all_file = _calculate_vmiss_plink(
            bed_prefix, plink2_path, tmpdir, "all", 
            sample_list=None, n_threads=n_threads
        )
        logger.info(f"Successfully calculated overall VMISS")
        
        # Calculate VMISS for 15X group
        if '15X' in sample_lists:
            logger.info("Calculating VMISS for 15X group...")
            vmiss_15x_file = _calculate_vmiss_plink(
                bed_prefix, plink2_path, tmpdir, "15x",
                sample_list=sample_lists['15X'], n_threads=n_threads
            )
            logger.info(f"Successfully calculated VMISS for 15X group")
        else:
            logger.warning("No 15X group samples found")
            vmiss_15x_file = None
        
        # Calculate VMISS for 30X group
        if '30X' in sample_lists:
            logger.info("Calculating VMISS for 30X group...")
            vmiss_30x_file = _calculate_vmiss_plink(
                bed_prefix, plink2_path, tmpdir, "30x",
                sample_list=sample_lists['30X'], n_threads=n_threads
            )
            logger.info(f"Successfully calculated VMISS for 30X group")
        else:
            logger.warning("No 30X group samples found")
            vmiss_30x_file = None
        
        # Calculate MAF
        logger.info(f"Calculating MAF for {maf_group} group...")
        maf_file = _calculate_maf_plink(
            bed_prefix, plink2_path, tmpdir, maf_group,
            sample_list=sample_lists.get(maf_group), n_threads=n_threads
        )
        logger.info(f"Successfully calculated MAF")
        
        # Merge results and write directly to file (streaming processing, no memory loading)
        logger.info("Merging results and writing to file (streaming mode)...")
        n_variants = _merge_and_write_results_streaming(
            vmiss_all_file, vmiss_15x_file, vmiss_30x_file, maf_file, output_path
        )
    logger.info(f"Successfully wrote {n_variants} variants to file")
    
    logger.info("="*60)
    logger.info(f"Calculation complete! Results saved to: {output_path}")
    logger.info("="*60)
    
    return output_path


def _check_input_files(bed_prefix: str, info_path: str, plink2_path: str):
    """Check if input files exist"""
    # Check PLINK files
    for ext in ['.bed', '.bim', '.fam']:
        file_path = bed_prefix + ext
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
    
    # Check info file
    if not os.path.exists(info_path):
        raise FileNotFoundError(f"Sample information file not found: {info_path}")
    
    # Check plink2 (should be resolved executable path at this point)
    if not os.path.exists(plink2_path):
        raise FileNotFoundError(f"plink2 executable not found: {plink2_path}")
    if not os.access(plink2_path, os.X_OK):
        raise PermissionError(f"plink2 is not executable: {plink2_path}")


def _read_sample_info(info_path: str, id_column: str, dp_column: str, pheno_column: str) -> pd.DataFrame:
    """Read sample information file and ensure data quality"""
    # Select reading method based on file extension
    try:
        if info_path.endswith('.xlsx'):
            df = pd.read_excel(info_path)
        elif info_path.endswith('.csv'):
            df = pd.read_csv(info_path)
        elif info_path.endswith('.tsv') or info_path.endswith('.txt'):
            df = pd.read_csv(info_path, sep='\t')
        else:
            raise ValueError(f"Unsupported file format: {info_path}")
    except Exception as e:
        raise RuntimeError(f"Failed to read sample info file {info_path}: {e}")
    
    # Check if required columns exist
    required_cols = [id_column, dp_column]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Sample info file missing required columns: {missing_cols}\nAvailable columns: {df.columns.tolist()}")
    
    # Strip whitespace from string columns to prevent mismatch errors
    if pd.api.types.is_string_dtype(df[id_column]):
        df[id_column] = df[id_column].str.strip()
    if pd.api.types.is_string_dtype(df[dp_column]):
        df[dp_column] = df[dp_column].str.strip()

    # Remove rows with missing/empty sample IDs to avoid invalid matching/duplicate noise
    valid_id_mask = df[id_column].notna()
    if pd.api.types.is_string_dtype(df[id_column]):
        valid_id_mask &= (df[id_column] != '')

    n_missing_ids = int((~valid_id_mask).sum())
    if n_missing_ids > 0:
        logger.warning(f"Found {n_missing_ids} rows with missing/empty sample IDs in info file; dropping them.")
        df = df.loc[valid_id_mask].copy()
        
    # Check for duplicate sample IDs
    if df[id_column].duplicated().any():
        n_dupes = df[id_column].duplicated().sum()
        logger.warning(f"Found {n_dupes} duplicate sample IDs in info file! Keeping first occurrence.")
        df = df.drop_duplicates(subset=[id_column], keep='first')
    
    # Check phenotype column if case/ctrl MAF calculation is needed
    if pheno_column not in df.columns:
        logger.warning(f"Phenotype column '{pheno_column}' not found, case-control grouping will not be possible")
    else:
        if pd.api.types.is_string_dtype(df[pheno_column]):
            df[pheno_column] = df[pheno_column].str.strip()
    
    return df


def _read_fam_file(bed_prefix: str) -> pd.DataFrame:
    """Read FAM file and ensure clean sample IDs"""
    fam_path = bed_prefix + '.fam'
    # FAM file format: FID IID FATHER MOTHER SEX PHENO
    # No header, use column indices
    try:
        fam_df = pd.read_csv(
            fam_path, 
            sep=r'\s+', 
            header=None,
            names=['FID', 'IID', 'FATHER', 'MOTHER', 'SEX', 'PHENO']
        )
        
        # Ensure IID is treated as pure stripped string to match info perfectly
        fam_df['IID'] = fam_df['IID'].astype(str).str.strip()
        return fam_df
    except Exception as e:
        raise RuntimeError(f"Failed to read FAM file {fam_path}: {e}")


def _merge_sample_info(fam_df: pd.DataFrame, info_df: pd.DataFrame, 
                       id_column: str, dp_column: str, pheno_column: str) -> pd.DataFrame:
    """Merge FAM file with sample information"""
    # Match FAM IID column with info file ID column
    merged = fam_df.merge(
        info_df[[id_column, dp_column] + ([pheno_column] if pheno_column in info_df.columns else [])],
        left_on='IID',
        right_on=id_column,
        how='left'
    )
    
    # Check how many samples were successfully matched
    n_matched = merged[dp_column].notna().sum()
    n_total = len(fam_df)
    if n_matched == 0:
        raise ValueError(f"No samples successfully matched between FAM ({n_total} samples) and info file! Please check ID column.")
    elif n_matched < n_total:
        n_missing = n_total - n_matched
        logger.warning(f"Failed to match {n_missing} samples out of {n_total}. Missing samples will be excluded from group statistics.")
    
    logger.info(f"Successfully matched {n_matched}/{n_total} samples' depth information")
    
    return merged


def _group_samples_by_dp(merged_df: pd.DataFrame, dp_column: str) -> dict:
    """Group samples by sequencing depth"""
    # Extract depth values (recognizes 15x/15X/30x/30X, case-insensitive)
    def extract_dp(dp_str):
        if pd.isna(dp_str):
            return None
        dp_str = str(dp_str).strip().upper()
        if '30' in dp_str or dp_str.startswith('30'):
            return '30X'
        elif '15' in dp_str or dp_str.startswith('15'):
            return '15X'
        else:
            return None
    
    merged_df['DP_GROUP'] = merged_df[dp_column].apply(extract_dp)
    
    # Count samples in each group
    group_counts = merged_df['DP_GROUP'].value_counts()
    
    return {
        'data': merged_df,
        'count': group_counts
    }


def _create_sample_list_files(tmpdir: str, merged_df: pd.DataFrame, 
                               dp_groups: dict, maf_group: str, 
                               pheno_column: str, case_value: str, ctrl_value: str) -> dict:
    """Create sample list files for each group"""
    sample_lists = {}
    
    merged_df = dp_groups['data']
    
    # Create 15X group sample list
    df_15x = merged_df[merged_df['DP_GROUP'] == '15X']
    if len(df_15x) > 0:
        list_15x = os.path.join(tmpdir, 'samples_15x.txt')
        df_15x[['FID', 'IID']].to_csv(list_15x, sep='\t', header=False, index=False)
        sample_lists['15X'] = list_15x
        logger.info(f"Created 15X group sample list: {len(df_15x)} samples")
    
    # Create 30X group sample list
    df_30x = merged_df[merged_df['DP_GROUP'] == '30X']
    if len(df_30x) > 0:
        list_30x = os.path.join(tmpdir, 'samples_30x.txt')
        df_30x[['FID', 'IID']].to_csv(list_30x, sep='\t', header=False, index=False)
        sample_lists['30X'] = list_30x
        logger.info(f"Created 30X group sample list: {len(df_30x)} samples")
    
    # Create sample lists for MAF calculation
    if maf_group == "ctrl" and pheno_column in merged_df.columns:
        # Control group: use user-specified ctrl_value
        df_ctrl = merged_df[
            merged_df[pheno_column].astype(str).str.upper() == ctrl_value.upper()
        ]
        if len(df_ctrl) > 0:
            list_ctrl = os.path.join(tmpdir, 'samples_ctrl.txt')
            df_ctrl[['FID', 'IID']].to_csv(list_ctrl, sep='\t', header=False, index=False)
            sample_lists['ctrl'] = list_ctrl
            logger.info(f"Created control group sample list ({ctrl_value}): {len(df_ctrl)} samples")
        else:
            logger.warning(f"No control group samples found ({ctrl_value}), will use all samples for MAF calculation")
    
    elif maf_group == "case" and pheno_column in merged_df.columns:
        # Case group: use user-specified case_value
        df_case = merged_df[
            merged_df[pheno_column].astype(str).str.upper() == case_value.upper()
        ]
        if len(df_case) > 0:
            list_case = os.path.join(tmpdir, 'samples_case.txt')
            df_case[['FID', 'IID']].to_csv(list_case, sep='\t', header=False, index=False)
            sample_lists['case'] = list_case
            logger.info(f"Created case group sample list ({case_value}): {len(df_case)} samples")
        else:
            logger.warning(f"No case group samples found ({case_value}), will use all samples for MAF calculation")
    
    # When maf_group is "all", no additional list is needed
    
    return sample_lists


def _calculate_vmiss_plink(bed_prefix: str, plink2_path: str, tmpdir: str,
                            suffix: str, sample_list: Optional[str] = None,
                            n_threads: int = 4) -> str:
    """
    Calculate VMISS using plink2
    Returns output file path instead of DataFrame to avoid loading large data into memory
    """
    output_prefix = os.path.join(tmpdir, f'vmiss_{suffix}')
    
    cmd = [
        plink2_path,
        '--bfile', bed_prefix,
        '--missing', 'variant-only',
        '--out', output_prefix,
        '--threads', str(n_threads)
    ]
    
    if sample_list is not None:
        cmd.extend(['--keep', sample_list])
    
    logger.info(f"Executing command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logger.debug(f"plink2 stdout:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"plink2 execution failed for vmiss_{suffix}:\nCommand: {' '.join(cmd)}\nStderr:\n{e.stderr}\nStdout:\n{e.stdout}")
        raise RuntimeError(f"plink2 failed during vmiss_{suffix} calculation.") from e
        
    # Return result file path instead of reading
    vmiss_file = output_prefix + '.vmiss'
    if not os.path.exists(vmiss_file):
        raise FileNotFoundError(f"plink2 output file not found: {vmiss_file}")
    
    # Only count lines, don't load data
    with open(vmiss_file, 'r') as f:
        n_lines = sum(1 for _ in f) - 1  # Minus header
    logger.info(f"VMISS file contains {n_lines} variants")
    
    return vmiss_file


def _calculate_maf_plink(bed_prefix: str, plink2_path: str, tmpdir: str,
                         maf_group: str, sample_list: Optional[str] = None,
                         n_threads: int = 4) -> str:
    """
    Calculate MAF using plink2
    Returns output file path instead of DataFrame to avoid loading large data into memory
    """
    output_prefix = os.path.join(tmpdir, f'freq_{maf_group}')
    
    cmd = [
        plink2_path,
        '--bfile', bed_prefix,
        '--freq',
        '--out', output_prefix,
        '--threads', str(n_threads)
    ]
    
    if sample_list is not None:
        cmd.extend(['--keep', sample_list])
    
    logger.info(f"Executing command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logger.debug(f"plink2 stdout:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"plink2 execution failed for MAF calculation:\nCommand: {' '.join(cmd)}\nStderr:\n{e.stderr}\nStdout:\n{e.stdout}")
        raise RuntimeError(f"plink2 failed during MAF calculation.") from e
    
    # Return result file path instead of reading
    freq_file = output_prefix + '.afreq'
    if not os.path.exists(freq_file):
        raise FileNotFoundError(f"plink2 output file not found: {freq_file}")
    
    # Only count lines, don't load data
    with open(freq_file, 'r') as f:
        n_lines = sum(1 for _ in f) - 1  # Minus header
    logger.info(f"Frequency file contains {n_lines} variants")
    
    return freq_file


def _merge_and_write_results_streaming(vmiss_all_file: str, vmiss_15x_file: Optional[str],
                                        vmiss_30x_file: Optional[str], maf_file: str,
                                        output_path: str) -> int:
    """
    Fully streaming processing: Read directly from plink2 output files and merge-write
    Extremely low memory footprint, suitable for processing millions of variants
    
    Parameters:
        vmiss_all_file: Overall VMISS file path
        vmiss_15x_file: 15X group VMISS file path (optional)
        vmiss_30x_file: 30X group VMISS file path (optional)
        maf_file: MAF file path
        output_path: Output file path
    
    Returns:
        int: Number of variants written
    """
    import gc
    
    logger.info("Starting streaming index dictionary construction...")
    
    # General function to build dictionary in chunks
    def build_dict_streaming(filepath, key_col, val_col, chunksize=50000):
        """Build dictionary in streaming mode, processing only a portion of data at a time"""
        result_dict = {}
        if filepath is None or not os.path.exists(filepath):
            return result_dict
        
        chunk_count = 0
        for chunk in pd.read_csv(filepath, sep=r'\s+', chunksize=chunksize, 
                                  usecols=[key_col, val_col], dtype=str):  # Only read required columns to reduce memory
            if key_col in chunk.columns and val_col in chunk.columns:
                chunk_dict = dict(zip(chunk[key_col].values, chunk[val_col].values))
                result_dict.update(chunk_dict)
                del chunk, chunk_dict
                chunk_count += 1
                if chunk_count % 5 == 0:  # More frequent garbage collection
                    gc.collect()
        
        return result_dict
    
    # Build each dictionary
    logger.info("Building 15X group index...")
    vmiss_15x_dict = build_dict_streaming(vmiss_15x_file, 'ID', 'F_MISS')
    logger.info(f"15X group index complete: {len(vmiss_15x_dict)} variants")
    
    logger.info("Building 30X group index...")
    vmiss_30x_dict = build_dict_streaming(vmiss_30x_file, 'ID', 'F_MISS')
    logger.info(f"30X group index complete: {len(vmiss_30x_dict)} variants")
    
    logger.info("Building MAF index...")
    maf_dict = build_dict_streaming(maf_file, 'ID', 'ALT_FREQS')
    logger.info(f"MAF index complete: {len(maf_dict)} variants")
    
    # Stream read main file and write results
    logger.info("Starting streaming merge and write...")
    n_variants = 0
    
    with open(output_path, 'w', buffering=8192*16) as outf:  # Increase write buffer
        # Write header
        outf.write('VARIANT_ID\tVMISS\tVMISS_15X\tVMISS_30X\tMAF\n')
        
        # Process main VMISS file in chunks, reduce chunksize to minimize memory usage
        chunksize = 20000  # Reduced from 50000 to 20000
        for chunk in pd.read_csv(vmiss_all_file, sep=r'\s+', chunksize=chunksize):
            # Use vectorized operations instead of iterrows(), faster and more memory efficient
            for variant_id, vmiss in zip(chunk['ID'].values, chunk['F_MISS'].values):
                vmiss_15x = vmiss_15x_dict.get(variant_id, 'NA')
                vmiss_30x = vmiss_30x_dict.get(variant_id, 'NA')
                
                # Get ALT frequency and convert to true MAF (Minor Allele Frequency)
                alt_freq = maf_dict.get(variant_id, 'NA')
                if alt_freq != 'NA':
                    try:
                        alt_freq_val = float(alt_freq)
                        # MAF = min(ALT_FREQ, 1 - ALT_FREQ)
                        maf = min(alt_freq_val, 1 - alt_freq_val)
                    except (ValueError, TypeError):
                        maf = 'NA'
                else:
                    maf = 'NA'
                
                outf.write(f"{variant_id}\t{vmiss}\t{vmiss_15x}\t{vmiss_30x}\t{maf}\n")
                n_variants += 1
            
            del chunk
            if n_variants % 50000 == 0:
                logger.info(f"Processed {n_variants} variants...")
                gc.collect()
    
    # Clean up dictionaries
    del vmiss_15x_dict, vmiss_30x_dict, maf_dict
    gc.collect()
    
    logger.info(f"Streaming processing complete, wrote {n_variants} variants")
    return n_variants


def _merge_and_write_results(vmiss_all: pd.DataFrame, vmiss_15x: pd.DataFrame,
                             vmiss_30x: pd.DataFrame, maf_data: pd.DataFrame,
                             output_path: str) -> int:
    """
    Fully streaming merge of results with minimal memory usage
    Does not save complete data in memory; reads and writes line by line
    
    Returns:
        int: Number of variants written
    """
    logger.info("Preparing data dictionaries (streaming processing)...")
    
    # Create dictionaries in chunks to reduce memory peaks
    # Use chunksize for chunked reading (if data already in memory, save to temp files first)
    import gc
    
    # Create temporary files to save each data
    tmpdir = os.path.dirname(output_path)
    tmp_vmiss_all = os.path.join(tmpdir, '.tmp_vmiss_all.txt')
    tmp_vmiss_15x = os.path.join(tmpdir, '.tmp_vmiss_15x.txt')
    tmp_vmiss_30x = os.path.join(tmpdir, '.tmp_vmiss_30x.txt')
    tmp_maf = os.path.join(tmpdir, '.tmp_maf.txt')
    
    # Save data to temporary files and release memory immediately
    logger.info("Saving intermediate results to temporary files...")
    vmiss_all[['ID', 'F_MISS']].to_csv(tmp_vmiss_all, sep='\t', index=False)
    del vmiss_all
    gc.collect()
    
    if len(vmiss_15x) > 0:
        vmiss_15x[['ID', 'F_MISS']].to_csv(tmp_vmiss_15x, sep='\t', index=False)
    del vmiss_15x
    gc.collect()
    
    if len(vmiss_30x) > 0:
        vmiss_30x[['ID', 'F_MISS']].to_csv(tmp_vmiss_30x, sep='\t', index=False)
    del vmiss_30x
    gc.collect()
    
    if len(maf_data) > 0:
        maf_data[['ID', 'ALT_FREQS']].to_csv(tmp_maf, sep='\t', index=False)
    del maf_data
    gc.collect()
    
    logger.info("Building index dictionaries (chunked processing)...")
    
    # Build dictionaries in chunks, processing only a portion of data at a time
    def build_dict_in_chunks(filepath, key_col, val_col, chunksize=100000):
        """Build dictionaries in chunks to avoid memory peaks"""
        result_dict = {}
        if not os.path.exists(filepath):
            return result_dict
        
        for chunk in pd.read_csv(filepath, sep='\t', chunksize=chunksize):
            chunk_dict = dict(zip(chunk[key_col], chunk[val_col]))
            result_dict.update(chunk_dict)
            del chunk, chunk_dict
            gc.collect()
        
        return result_dict
    
    vmiss_15x_dict = build_dict_in_chunks(tmp_vmiss_15x, 'ID', 'F_MISS')
    vmiss_30x_dict = build_dict_in_chunks(tmp_vmiss_30x, 'ID', 'F_MISS')
    maf_dict = build_dict_in_chunks(tmp_maf, 'ID', 'ALT_FREQS')
    
    logger.info(f"Index dictionary construction complete: 15X={len(vmiss_15x_dict)}, 30X={len(vmiss_30x_dict)}, MAF={len(maf_dict)}")
    
    # Streaming processing: read, merge, write line by line
    logger.info("Starting streaming merge and writing results...")
    n_variants = 0
    
    with open(output_path, 'w') as outf:
        # Write header
        outf.write('VARIANT_ID\tVMISS\tVMISS_15X\tVMISS_30X\tMAF\n')
        
        # Read and process in chunks
        chunksize = 50000  # Process 50,000 rows at a time
        for chunk in pd.read_csv(tmp_vmiss_all, sep='\t', chunksize=chunksize):
            # Process current chunk
            for _, row in chunk.iterrows():
                variant_id = row['ID']
                vmiss = row['F_MISS']
                vmiss_15x = vmiss_15x_dict.get(variant_id, '')
                vmiss_30x = vmiss_30x_dict.get(variant_id, '')
                maf = maf_dict.get(variant_id, '')
                
                # Write one row
                outf.write(f"{variant_id}\t{vmiss}\t{vmiss_15x}\t{vmiss_30x}\t{maf}\n")
                n_variants += 1
            
            # Clean up after each chunk
            del chunk
            gc.collect()
            
            if n_variants % 100000 == 0:
                logger.info(f"Processed {n_variants} variants...")
    
    # Cleaning dictionaries and temporary files
    del vmiss_15x_dict, vmiss_30x_dict, maf_dict
    gc.collect()
    
    # Deleting temporary files
    for tmp_file in [tmp_vmiss_all, tmp_vmiss_15x, tmp_vmiss_30x, tmp_maf]:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
    
    logger.info(f"Streaming processing complete, wrote {n_variants} variants")
    
    return n_variants


# new function

def analyze_vmiss_thresholds(
    variant_metrics_file: str,
    output_dir: str,
    analysis_mode: Literal[1, 2, 3, 4] = 4,
    maf_thresholds: tuple = (0.01, 0.05),
    hist_step: float = 0.01,
    cdf_step: Optional[float] = None,
    cdf_step_15x: Optional[float] = None,
    cdf_step_30x: Optional[float] = None,
    bin_position: Literal['left', 'right'] = 'left',
    knee_curve: str = 'concave',
    knee_direction: str = 'increasing',
    knee_S: float = 1.0,
    knee_weight_x: float = 1.0,
    knee_weight_y: float = 1.0,
    threshold_decimals: int = 2,
    dpi: int = 600,
    chunksize: int = 50000,
    n_threads: int = 4
) -> dict:
    """
    Analyze variant missingness (VMISS) thresholds to find optimal filtering parameters
    
    Parameters:
        variant_metrics_file: Variant metrics file path (output from calculate_variant_metrics)
        output_dir: Output directory path
        analysis_mode: Analysis mode
            1: Identify VMISS for all variants (1 parameter recommendation)
            2: Identify VMISS_15X and VMISS_30X for all variants (2 parameter recommendations)
            3: Identify VMISS for three MAF groups (3 parameter recommendations)
            4: Identify VMISS_15X and VMISS_30X for three MAF groups (6 parameter recommendations)
        maf_thresholds: MAF grouping thresholds, default (0.01, 0.05)
        hist_step: Histogram bin step size, default 0.01 (fixed, for distribution plotting)
        cdf_step: CDF calculation step size, default None (auto-set based on mode)
            - Mode 1/3: defaults to 0.001
            - Mode 2/4: specified by cdf_step_15x and cdf_step_30x separately
        cdf_step_15x: CDF step size for 15X group, default 0.001 (only for modes 2/4)
        cdf_step_30x: CDF step size for 30X group, default 0.01 (only for modes 2/4)
        bin_position: CDF x-axis coordinate definition, default 'left'
            'left': Use bin left boundary (CDF represents ≤x cumulative proportion)
            'right': Use bin right boundary (CDF represents <x cumulative proportion, matches plink2 --geno behavior)
        knee_curve: kneed curve type, default 'concave'
        knee_direction: kneed direction, default 'increasing'
        knee_S: kneed sensitivity parameter
        knee_weight_x: kneed x-axis weight
        knee_weight_y: kneed y-axis weight
        threshold_decimals: Number of decimal places for threshold rounding (default: 2)
        dpi: Image resolution, default 600
        chunksize: Chunk size for streaming reads
        n_threads: Number of parallel threads
    
    Returns:
        dict: Dictionary containing all analysis results, including knee point info, image paths, statistics, etc.
    """
    import gc
    import json
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    from kneed import KneeLocator 

    n_threads = _normalize_thread_count(n_threads)

    if not os.path.exists(variant_metrics_file):
        raise FileNotFoundError(f"Variant metrics file not found: {variant_metrics_file}")
    if hist_step <= 0:
        raise ValueError(f"hist_step must be > 0, got: {hist_step}")
    if chunksize <= 0:
        raise ValueError(f"chunksize must be > 0, got: {chunksize}")
    if dpi <= 0:
        raise ValueError(f"dpi must be > 0, got: {dpi}")
    if bin_position not in ('left', 'right'):
        raise ValueError(f"bin_position must be 'left' or 'right', got: {bin_position}")

    maf_low, maf_mid = maf_thresholds
    if not (0 <= maf_low <= 1 and 0 <= maf_mid <= 1):
        raise ValueError(f"maf_thresholds must be in [0, 1], got: {maf_thresholds}")
    if maf_low >= maf_mid:
        raise ValueError(f"maf_thresholds must satisfy low < mid, got: {maf_thresholds}")
    
    logger.info("="*60)
    logger.info("Starting VMISS threshold analysis")
    logger.info(f"Input file: {variant_metrics_file}")
    logger.info(f"Analysis mode: {analysis_mode}")
    logger.info(f"MAF grouping thresholds: {maf_thresholds}")
    logger.info(f"Histogram step (hist_step): {hist_step}")
    logger.info(f"Number of threads: {n_threads}")
    
    # Auto-set CDF step size based on mode
    if analysis_mode in [1, 3]:
        # Mode 1 and 3: Don't separate 15X and 30X, use unified cdf_step
        if cdf_step is None:
            cdf_step = 0.001  # Default 0.001
        logger.info(f"CDF step: {cdf_step}")
        cdf_step_dict = {'default': cdf_step}
    elif analysis_mode in [2, 4]:
        # Mode 2 and 4: Separate 15X and 30X, set step sizes separately
        if cdf_step_15x is None:
            cdf_step_15x = 0.001  # 15X default 0.001
        if cdf_step_30x is None:
            cdf_step_30x = 0.001   # 30X default 0.001
        if cdf_step_15x <= 0 or cdf_step_30x <= 0:
            raise ValueError(f"cdf_step_15x and cdf_step_30x must be > 0, got: {cdf_step_15x}, {cdf_step_30x}")
        logger.info(f"CDF step (15X): {cdf_step_15x}")
        logger.info(f"CDF step (30X): {cdf_step_30x}")
        cdf_step_dict = {'15X': cdf_step_15x, '30X': cdf_step_30x}
    else:
        raise ValueError(f"Unsupported analysis_mode: {analysis_mode}, only 1/2/3/4 supported")
    
    logger.info("="*60)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Stream read data and group statistics (memory optimized)
    logger.info("Step 1: Streaming read and group data...")
    group_stats = _streaming_group_variants(
        variant_metrics_file, 
        analysis_mode, 
        maf_thresholds, 
        chunksize
    )
    
    # Step 2: Calculate cumulative distributions
    logger.info("Step 2: Calculate VMISS cumulative distributions...")
    cumulative_results = _calculate_cumulative_distribution(
        group_stats, 
        hist_step,
        cdf_step_dict,
        analysis_mode,
        bin_position
    )

    if not cumulative_results:
        raise ValueError("No valid cumulative distribution data generated. Please check input metrics and grouping parameters.")
    
    # Step 3: Find knee points in cumulative distributions
    logger.info("Step 3: Find cumulative distribution knee points...")
    knee_results = _find_cumulative_knee_points(
        cumulative_results,
        group_stats,
        knee_curve,
        knee_direction,
        knee_S,
        knee_weight_x,
        knee_weight_y,
        threshold_decimals
    )
    
    # Step 4: Plot distribution and cumulative distribution charts
    logger.info("Step 4: Plot distributions and cumulative distribution curves...")
    plot_path = _plot_vmiss_distributions(
        cumulative_results,
        knee_results,
        analysis_mode,
        output_dir,
        dpi,
        threshold_decimals
    )
    
    # Step 5: Save results to JSON (excluding temporary array data)
    logger.info("Step 5: Save results to JSON...")
    json_path = os.path.join(output_dir, 'vmiss_threshold_analysis.json')
    
    # Clean temporary arrays in group_stats (keys starting with _)
    group_stats_for_json = {}
    for key, value in group_stats.items():
        group_stats_for_json[key] = {k: v for k, v in value.items() if not k.startswith('_')}
    
    results = {
        'analysis_mode': analysis_mode,
        'maf_thresholds': maf_thresholds,
        'hist_step': hist_step,
        'cdf_step_dict': cdf_step_dict,
        'bin_position': bin_position,
        'knee_parameters': {
            'curve': knee_curve,
            'direction': knee_direction,
            'S': knee_S,
            'weight_x': knee_weight_x,
            'weight_y': knee_weight_y,
            'threshold_decimals': threshold_decimals
        },
        'group_statistics': group_stats_for_json,
        'knee_points': knee_results,
        'plot_path': plot_path,
        'json_path': json_path
    }
    
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    logger.info("="*60)
    logger.info(f"Analysis complete!")
    logger.info(f"Results JSON: {json_path}")
    logger.info(f"Analysis plot: {plot_path}")
    logger.info("="*60)
    
    # Print knee point recommendations
    _print_knee_recommendations(knee_results, analysis_mode)

    recommendation_tsv = os.path.join(output_dir, 'vmiss_threshold_recommendations.tsv')
    _write_threshold_recommendations_tsv(knee_results, recommendation_tsv)
    logger.info(f"Recommendation table: {recommendation_tsv}")
    
    return results


def _streaming_group_variants(
    variant_file: str,
    analysis_mode: int,
    maf_thresholds: tuple,
    chunksize: int
) -> dict:
    """
    Stream read variant data and gather group statistics
    Memory efficient: does not load all data at once
    """
    import gc
    
    # Initialize group statistics
    group_stats = {}
    
    # Determine required columns based on analysis mode
    required_cols = ['VARIANT_ID', 'MAF']
    vmiss_cols = []
    
    if analysis_mode in [1, 3]:
        vmiss_cols = ['VMISS']
    elif analysis_mode in [2, 4]:
        vmiss_cols = ['VMISS_15X', 'VMISS_30X']
    
    required_cols.extend(vmiss_cols)
    
    # Define MAF grouping (if needed)
    if analysis_mode in [3, 4]:
        maf_low, maf_mid = maf_thresholds
        maf_groups = {
            f'MAF<{maf_low}': lambda x: x < maf_low,
            f'MAF>={maf_low}&<={maf_mid}': lambda x: (x >= maf_low) & (x <= maf_mid),
            f'MAF>{maf_mid}': lambda x: x > maf_mid
        }
    else:
        maf_groups = {'All': lambda x: pd.Series([True] * len(x))}
    
    # Initialize data collectors for each group
    for maf_group in maf_groups:
        for vmiss_col in vmiss_cols:
            group_key = f"{maf_group}_{vmiss_col}" if len(maf_groups) > 1 else vmiss_col
            group_stats[group_key] = {
                'vmiss_values': [],  # Store VMISS values (for subsequent threshold analysis)
                'total_variants': 0,
                'vmiss_col': vmiss_col,
                'maf_group': maf_group
            }
    
    # Stream read and group
    logger.info(f"Starting streaming variant file read (chunksize={chunksize}）...")
    chunk_count = 0
    
    for chunk in pd.read_csv(variant_file, sep='\t', chunksize=chunksize, 
                              usecols=required_cols):
        chunk_count += 1
        
        # Convert MAF to numeric type, handle 'NA'
        chunk['MAF'] = pd.to_numeric(chunk['MAF'], errors='coerce')
        
        # Process each MAF group
        for maf_group_name, maf_filter in maf_groups.items():
            if maf_group_name == 'All':
                # For 'All' group, use entire chunk directly without filtering
                group_data = chunk
            else:
                # For specific MAF groups, use filter function
                maf_mask = maf_filter(chunk['MAF'])
                group_data = chunk[maf_mask]
            
            # Process each VMISS column
            for vmiss_col in vmiss_cols:
                group_key = f"{maf_group_name}_{vmiss_col}" if len(maf_groups) > 1 else vmiss_col
                
                # Convert VMISS to numeric, filter NA
                vmiss_series = pd.to_numeric(group_data[vmiss_col], errors='coerce')
                valid_vmiss = vmiss_series.dropna().values
                
                # Accumulate statistics
                group_stats[group_key]['vmiss_values'].extend(valid_vmiss.tolist())
                group_stats[group_key]['total_variants'] += len(valid_vmiss)
        
        # Periodically clean memory
        del chunk
        if chunk_count % 10 == 0:
            logger.info(f"Processed {chunk_count} data chunks...")
            gc.collect()
    
    # Convert to numpy arrays to accelerate subsequent calculations
    logger.info("Converting data to numpy arrays...")
    for group_key in group_stats:
        vmiss_list = group_stats[group_key]['vmiss_values']
        vmiss_arr = np.array(vmiss_list, dtype=np.float32)
        
        # Store arrays for subsequent calculations but don't save to final results
        group_stats[group_key]['_vmiss_array'] = vmiss_arr  # Underscore indicates temporary data
        # Delete lists to free memory
        del group_stats[group_key]['vmiss_values']
        
        # Calculate statistics
        if len(vmiss_arr) > 0:
            # Use actual min and max values without excessive rounding to preserve distribution details
            group_stats[group_key]['vmiss_min'] = float(vmiss_arr.min())
            group_stats[group_key]['vmiss_max'] = float(vmiss_arr.max())
            group_stats[group_key]['vmiss_mean'] = float(vmiss_arr.mean())
            group_stats[group_key]['vmiss_median'] = float(np.median(vmiss_arr))
        else:
            group_stats[group_key]['vmiss_min'] = 0
            group_stats[group_key]['vmiss_max'] = 0
            group_stats[group_key]['vmiss_mean'] = 0
            group_stats[group_key]['vmiss_median'] = 0
        
        logger.info(f"Group [{group_key}]: {group_stats[group_key]['total_variants']} variants, "
                   f"VMISS range [{group_stats[group_key]['vmiss_min']:.2f}, "
                   f"{group_stats[group_key]['vmiss_max']:.2f}]")
    
    gc.collect()
    return group_stats


def _calculate_cumulative_distribution(
    group_stats: dict,
    hist_step: float,
    cdf_step_dict: dict,
    analysis_mode: int,
    bin_position: str = 'left'
) -> dict:
    """
    Calculate distribution and cumulative distribution for each VMISS group
    
    Parameters:
        group_stats: Group statistics data
        hist_step: Histogram bin step size (fixed, used for plotting distributions)
        cdf_step_dict: CDF calculation step dictionary
            - Mode 1/3: {'default': 0.001}
            - Mode 2/4: {'15X': 0.001, '30X': 0.01}
        analysis_mode: Analysis mode (determines how to select cdf_step)
        bin_position: CDF x-axis coordinate definition method ('left', 'right')
    """
    cumulative_results = {}
    
    for group_key, stats in group_stats.items():
        logger.info(f"Calculating group [{group_key}] cumulative distribution...")
        
        vmiss_arr = stats['_vmiss_array']
        if len(vmiss_arr) == 0:
            logger.warning(f"Group [{group_key}] No valid data, skipping")
            continue
        
        vmiss_min = stats['vmiss_min']
        vmiss_max = stats['vmiss_max']
        
        # Determine CDF step size based on group name
        if analysis_mode in [1, 3]:
            # Mode 1/3: All groups use same step size
            cdf_step = cdf_step_dict['default']
        elif analysis_mode in [2, 4]:
            # Mode 2/4: Determine if 15X or 30X based on group name
            vmiss_col = stats.get('vmiss_col', '')
            if '15X' in vmiss_col.upper():
                cdf_step = cdf_step_dict['15X']
            elif '30X' in vmiss_col.upper():
                cdf_step = cdf_step_dict['30X']
            else:
                # Fallback: use default value if cannot determine
                cdf_step = cdf_step_dict.get('15X', 0.001)
                logger.warning(f"Cannot determine step for Group [{group_key}], using default value {cdf_step}")
        else:
            raise ValueError(f"Unsupported analysis_mode: {analysis_mode}, only 1/2/3/4 supported")
        
        logger.info(f"  Histogram step: {hist_step}, CDF step: {cdf_step}")
        
        # Align min and max to multiples of cdf_step (CDF calculation uses cdf_step)
        # Round down to nearest cdf_step multiple
        vmiss_min_aligned = np.floor(vmiss_min / cdf_step) * cdf_step
        # Round up to nearest cdf_step multiple
        vmiss_max_aligned = np.ceil(vmiss_max / cdf_step) * cdf_step
        
        # For right boundary, need to extend by one step to ensure last bin right edge includes all data
        if bin_position == 'right':
            vmiss_max_aligned += cdf_step
        
        # Ensure range contains at least 3 bins
        if vmiss_max_aligned - vmiss_min_aligned < cdf_step * 3:
            # Extend both sides to ensure at least 3 bins
            center = (vmiss_min_aligned + vmiss_max_aligned) / 2
            vmiss_min_aligned = center - cdf_step * 1.5
            vmiss_max_aligned = center + cdf_step * 1.5
        
        # Ensure within reasonable range (VMISS should be in [0, 1])
        vmiss_min_aligned = max(0, vmiss_min_aligned)
        vmiss_max_aligned = min(1, vmiss_max_aligned)
        
        # Generate VMISS value sequence (for CDF calculation, using cdf_step)
        # Use arange to ensure each bin width is cdf_step
        vmiss_bins = np.arange(vmiss_min_aligned, vmiss_max_aligned + cdf_step/2, cdf_step)
        vmiss_bins = np.round(vmiss_bins, 6)  # Increase precision to avoid floating point errors
        
        # Calculate histogram (VMISS distribution)
        hist_counts, _ = np.histogram(vmiss_arr, bins=vmiss_bins)
        
        # Calculate cumulative distribution (accumulate from small to large)
        cumulative_counts = np.cumsum(hist_counts)
        cumulative_percentage = (cumulative_counts / len(vmiss_arr)) * 100
        
        # Determine x-axis coordinates based on bin_position parameter
        if bin_position == 'left':
            # Left boundary: use left edge of bin (lower bound)
            # CDF meaning: proportion of variants ≤ x
            vmiss_x_values = vmiss_bins[:-1]
        elif bin_position == 'right':
            # Right boundary: use right edge of bin (upper bound)
            # CDF meaning: proportion of variants < x (consistent with plink2 --geno)
            vmiss_x_values = vmiss_bins[1:]
        else:
            raise ValueError(f"Invalid bin_position: {bin_position}. Must be 'left' or 'right'.")
        
        # Do not force starting point, directly use actual data CDF
        hist_counts_padded = hist_counts
        
        cumulative_results[group_key] = {
            'vmiss_values': vmiss_x_values.tolist(),  # x-axis：VMISS values (determined by bin_position)
            'histogram': hist_counts_padded.tolist(),     # Histogram counts
            'cumulative_counts': cumulative_counts.tolist(),  # Cumulative counts
            'cumulative_percentage': cumulative_percentage.tolist(),  # Cumulative percentage
            'total_variants': stats['total_variants'],
            'vmiss_col': stats['vmiss_col'],
            'maf_group': stats['maf_group'],
            'vmiss_min': vmiss_min,
            'vmiss_max': vmiss_max,
            'vmiss_mean': stats['vmiss_mean'],
            'vmiss_median': stats['vmiss_median'],
            'hist_step': hist_step,  # Record histogram step
            'cdf_step': cdf_step,    # Record CDF step
            'bin_position': bin_position,  # Record bin position definition method
            'n_bins': len(hist_counts)  # Record number of bins
        }
        
        logger.info(f"  Complete! Original range: [{vmiss_min:.4f}, {vmiss_max:.4f}], "
                   f"Aligned range: [{vmiss_min_aligned:.4f}, {vmiss_max_aligned:.4f}], "
                   f"{len(hist_counts)} bins (CDF step={cdf_step})")
    
    return cumulative_results


def _find_cumulative_knee_points(
    cumulative_results: dict,
    group_stats: dict,
    knee_curve: str,
    knee_direction: str,
    knee_S: float,
    knee_weight_x: float,
    knee_weight_y: float,
    threshold_decimals: int = 2
) -> dict:
    """
    Find knee points on cumulative distribution curves
    
    Optimization strategies：
    1. Adaptively adjust knee_S based on CDF step (smaller step size requires larger knee_S)
    2. Filter low cumulative percentage regions (avoid false positives at steep starts)
    3. Use stronger smoothing parameters to reduce noise
    4. Round threshold to specified decimal places and recalculate cumulative metrics
    """
    from kneed import KneeLocator
    
    knee_results = {}
    
    for group_key, result in cumulative_results.items():
        logger.info(f"Finding knee point for group [{group_key}] cumulative distribution...")
        
        vmiss_values = np.array(result['vmiss_values'])
        cumulative_pct = np.array(result['cumulative_percentage'])
        total = result['total_variants']
        cdf_step = result.get('cdf_step', 0.01)
        
        logger.info(f"  CDF data range: VMISS [{vmiss_values[0]:.3f}, {vmiss_values[-1]:.3f}], "
                   f"cumulative% [{cumulative_pct[0]:.1f}%, {cumulative_pct[-1]:.1f}%], "
                   f"{len(vmiss_values)} data points, CDF step={cdf_step}")
        
        # Check if data points are sufficient
        if len(vmiss_values) < 3:
            logger.warning(f"  Group [{group_key}] Insufficient data points (<3), skipping knee detection")
            knee_results[group_key] = {
                'knee_found': False,
                'message': 'Insufficient data points'
            }
            continue
        
        # Adaptively adjust knee_S: smaller CDF step requires larger knee_S to reduce noise
        adaptive_knee_S = knee_S
        if cdf_step <= 0.001:
            # Ultra-high precision: significantly increase knee_S
            adaptive_knee_S = knee_S * 3.0
            logger.info(f"  Detected ultra-high precision CDF (step={cdf_step}), adaptively increasing knee_S: {knee_S} → {adaptive_knee_S}")
        elif cdf_step <= 0.005:
            # High precision: moderately increase knee_S
            adaptive_knee_S = knee_S * 2.0
            logger.info(f"  Detected high precision CDF (step={cdf_step}), adaptively increasing knee_S: {knee_S} → {adaptive_knee_S}")
        
        # Strategy 1: Filter low cumulative percentage regions to avoid false positives at steep starts
        # Only keep data points with cumulative percentage >=5% for knee detection
        min_cumulative_threshold = 5.0  # Minimum cumulative percentage threshold
        valid_mask = cumulative_pct >= min_cumulative_threshold
        
        if valid_mask.sum() < 3:
            # If insufficient data points after filtering, lower threshold and retry
            min_cumulative_threshold = 1.0
            valid_mask = cumulative_pct >= min_cumulative_threshold
            logger.info(f"  Lowered cumulative threshold to {min_cumulative_threshold}% to preserve sufficient data points")
        
        if valid_mask.sum() < 3:
            logger.warning(f"  Group [{group_key}] Still insufficient data points after filtering, skipping knee detection")
            knee_results[group_key] = {
                'knee_found': False,
                'message': f'Insufficient data points after filtering (cumulative<{min_cumulative_threshold}%)'
            }
            continue
        
        # Use filtered data
        vmiss_filtered = vmiss_values[valid_mask]
        cumulative_pct_filtered = cumulative_pct[valid_mask]
        
        logger.info(f"  Filtered data: {len(vmiss_filtered)} points "
                   f"(VMISS range [{vmiss_filtered[0]:.3f}, {vmiss_filtered[-1]:.3f}], "
                   f"cumulative range [{cumulative_pct_filtered[0]:.1f}%, {cumulative_pct_filtered[-1]:.1f}%])")
        
        try:
            # Find knee point on filtered cumulative distribution curve
            # Use adaptive knee_S parameter
            kneedle = KneeLocator(
                vmiss_filtered,
                cumulative_pct_filtered,
                S=adaptive_knee_S,
                curve=knee_curve,
                direction=knee_direction,
                weight_x=knee_weight_x,
                weight_y=knee_weight_y
            )
            
            if kneedle.knee is not None:
                knee_vmiss_raw = kneedle.knee
                
                # Round threshold to specified decimal places
                knee_vmiss = round(knee_vmiss_raw, threshold_decimals)

                # Exact cumulative stats at rounded threshold: count variants with VMISS <= threshold
                # (not nearest-bin approximation on CDF grid).
                raw_vmiss = group_stats[group_key].get('_vmiss_array')
                if raw_vmiss is None:
                    raise KeyError(f"Missing _vmiss_array for group [{group_key}]")

                knee_cumulative_count = int(np.count_nonzero(raw_vmiss <= knee_vmiss))
                knee_cumulative_pct = float((knee_cumulative_count / total) * 100) if total > 0 else 0.0
                
                knee_results[group_key] = {
                    'knee_found': True,
                    'knee_vmiss': float(knee_vmiss),  # VMISS value at knee point (rounded)
                    'knee_vmiss_raw': float(knee_vmiss_raw),  # Original unrounded value
                    'knee_cumulative_percentage': knee_cumulative_pct,  # Cumulative percentage
                    'knee_cumulative_count': knee_cumulative_count,  # Cumulative variant count
                    'total_variants': total,
                    'vmiss_col': result['vmiss_col'],
                    'maf_group': result['maf_group'],
                    'adaptive_knee_S': adaptive_knee_S,  # Record used knee_S
                    'min_cumulative_filter': min_cumulative_threshold,  # Record filter threshold
                    'threshold_decimals': threshold_decimals  # Record decimal places used
                }
                
                logger.info(f"  ✓ Found knee point: VMISS={knee_vmiss_raw:.4f} → {knee_vmiss:.{threshold_decimals}f} (rounded), "
                           f"cumulative {knee_cumulative_pct:.1f}% ({knee_cumulative_count:,} variants), "
                           f"using knee_S={adaptive_knee_S:.1f}")
            else:
                logger.warning(f"  No obvious knee point found")
                knee_results[group_key] = {
                    'knee_found': False,
                    'message': 'No obvious knee point detected'
                }
        
        except Exception as e:
            logger.error(f"  Knee detection failed: {str(e)}")
            knee_results[group_key] = {
                'knee_found': False,
                'message': f'Detection failed: {str(e)}'
            }
    
    return knee_results


def _plot_vmiss_distributions(
    cumulative_results: dict,
    knee_results: dict,
    analysis_mode: int,
    output_dir: str,
    dpi: int,
    threshold_decimals: int = 2
) -> str:
    """
    Plot VMISS distribution charts and cumulative distribution curves (publication-quality).
    Each subplot contains:
    - Left Y-axis: Histogram (VMISS distribution) - always displayed with hist_step=0.01 width
    - Right Y-axis: Cumulative distribution curve - calculated using cdf_step precision
    - Mark knee point positions with publication-quality styling
    - Display threshold with specified decimal places
    """
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    import matplotlib as mpl
    
    # Clear all current figures to avoid memory leaks
    plt.close('all')
    
    # Publication-quality style configuration
    rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    rcParams['font.size'] = 14
    rcParams['axes.unicode_minus'] = False
    rcParams['axes.labelsize'] = 18
    rcParams['axes.titlesize'] = 22
    rcParams['xtick.labelsize'] = 14
    rcParams['ytick.labelsize'] = 14
    rcParams['legend.fontsize'] = 14
    rcParams['figure.titlesize'] = 26
    rcParams['axes.linewidth'] = 2.0
    rcParams['xtick.major.width'] = 2.0
    rcParams['ytick.major.width'] = 2.0
    rcParams['xtick.major.size'] = 8
    rcParams['ytick.major.size'] = 8
    rcParams['grid.linewidth'] = 1.0
    
    # Determine subplot layout
    n_groups = len(cumulative_results)
    if n_groups == 0:
        logger.warning("No data to plot")
        return ""
    
    # Determine layout based on analysis mode
    if analysis_mode == 1:
        nrows, ncols = 1, 1
    elif analysis_mode == 2:
        nrows, ncols = 1, 2
    elif analysis_mode == 3:
        nrows, ncols = 3, 1
    elif analysis_mode == 4:
        nrows, ncols = 3, 2
    else:
        raise ValueError(f"Unsupported analysis_mode: {analysis_mode}, only 1/2/3/4 supported")
    
    # Create figure with larger size for publication quality
    fig_width = 14 * ncols
    fig_height = 10 * nrows
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), dpi=dpi)
    fig.patch.set_facecolor('white')
    
    if n_groups == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Plot each group
    for idx, (group_key, result) in enumerate(cumulative_results.items()):
        ax1 = axes[idx]
        
        # Convert to numpy arrays to ensure correct plotting
        vmiss_values = np.array(result['vmiss_values'])  # CDF x-axis coordinates (cdf_step precision)
        histogram = np.array(result['histogram'])  # Histogram counts (cdf_step precision)
        cumulative_pct = np.array(result['cumulative_percentage'])
        total = result['total_variants']
        
        # Use stored step size parameters
        hist_step = result.get('hist_step', 0.01)  # Histogram step（Fixed at 0.01 for display）
        cdf_step = result.get('cdf_step', 0.001)   # CDF step（For calculation precision）
        
        # Critical fix: if cdf_step != hist_step, need to re-aggregate histogram data
        # This ensures visual consistency of histograms across all modes
        if abs(cdf_step - hist_step) > 1e-6:
            # cdf_step is finer, need to re-aggregate by hist_step
            vmiss_min = vmiss_values[0]
            vmiss_max = vmiss_values[-1]
            
            # Generate bin boundaries for hist_step
            hist_bins = np.arange(
                np.floor(vmiss_min / hist_step) * hist_step,
                np.ceil(vmiss_max / hist_step) * hist_step + hist_step/2,
                hist_step
            )
            
            # Aggregate fine-grained histogram to coarse-grained
            hist_aggregated = []
            hist_x_values = []
            
            for i in range(len(hist_bins) - 1):
                bin_left = hist_bins[i]
                bin_right = hist_bins[i + 1]
                
                # # Find all cdf data points in this bin range
                mask = (vmiss_values >= bin_left) & (vmiss_values < bin_right)
                bin_count = histogram[mask].sum()
                
                hist_aggregated.append(bin_count)
                hist_x_values.append(bin_left)
            
            # Use aggregated data for plotting
            hist_x_values = np.array(hist_x_values)
            hist_aggregated = np.array(hist_aggregated)
            
            logger.info(f"  [{group_key}] Re-aggregate histogram: "
                       f"{len(histogram)}points(step{cdf_step}) → {len(hist_aggregated)}points(step{hist_step})")
        else:
            # cdf_step == hist_step, use original data directly
            hist_x_values = vmiss_values
            hist_aggregated = histogram
        
        # Left Y-axis: plot histogram (distribution)
        # Now histogram data and bar width are based on hist_step, visually consistent
        color_hist = '#2E86AB'  # Professional blue
        ax1.bar(hist_x_values, hist_aggregated, width=hist_step * 0.8,
               color=color_hist, alpha=0.7, label='Distribution', edgecolor='white', linewidth=0.5)
        ax1.set_xlabel('Variant Missingness (VMISS)', fontsize=18, fontweight='bold', labelpad=14)
        ax1.set_ylabel('Variant Count', fontsize=18, fontweight='bold', color=color_hist, labelpad=14)
        ax1.tick_params(axis='y', labelcolor=color_hist, labelsize=14)
        ax1.tick_params(axis='x', labelsize=14)
        
        # Format y-axis with thousand separators
        from matplotlib.ticker import FuncFormatter
        ax1.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x):,}'))
        
        # Right Y-axis: plot cumulative distribution curve
        ax2 = ax1.twinx()
        color_cum = '#A23B72'  # Professional purple - Opaque for clarity
        # Use simple line plot without dense markers to ensure continuity
        # Only add markers at intervals if needed, but a solid line is standard for CDF
        ax2.plot(vmiss_values, cumulative_pct, color=color_cum, linewidth=4.0, 
                label='Cumulative %', linestyle='-', alpha=1.0)
        ax2.set_ylabel('Cumulative Percentage (%)', fontsize=18, fontweight='bold', 
                      color=color_cum, labelpad=14)
        ax2.tick_params(axis='y', labelcolor=color_cum, labelsize=14)
        ax2.set_ylim([0, 105])
        
        # Mark elbow point (academic style)
        if group_key in knee_results and knee_results[group_key].get('knee_found'):
            knee_info = knee_results[group_key]
            knee_vmiss = knee_info['knee_vmiss']
            knee_cum_pct = knee_info['knee_cumulative_percentage']
            knee_cum_count = knee_info['knee_cumulative_count']
            
            # Mark elbow point on cumulative distribution curve (solid circle, academic style)
            ax2.scatter([knee_vmiss], [knee_cum_pct], c='#C41E3A', marker='D',  # Red diamond
                       s=350, zorder=10, label='Optimal Threshold', 
                       edgecolors='white', linewidths=2.8)
            ax2.axvline(knee_vmiss, color='#C41E3A', linestyle='--', linewidth=2.5, 
                       alpha=0.8, zorder=5)
            
            # Add annotation (academic style: concise white background)
            # Use dynamic formatting based on threshold_decimals
            annotation_text = (f"Threshold = {knee_vmiss:.{threshold_decimals}f}\n"
                             f"Cumulative = {knee_cum_pct:.2f}%\n"
                             f"Retained N = {knee_cum_count:,}")
            
            ax2.annotate(annotation_text,
                        xy=(knee_vmiss, knee_cum_pct),
                        xytext=(30, -50),
                        textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=1.0', fc='white', ec='#C41E3A', 
                                 alpha=1.0, linewidth=2.2),
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.25',
                                      color='#C41E3A', lw=2.2),
                        fontsize=14,
                        fontweight='bold',
                        ha='left')
        
        # Set title with improved styling
        title = f"{result['maf_group']} — {result['vmiss_col']}"
        ax1.set_title(title, fontsize=22, fontweight='bold', pad=20)
        ax1.grid(True, alpha=0.35, axis='both', linestyle='--', linewidth=1.0)
        ax1.set_axisbelow(True)
        
        # Improve spine visibility
        for spine in ax1.spines.values():
            spine.set_linewidth(2.0)
        for spine in ax2.spines.values():
            spine.set_linewidth(2.0)
        
        # Combine legends, place in lower right corner with better styling
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='lower right', 
                  fontsize=14, framealpha=0.98, edgecolor='black', 
                  fancybox=True, shadow=True, frameon=True, borderpad=1.5)
    
    # Hide extra subplots
    for idx in range(n_groups, len(axes)):
        axes[idx].set_visible(False)
    
    # Adjust layout with enhanced spacing for professional appearance
    plt.tight_layout(pad=3.0, h_pad=2.5, w_pad=2.5)
    
    # Save figure with publication-quality settings
    plot_filename = f'vmiss_distribution_analysis_mode{analysis_mode}.png'
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path, dpi=max(dpi, 300), bbox_inches='tight', 
                facecolor='white', edgecolor='none', pad_inches=0.3)
    
    # Close figure object and release memory
    plt.close(fig)
    plt.close('all')
    
    logger.info(f"Publication-quality distribution plot saved: {plot_path} (DPI: {max(dpi, 300)})")
    return plot_path


def _print_knee_recommendations(knee_results: dict, analysis_mode: int):
    """
    Print VMISS threshold recommendations summary based on cumulative distribution elbow points
    """
    logger.info("\n" + "="*60)


def _write_threshold_recommendations_tsv(knee_results: dict, output_path: str):
    """Write knee-point recommendations to a machine-friendly TSV file."""
    rows = []
    for group_key, result in knee_results.items():
        rows.append({
            'group_key': group_key,
            'knee_found': bool(result.get('knee_found', False)),
            'recommended_vmiss': result.get('knee_vmiss', 'NA'),
            'knee_cumulative_percentage': result.get('knee_cumulative_percentage', 'NA'),
            'knee_cumulative_count': result.get('knee_cumulative_count', 'NA'),
            'total_variants': result.get('total_variants', 'NA'),
            'message': result.get('message', '')
        })

    pd.DataFrame(rows).to_csv(output_path, sep='\t', index=False)
    logger.info("VMISS threshold recommendations summary (based on cumulative distribution elbow points)")
    logger.info("="*60)
    
    for group_key, result in knee_results.items():
        if result.get('knee_found'):
            logger.info(f"\nGroup: {group_key}")
            logger.info(f"  Recommended VMISS threshold: {result['knee_vmiss']:.3f}")
            logger.info(f"  Cumulative percentage: {result['knee_cumulative_percentage']:.1f}%")
            logger.info(f"  Cumulative variant count: {result['knee_cumulative_count']:,}")
            logger.info(f"  Total variant count: {result['total_variants']:,}")
            logger.info(f"  Explanation: at this VMISS value, a clear elbow point appears in the cumulative distribution curve")
        else:
            logger.info(f"\nGroup: {group_key}")
            logger.info(f"  Status: {result.get('message', 'No knee point found')}")
    
    logger.info("\n" + "="*60)