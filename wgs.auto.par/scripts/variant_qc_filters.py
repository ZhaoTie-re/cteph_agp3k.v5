"""Filtering and visualization utilities for variant-level QC.

This module provides:
- VMISS distribution plotting and pass-list export
- HWE scatter/KDE visualization by MAF category
- Intersection extraction for VMISS/HWE pass variants
- Optional PLINK2 subset extraction from intersected variant IDs
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import logging
from matplotlib.patches import Patch

plt.style.use('default')

LOGGER = logging.getLogger("VariantQC.Filters")


def _log_info(phase: str, message: str) -> None:
    LOGGER.info(f"[{phase}] {message}")


def _log_warn(phase: str, message: str) -> None:
    LOGGER.warning(f"[{phase}] {message}")

def plot_vmiss_distribution(
    variant_qc_summary: str,
    vmiss_threshold: float = 0.05,
    variant_id_col: str = "VARIANT_ID",
    output_tsv: str = "vmiss_pass_variants.tsv",
    output_prefix: str = "cteph_agp3k"
) -> str:
    """Plot VMISS distribution and export variants passing the threshold."""
    data_chunks = []
    total_na_rows = 0

    cols_needed = ["VMISS", variant_id_col]
    reader = pd.read_csv(variant_qc_summary, sep="	", usecols=cols_needed, chunksize=500000)

    for i, chunk in enumerate(reader):
        before = len(chunk)
        chunk = chunk[chunk["VMISS"].notna()]
        na_dropped = before - len(chunk)
        total_na_rows += na_dropped
        if na_dropped > 0:
            _log_warn("PHASE-2", f"Chunk {i}: dropped {na_dropped} rows with NA VMISS")
        data_chunks.append(chunk)

    df = pd.concat(data_chunks, ignore_index=True)

    x_min, x_max = 0.0, 1.0
    bin_edges = np.arange(x_min, x_max + 0.02, 0.02)

    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(8, 6))

    total_variant_count = len(df)
    total_below_threshold = (df["VMISS"] < vmiss_threshold).sum()
    percent_pass = (100 * total_below_threshold / total_variant_count) if total_variant_count > 0 else 0

    passed_variants = df[df["VMISS"] < vmiss_threshold][variant_id_col].dropna().unique()
    pd.Series(passed_variants).to_csv(output_tsv, sep="	", index=False, header=False)

    sns.histplot(df["VMISS"], bins=bin_edges, kde=True, color="#1f77b4", edgecolor="black", ax=ax)
    ax.axvline(vmiss_threshold, linestyle="--", color="red", label=f"Threshold = {vmiss_threshold}")
    ax.set_xlim(x_min, x_max)
    ax.set_title("VMISS Distribution")
    ax.set_xlabel("VMISS")

    summary_label = (
        f"Total: {total_variant_count:,}\n"
        f"Pass: {total_below_threshold:,} ({percent_pass:.1f}%)"
    )
    dummy_patch = Patch(color='none', label=summary_label)

    handles, labels = ax.get_legend_handles_labels()
    handles.append(dummy_patch)
    ax.legend(handles=handles, loc='upper right', frameon=True, fontsize=10)

    fig.suptitle(f"VMISS Distribution\nPass Threshold (< {vmiss_threshold})", fontsize=16)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.vmiss.png", dpi=600)
    plt.show()
    _log_info("PHASE-2", f"VMISS plot written: {output_prefix}.vmiss.png")
    _log_info(
        "PHASE-2",
        f"VMISS summary: threshold={vmiss_threshold}, total_variants={total_variant_count:,}, "
        f"pass_variants={total_below_threshold:,} ({percent_pass:.1f}%)"
    )

    if total_na_rows > 0:
        _log_warn("PHASE-2", f"Total rows dropped due to NA VMISS: {total_na_rows:,}")

    return output_tsv


def plot_hwe_scatter_by_maf_category(
    variant_qc_summary: str,
    hwe_thresholds: dict = None,  # type: ignore
    maf_col: str = "CTRL_MAF",
    maf_threshold: float = 0.01,
    variant_id_col: str = "VARIANT_ID",
    output_tsv: str = "hwe_pass_variants.tsv",
    output_prefix: str = "cteph_agp3k",
    pass_vmiss_path: str | None = None,
) -> str:
    """Plot CTRL_HWE vs CASE_HWE by MAF bin and export pass variants.

    `hwe_thresholds` accepts per-bin thresholds such as:
    {
      "MAF < 0.01": {"CTRL_HWE": 1e-6, "CASE_HWE": 1e-6},
      "MAF >= 0.01": {"CTRL_HWE": 1e-6, "CASE_HWE": None}
    }
    """
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    import gc

    # Read header once to determine whether PASS_VMISS has been pre-encoded
    # in the summary table. When present, we avoid loading a large external
    # VMISS pass list and simply filter on this boolean column.
    header_df = pd.read_csv(variant_qc_summary, sep="\t", nrows=0)
    has_pass_vmiss = "PASS_VMISS" in header_df.columns

    # Optionally restrict to variants that passed VMISS QC.
    pass_vmiss_ids: set[str] | None = None
    if pass_vmiss_path is not None and not has_pass_vmiss:
        pass_vmiss_ids = set()
        with open(pass_vmiss_path) as f:
            for line in f:
                vid = line.strip()
                if vid:
                    pass_vmiss_ids.add(vid)
        _log_info(
            "PHASE-3",
            f"Restricting HWE analysis to VMISS-passed variants using external list: {len(pass_vmiss_ids):,} IDs from {pass_vmiss_path}"
        )
    elif has_pass_vmiss:
        _log_info("PHASE-3", "Restricting HWE analysis using PASS_VMISS flag encoded in variant_qc_summary")

    # Read required columns only, using chunking for memory efficiency.
    usecols = [variant_id_col, maf_col, "CTRL_HWE", "CASE_HWE"]
    if has_pass_vmiss:
        usecols.append("PASS_VMISS")

    def classify_maf(maf):
        if pd.isna(maf) or maf < 0:
            return None
        if maf < maf_threshold:
            return f"MAF < {maf_threshold}"
        else:
            return f"MAF >= {maf_threshold}"

    data_chunks: list[pd.DataFrame] = []
    reader = pd.read_csv(variant_qc_summary, sep="\t", usecols=usecols, chunksize=500000)

    raw_total = 0
    dropped_not_vmiss = 0
    dropped_invalid_maf = 0

    for i, chunk in enumerate(reader):
        raw_total += len(chunk)

        # First, restrict to VMISS-passed variants.
        if has_pass_vmiss:
            before_vmiss = len(chunk)
            # Use .loc and .copy() to avoid SettingWithCopyWarning when
            # assigning new columns on this filtered view.
            chunk = chunk.loc[chunk["PASS_VMISS"]].copy()
            removed_vmiss = before_vmiss - len(chunk)
            if removed_vmiss > 0:
                dropped_not_vmiss += removed_vmiss
                _log_info(
                    "PHASE-3",
                    f"Chunk {i}: removed {removed_vmiss} variants failing VMISS (PASS_VMISS == False)"
                )
        elif pass_vmiss_ids is not None:
            before_vmiss = len(chunk)
            # Same here: ensure we operate on a fresh copy after filtering.
            chunk = chunk.loc[chunk[variant_id_col].isin(pass_vmiss_ids)].copy()
            removed_vmiss = before_vmiss - len(chunk)
            if removed_vmiss > 0:
                dropped_not_vmiss += removed_vmiss
                _log_info(
                    "PHASE-3",
                    f"Chunk {i}: removed {removed_vmiss} variants failing VMISS (not in pass_vmiss list)"
                )

        # Then classify by MAF and drop only rows that cannot be assigned to a bin.
        chunk["Category"] = chunk[maf_col].apply(classify_maf)
        before_maf = len(chunk)
        chunk = chunk[chunk["Category"].notna()]
        removed_maf = before_maf - len(chunk)
        if removed_maf > 0:
            dropped_invalid_maf += removed_maf
            _log_warn(
                "PHASE-3",
                f"Chunk {i}: dropped {removed_maf} variants with invalid {maf_col} for MAF binning"
            )

        if len(chunk) == 0:
            del chunk
            gc.collect()
            continue

        data_chunks.append(chunk)
        del chunk
        gc.collect()

    if data_chunks:
        df = pd.concat(data_chunks, ignore_index=True)
    else:
        df = pd.DataFrame(columns=usecols + ["Category"])  # type: ignore[assignment]

    _log_info(
        "PHASE-3",
        f"HWE input loaded from {variant_qc_summary}: raw_total_variants={raw_total:,}, "
        f"after_vmiss_filter={raw_total - dropped_not_vmiss:,} (removed_by_vmiss={dropped_not_vmiss:,}), "
        f"used_for_HWE={len(df):,} (dropped_invalid_maf={dropped_invalid_maf:,}), "
        f"maf_col={maf_col}, maf_threshold={maf_threshold}"
    )

    pass_variant_ids = []

    categories = [
        f"MAF < {maf_threshold}",
        f"MAF >= {maf_threshold}"
    ]
    colors = ["#1f77b4", "#ff7f0e"]

    def plot_one_category(ax_main, ax_top, ax_right, sub_df, label, color, threshold_control, threshold_case):
        # Apply pass/fail thresholds on raw values; clean data only for plotting.
        x = sub_df["CTRL_HWE"]
        y = sub_df["CASE_HWE"]

        # Compute pass counts and quadrant statistics from raw values.
        if threshold_control is not None and threshold_case is not None:
            q1 = ((x >= threshold_control) & (y >= threshold_case)).sum()
            pass_ids = sub_df[(x >= threshold_control) & (y >= threshold_case)].index
            pass_variant_ids.extend(pass_ids)
            q2 = ((x < threshold_control) & (y >= threshold_case)).sum()
            q3 = ((x < threshold_control) & (y < threshold_case)).sum()
            q4 = ((x >= threshold_control) & (y < threshold_case)).sum()
        elif threshold_control is not None and threshold_case is None:
            q1 = (x >= threshold_control).sum()
            pass_ids = sub_df[(x >= threshold_control)].index
            pass_variant_ids.extend(pass_ids)
            q2 = (x < threshold_control).sum()
            q3 = 0
            q4 = 0
        else:
            q1 = len(sub_df)
            pass_ids = sub_df.index
            pass_variant_ids.extend(pass_ids)
            q2 = q3 = q4 = 0

        _log_info(
            "PHASE-3",
            f"HWE category={label}, total_variants={len(sub_df):,}, "
            f"pass_variants={len(pass_ids):,}, "
            f"CTRL_HWE_threshold={threshold_control if threshold_control is not None else 'None'}, "
            f"CASE_HWE_threshold={threshold_case if threshold_case is not None else 'None'}"
        )

        # ---- Scatter view (remove non-positive values for log scales) ----
        x_scatter = pd.to_numeric(x, errors='coerce')
        y_scatter = pd.to_numeric(y, errors='coerce')
        valid_scatter = np.isfinite(x_scatter) & np.isfinite(y_scatter) & (x_scatter > 0) & (y_scatter > 0)
        ax_main.scatter(x_scatter[valid_scatter], y_scatter[valid_scatter], alpha=0.3, c=color, s=20)
        if threshold_control is not None:
            ax_main.axvline(x=threshold_control, color='red', linestyle='--')
        if threshold_case is not None:
            ax_main.axhline(y=threshold_case, color='blue', linestyle='--')
        ax_main.set_xscale('log')
        ax_main.set_yscale('log')
        ax_main.set_xlabel('CTRL_HWE')
        ax_main.set_ylabel('CASE_HWE')
        ax_main.grid(True)

        # ---- Top KDE (CTRL_HWE) ----
        x_plot = pd.to_numeric(x, errors='coerce')
        x_plot = x_plot[np.isfinite(x_plot) & (x_plot > 0)]
        if x_plot.nunique() >= 2:
            sns.kdeplot(x_plot, ax=ax_top, fill=True, color='gray', linewidth=1.5, cut=0)
            ax_top.set_xscale('log')
            if threshold_control is not None:
                ax_top.axvline(x=threshold_control, color='red', linestyle='--')
        else:
            ax_top.text(0.5, 0.5, 'KDE skipped\n(non-positive or low variance)',
                        ha='center', va='center', transform=ax_top.transAxes, fontsize=8)
        ax_top.set_xlabel('')
        ax_top.set_ylabel('')
        ax_top.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                           labelbottom=False, labelleft=False)
        ax_top.grid(False)

        # ---- Right KDE (CASE_HWE) ----
        y_plot = pd.to_numeric(y, errors='coerce')
        y_plot = y_plot[np.isfinite(y_plot) & (y_plot > 0)]
        if y_plot.nunique() >= 2:
            sns.kdeplot(y=y_plot, ax=ax_right, fill=True, color='gray', linewidth=1.5, cut=0)
            ax_right.set_yscale('log')
            if threshold_case is not None:
                ax_right.axhline(y=threshold_case, color='blue', linestyle='--')
        else:
            ax_right.text(0.5, 0.5, 'KDE skipped\n(non-positive or low variance)',
                          ha='center', va='center', transform=ax_right.transAxes, fontsize=8)
        ax_right.set_xlabel('')
        ax_right.set_ylabel('')
        ax_right.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                             labelbottom=False, labelleft=False)
        ax_right.grid(False)
        # Inset panel for quadrant counts / totals.
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

        if threshold_control is not None and threshold_case is not None:
            # Draw quadrant separators and counts.
            inset_ax.axvline(x=xmid, linestyle='--', color='red', linewidth=1.5)
            inset_ax.axhline(y=ymid, linestyle='--', color='blue', linewidth=1.5)
            inset_ax.text((xmid + right_xlim[1]) / 2, (ymid + top_ylim[1]) / 2, f'{q1:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((right_xlim[0] + xmid) / 2, (ymid + top_ylim[1]) / 2, f'{q2:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((right_xlim[0] + xmid) / 2, (top_ylim[0] + ymid) / 2, f'{q3:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((xmid + right_xlim[1]) / 2, (top_ylim[0] + ymid) / 2, f'{q4:,}', ha='center', va='center', fontsize=7, fontweight='bold')
        elif threshold_control is not None and threshold_case is None:
            # Draw only CTRL_HWE separator and left/right counts.
            inset_ax.axvline(x=xmid, linestyle='--', color='red', linewidth=1.5)
            left = q2
            right = q1
            inset_ax.text((right_xlim[0] + xmid) / 2, ymid, f'{left:,}', ha='center', va='center', fontsize=7, fontweight='bold')
            inset_ax.text((xmid + right_xlim[1]) / 2, ymid, f'{right:,}', ha='center', va='center', fontsize=7, fontweight='bold')
        else:
            # No threshold configured: show total count only.
            inset_ax.text(xmid, ymid, f'{q1:,}', ha='center', va='center', fontsize=7, fontweight='bold')

    # Build 1x3 layout: two plots + one summary table panel.
    fig = plt.figure(figsize=(21, 7))
    outer_gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
    subplots = [outer_gs[0, 0], outer_gs[0, 1]]

    summary_rows = []

    # Render each MAF bin.
    for subplot_spec, category, color in zip(subplots, categories, colors):
        sub_df = df[df["Category"] == category]
        threshold = hwe_thresholds.get(category, {}) if hwe_thresholds else {}
        threshold_control = threshold.get("CTRL_HWE")
        threshold_case = threshold.get("CASE_HWE")

        # Nested layout: main scatter + top/right KDE panels.
        inner_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=subplot_spec,
                                                    width_ratios=[6, 1.5], height_ratios=[1.5, 6],
                                                    wspace=0.05, hspace=0.05)
        ax_main = plt.Subplot(fig, inner_gs[1, 0]) # type: ignore
        ax_top = plt.Subplot(fig, inner_gs[0, 0], sharex=ax_main) # type: ignore
        ax_right = plt.Subplot(fig, inner_gs[1, 1], sharey=ax_main) # type: ignore
        fig.add_subplot(ax_main)
        fig.add_subplot(ax_top)
        fig.add_subplot(ax_right)

        # Summary values used in the right-side table.
        summary_rows.append({
            f"{maf_col} Category": category.replace("MAF", maf_col).replace(" (", "\n("),
            "raw_category": category,
            "CTRL_HWE Threshold": f"{threshold_control:.1e}" if threshold_control is not None else "None",
            "CASE_HWE Threshold": f"{threshold_case:.1e}" if threshold_case is not None else "None",
        })

        plot_one_category(ax_main, ax_top, ax_right, sub_df, category, color, threshold_control, threshold_case)

    # Add summary table on the right panel.
    ax_legend = fig.add_subplot(outer_gs[0, 2])
    ax_legend.axis('off')
    summary_df = pd.DataFrame(summary_rows)
    summary_df_display = summary_df.drop(columns=["raw_category"])
    table = ax_legend.table(cellText=summary_df_display.values, # type: ignore
                            colLabels=summary_df_display.columns, # type: ignore
                            cellLoc='center',
                            colWidths=[0.35, 0.325, 0.325],
                            loc='center')
    table.scale(1.2, 1.6)
    table.auto_set_font_size(False)
    table.set_fontsize(9)

    # Style table header.
    for col_idx in range(len(summary_df_display.columns)):
        cell = table[(0, col_idx)]
        cell.set_facecolor('#f0f0f0')
        cell.set_text_props(color='black', weight='bold')

    # Highlight category rows by color.
    for row_idx, raw_category in enumerate(summary_df["raw_category"]):
        color = colors[categories.index(raw_category)]
        facecolor_rgba = plt.matplotlib.colors.to_rgba(color, alpha=0.35) # type: ignore
        cell = table[(row_idx + 1, 0)]
        cell.set_facecolor(facecolor_rgba)
        cell.set_text_props(color='black', weight='bold')

    fig.suptitle(f"CTRL_HWE vs CASE_HWE by {maf_col} Category", fontsize=16)
    
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="This figure includes Axes that are not compatible with tight_layout")
        plt.tight_layout(rect=[0, 0, 1, 0.96]) # type: ignore

    plt.savefig(f"{output_prefix}.hwe.png", dpi=600)
    plt.show()
    _log_info("PHASE-3", f"HWE plot written: {output_prefix}.hwe.png")

    # Export passed variant IDs when variant identifier column exists.
    if variant_id_col in df.columns:
        pass_variants = df.loc[pass_variant_ids, variant_id_col].dropna().unique()
        pd.Series(pass_variants).to_csv(output_tsv, sep="\t", index=False, header=False)
        _log_info(
            "PHASE-3",
            f"HWE pass list exported: {output_tsv} (n={len(pass_variants):,})"
        )
    else:
        _log_warn("PHASE-3", f"Column not found, skipping pass export: {variant_id_col}")
    
    return output_tsv


def extract_pass_variants_by_intersection(
    pass_vmiss_path,
    pass_hwe_path,
    output_path="pass_variants.tsv",
    bed_prefix=None,
    output_prefix="filtered",
    threads=4,
    plink2="plink2"
) -> str:
    """Finalize pass variants from the HWE pass list and export sorted IDs.

    In the current workflow, the HWE pass list is evaluated only on
    variants that already passed VMISS, so the final VMISS∩HWE set is
    semantically equivalent to the HWE pass list itself.

    This function:
    - Reads HWE-pass variant IDs;
    - Sorts by chromosome:position and writes the final pass list;
    - Optionally runs PLINK2 to extract a subset based on this list.
    """

    import os
    import subprocess

    def parse_variant_id(vid) -> tuple[str, int, str]:
        fields = str(vid).split(":")
        if len(fields) < 2:
            return ("", -1, vid)
        chr_str = fields[0].replace("chr", "")
        try:
            pos = int(fields[1])
        except Exception:
            pos = -1
        return (chr_str, pos, vid)

    chr_order = [str(i) for i in range(1, 23)]

    # Step 1: load HWE pass list; this is treated as the final VMISS∩HWE set.
    hwe_ids_list: list[str] = []
    with open(pass_hwe_path) as f:
        for line in f:
            vid = line.strip()
            if vid:
                hwe_ids_list.append(vid)
    hwe_ids_set = set(hwe_ids_list)
    _log_info(
        "PHASE-4",
        f"Loaded HWE pass list: {pass_hwe_path} (n={len(hwe_ids_list):,}); "
        "treating as final VMISS∩HWE pass set."
    )
    # Step 2: sort HWE IDs by chromosome and position and write final output.
    parsed = []
    for vid in hwe_ids_list:
        chr_str, pos, _ = parse_variant_id(vid)
        if chr_str in chr_order:
            parsed.append((chr_order.index(chr_str), pos, vid))

    parsed.sort()
    with open(output_path, "w") as f_out:
        for _, _, vid in parsed:
            f_out.write(f"{vid}\n")

    _log_info(
        "PHASE-4",
        f"Final pass list written: {output_path} (n={len(parsed):,})"
    )

    # Optional: run PLINK2 extraction on the final pass list.
    if bed_prefix is not None:
        plink_cmd = [
            plink2,
            "--bfile", bed_prefix,
            "--extract", output_path,
            "--make-bed",
            "--out", output_prefix,
            "--threads", str(threads)
        ]
        _log_info("PHASE-4", f"Running PLINK2 extraction for final pass list: {output_path}")
        subprocess.run(plink_cmd, check=True)

    return output_path


def stratify_by_maf(
    variant_qc_summary: str,
    pass_variants_path: str,
    maf_col: str = "CTRL_MAF",
    maf_threshold: float = 0.01,
    bed_prefix: str | None = None,
    output_prefix: str = "stratified",
    threads: int = 4,
    plink2: str = "plink2",
    variant_id_col: str = "VARIANT_ID"
) -> tuple[str, str]:
    """Stratify QC-passed variants into two MAF-based groups and export PLINK files.

    Output files are named dynamically based on MAF source and threshold:
    - Below threshold: {output_prefix}.{maf_source}_lt{threshold}.{bed,bim,fam,variants.tsv}
    - Above threshold: {output_prefix}.{maf_source}_ge{threshold}.{bed,bim,fam,variants.tsv}
    
    Example with CTRL_MAF and threshold 0.01:
    - output.ctrl_maf_lt001.bed (variants with CTRL_MAF < 0.01)
    - output.ctrl_maf_ge001.bed (variants with CTRL_MAF >= 0.01)

    Args:
        variant_qc_summary: Path to variant QC summary TSV file
        pass_variants_path: Path to pass variants list (one variant ID per line)
        maf_col: MAF column name to use for stratification (CTRL_MAF, CASE_MAF, or MAF)
        maf_threshold: MAF threshold for stratification
        bed_prefix: Prefix of input PLINK binary files
        output_prefix: Prefix for output files
        threads: Number of CPU threads
        plink2: Path to plink2 executable
        variant_id_col: Column name for variant IDs

    Returns:
        Tuple of (below_threshold_list_path, above_threshold_list_path)
    """
    import os
    import subprocess

    _log_info("PHASE-5", f"Stratifying variants by {maf_col} with threshold {maf_threshold}")

    # Derive MAF source identifier from column name
    maf_source_map = {
        "CTRL_MAF": "ctrl_maf",
        "CASE_MAF": "case_maf",
        "MAF": "all_maf"
    }
    maf_source = maf_source_map.get(maf_col, maf_col.lower().replace("_", ""))
    
    # Build dynamic suffix based on MAF source and threshold
    threshold_str = str(maf_threshold).replace(".", "")
    suffix_below = f"{maf_source}_lt{threshold_str}"
    suffix_above = f"{maf_source}_ge{threshold_str}"

    # Step 1: Load pass variants into a set for fast lookup
    pass_ids = set()
    with open(pass_variants_path) as f:
        for line in f:
            vid = line.strip()
            if vid:
                pass_ids.add(vid)
    _log_info("PHASE-5", f"Loaded {len(pass_ids):,} pass variants from {pass_variants_path}")

    # Step 2: Read variant QC summary and filter to pass variants only
    usecols = [variant_id_col, maf_col]
    df = pd.read_csv(variant_qc_summary, sep="\t", usecols=usecols)
    df = df[df[variant_id_col].isin(pass_ids)]
    _log_info("PHASE-5", f"Matched {len(df):,} variants in QC summary")

    # Step 3: Classify variants by MAF threshold
    df_below = df[df[maf_col] < maf_threshold]
    df_above = df[df[maf_col] >= maf_threshold]

    _log_info("PHASE-5", f"Below threshold ({maf_col} < {maf_threshold}): {len(df_below):,} variants")
    _log_info("PHASE-5", f"Above threshold ({maf_col} >= {maf_threshold}): {len(df_above):,} variants")

    # Step 4: Export variant ID lists
    below_list = f"{output_prefix}.{suffix_below}.variants.tsv"
    above_list = f"{output_prefix}.{suffix_above}.variants.tsv"

    df_below[variant_id_col].to_csv(below_list, sep="\t", index=False, header=False)
    df_above[variant_id_col].to_csv(above_list, sep="\t", index=False, header=False)

    _log_info("PHASE-5", f"Exported variant list ({maf_col} < {maf_threshold}): {below_list}")
    _log_info("PHASE-5", f"Exported variant list ({maf_col} >= {maf_threshold}): {above_list}")

    # Step 5: Run PLINK2 extraction for each group
    if bed_prefix is not None:
        # Below threshold genotypes
        plink_cmd_below = [
            plink2,
            "--bfile", bed_prefix,
            "--extract", below_list,
            "--make-bed",
            "--out", f"{output_prefix}.{suffix_below}",
            "--threads", str(threads)
        ]
        _log_info("PHASE-5", f"Extracting genotypes ({maf_col} < {maf_threshold}): {output_prefix}.{suffix_below}")
        subprocess.run(plink_cmd_below, check=True)

        # Above threshold genotypes
        plink_cmd_above = [
            plink2,
            "--bfile", bed_prefix,
            "--extract", above_list,
            "--make-bed",
            "--out", f"{output_prefix}.{suffix_above}",
            "--threads", str(threads)
        ]
        _log_info("PHASE-5", f"Extracting genotypes ({maf_col} >= {maf_threshold}): {output_prefix}.{suffix_above}")
        subprocess.run(plink_cmd_above, check=True)

        _log_info("PHASE-5", "MAF stratification completed successfully")

    return (below_list, above_list)
