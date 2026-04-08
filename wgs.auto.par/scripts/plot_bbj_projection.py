#!/usr/bin/env python3
"""Create publication-style projection PCA figures.

Outputs:
1) A single PNG summarizing eigenvalue variance, explained variance ratio,
   and cumulative explained variance.
2) Pairwise PC scatter plots (PC1 vs PC2, PC3 vs PC4, ...), each saved as
    a PNG and also collected into a multi-page PDF (one page per pair).

Data handling uses a single merged .sscore source:
- BBJ and case-control samples are scored together in one .sscore.
- Group labels are assigned by IID prefix (BBJ) and sample_info mapping.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot BBJ PCA variance summary and pairwise projection scatter plots."
    )
    parser.add_argument("--bbj-eigenval", required=True, help="Path to BBJ eigenval file")
    parser.add_argument("--projected-sscore", required=True, help="Path to projected .sscore file")
    parser.add_argument("--sample-info", required=True, help="Sample info table (xlsx/csv/tsv)")

    parser.add_argument("--sample-id-col", required=True, help="Sample ID column in sample info")
    parser.add_argument("--phenotype-col", required=True, help="Phenotype column in sample info")
    parser.add_argument("--phenotype-case-value", required=True, help="Case value in phenotype column")
    parser.add_argument("--phenotype-ctrl-value", required=True, help="Control value in phenotype column")

    parser.add_argument("--bbj-id-prefix", default="bbj_", help="Prefix of BBJ sample IDs")
    parser.add_argument("--bbj-label", default="BBJ", help="Label for BBJ samples")
    parser.add_argument("--case-label", default="CTEPH", help="Label for case samples")
    parser.add_argument("--ctrl-label", default="AGP3K", help="Label for control samples")

    parser.add_argument("--max-pcs", type=int, default=20, help="Number of PCs to use")
    parser.add_argument(
        "--keep-non-bbj-iids",
        default=None,
        help=(
            "Optional keep list file (FID IID without header, or IID-only). "
            "When provided, non-BBJ samples in projected .sscore are filtered "
            "to IIDs in this list; BBJ-prefix samples are always retained."
        ),
    )
    parser.add_argument("--out-prefix", required=True, help="Output prefix for all figures")
    return parser.parse_args()


def _read_table_auto(path: str) -> pd.DataFrame:
    lower = path.lower()
    if lower.endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    if lower.endswith(".tsv"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _sanitize_sample_info(sample_info: pd.DataFrame, sample_id_col: str) -> pd.DataFrame:
    """Drop rows with invalid sample IDs in sample_id_col.

    Invalid includes NA/NaN, empty strings, and literal 'nan' after stripping.
    """
    if sample_id_col not in sample_info.columns:
        raise ValueError(f"sample_id_col not found: {sample_id_col}")

    out = sample_info.copy()
    sid = out[sample_id_col]
    sid_str = sid.astype(str).str.strip()

    valid_mask = sid.notna() & sid_str.ne("") & sid_str.str.lower().ne("nan")
    dropped_n = int((~valid_mask).sum())
    if dropped_n > 0:
        print(f"[INFO] Dropped {dropped_n:,} sample_info rows with invalid {sample_id_col}")

    return out.loc[valid_mask].copy()


def _load_eigenval(path: str, max_pcs: int) -> np.ndarray:
    vals: list[float] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            for token in line.split():
                vals.append(float(token))
    if not vals:
        raise ValueError(f"No eigenvalues found: {path}")
    arr = np.asarray(vals, dtype=float)
    return arr[: min(max_pcs, arr.shape[0])]


def _load_projected_sscore(path: str, max_pcs: int) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", engine="python")

    iid_col = "IID" if "IID" in df.columns else None
    if iid_col is None:
        for cand in ["#IID", "ID", "SAMPLE", "IID1"]:
            if cand in df.columns:
                iid_col = cand
                break
    if iid_col is None:
        raise ValueError(f"Cannot find sample ID column in sscore: {path}")

    score_cols: list[tuple[int, str]] = []

    # Preferred patterns in order. The projection in this pipeline currently
    # yields PC1_AVG..PC20_AVG, while some PLINK settings can yield SCORE1_AVG...
    patterns = [
        re.compile(r"^PC(\d+)_AVG$"),
        re.compile(r"^SCORE(\d+)_AVG$"),
        re.compile(r"^PC(\d+)_SUM$"),
        re.compile(r"^SCORE(\d+)_SUM$"),
    ]

    for pat in patterns:
        score_cols.clear()
        for c in df.columns:
            m = pat.match(c)
            if m:
                score_cols.append((int(m.group(1)), c))
        if score_cols:
            break

    if not score_cols:
        raise ValueError(
            "No supported score columns found in sscore. "
            f"Expected PC*_AVG / SCORE*_AVG / PC*_SUM / SCORE*_SUM, file={path}"
        )

    score_cols.sort(key=lambda x: x[0])
    score_cols = score_cols[:max_pcs]

    out = pd.DataFrame({"IID": df[iid_col].astype(str)})
    for idx, col in score_cols:
        out[f"PC{idx}"] = pd.to_numeric(df[col], errors="coerce")

    required = [f"PC{i}" for i in range(1, len(score_cols) + 1)]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(f"Missing projected PC columns after parsing: {missing}")

    return out


def _load_keep_iids(path: str) -> set[str]:
    """Load keep IID set from a FID/IID file (no header) or IID-only file."""
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    if df.shape[1] >= 2:
        series = df.iloc[:, 1]
    else:
        series = df.iloc[:, 0]

    iids = (
        series.astype(str)
        .str.strip()
        .loc[lambda s: s.ne("") & s.str.lower().ne("nan")]
    )
    return set(iids.tolist())


def _build_group_map(
    sample_info: pd.DataFrame,
    sample_id_col: str,
    phenotype_col: str,
    case_value: str,
    ctrl_value: str,
    case_label: str,
    ctrl_label: str,
) -> dict[str, str]:
    if sample_id_col not in sample_info.columns:
        raise ValueError(f"sample_id_col not found: {sample_id_col}")
    if phenotype_col not in sample_info.columns:
        raise ValueError(f"phenotype_col not found: {phenotype_col}")

    sid = sample_info[sample_id_col].astype(str).str.strip()
    pheno = sample_info[phenotype_col].astype(str).str.strip()

    out: dict[str, str] = {}
    for s, p in zip(sid, pheno):
        if not s or s.lower() == "nan":
            continue
        if p == case_value:
            out[s] = case_label
        elif p == ctrl_value:
            out[s] = ctrl_label
    return out


def _configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 14,
            "axes.titlesize": 16,
            "axes.labelsize": 16,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 12,
            "axes.linewidth": 1.2,
            "grid.linewidth": 0.5,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _style_axes(ax: Axes) -> None:
    """Apply a clean, publication-style look to an axis."""

    for spine in ["top", "right"]:
        if spine in ax.spines:
            ax.spines[spine].set_visible(False)
    ax.tick_params(axis="both", which="both", direction="out", top=False, right=False)
    
    # Remove all background grids
    ax.grid(False)


def plot_variance_summary(eigenvals: np.ndarray, out_png: Path) -> tuple[np.ndarray, np.ndarray]:
    pcs = np.arange(1, len(eigenvals) + 1)
    explained = eigenvals / np.sum(eigenvals)
    explained_pct = explained * 100.0
    cumulative_pct = np.cumsum(explained_pct)

    # Single, publication-style panel:
    #  - Bars: explained variance (%), left y-axis
    #  - Line: cumulative variance (%), right y-axis
    fig, ax1 = plt.subplots(figsize=(8, 4.5))

    xtick_labels = [f"PC{i}" for i in pcs]

    # Explained variance (left y-axis) and cumulative variance (right y-axis)
    bars = ax1.bar(
        pcs,
        explained_pct,
        color="#55A868",
        edgecolor="black",
        linewidth=0.6,
        label="Explained variance (%)",
    )
    ax1.set_ylabel("Explained variance (%)")
    ax1.grid(axis="y", alpha=0.35, linestyle=":")
    ax1.set_xticks(pcs)
    ax1.set_xticklabels(xtick_labels, rotation=45, ha="right")
    ax1.set_ylim(bottom=0)

    ax2 = ax1.twinx()
    line = ax2.plot(
        pcs,
        cumulative_pct,
        marker="o",
        color="#C44E52",
        linewidth=2.0,
        markersize=4,
        label="Cumulative variance (%)",
    )[0]
    ax2.set_ylabel("Cumulative variance (%)", labelpad=18)
    ax2.set_ylim(0, 102)
    ax2.grid(False)
    _style_axes(ax1)
    _style_axes(ax2)

    # Combined legend for bars + line
    handles = [bars, line]
    labels = [h.get_label() for h in handles]
    ax1.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(1.20, 0.5),
        frameon=False,
        borderaxespad=0.0,
    )

    # Leave extra room on the right for the vertical y-label and legend
    fig.subplots_adjust(left=0.10, right=0.78, bottom=0.22, top=0.95)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

    return explained_pct, cumulative_pct


def plot_pairwise_scatter(
    plot_df: pd.DataFrame,
    explained_pct: np.ndarray,
    out_prefix: Path,
    bbj_label: str,
    case_label: str,
    ctrl_label: str,
) -> None:
    color_map = {
        bbj_label: "#B0B0B0",  # solid neutral grey
        ctrl_label: "#1F78B4",  # blue controls
        case_label: "#E31A1C",  # deep red cases
        "OTHER": "#33A02C",   # any unlabelled samples
    }
    order = [bbj_label, ctrl_label, case_label, "OTHER"]

    pdf_path = Path(str(out_prefix) + ".pc_pairs.pdf")

    pc_count = len(explained_pct)
    with PdfPages(pdf_path) as pdf:
        for i in range(1, pc_count + 1, 2):
            j = i + 1
            if j > pc_count:
                break
            x_col = f"PC{i}"
            y_col = f"PC{j}"

            fig, ax = plt.subplots(figsize=(8.5, 7.5))
            fig.subplots_adjust(left=0.12, right=0.78, bottom=0.15, top=0.9)
            _style_axes(ax)

            size_map = {
                bbj_label: 4.0,
                ctrl_label: 14.0,
                case_label: 18.0,
                "OTHER": 8.0,
            }
            alpha_map = {
                bbj_label: 0.60,
                ctrl_label: 0.80,
                case_label: 0.90,
                "OTHER": 0.60,
            }
            zorder_map = {
                bbj_label: 1,
                ctrl_label: 2,
                case_label: 3,
                "OTHER": 2,
            }

            group_counts: dict[str, int] = {}

            for group in order:
                sub = plot_df[plot_df["GROUP"] == group]
                if sub.empty:
                    continue
                group_counts[group] = len(sub)
                point_size = size_map.get(group, 8.0)
                point_alpha = alpha_map.get(group, 0.6)
                zorder = zorder_map.get(group, 2)
                ax.scatter(
                    sub[x_col],
                    sub[y_col],
                    s=point_size,
                    alpha=point_alpha,
                    c=color_map[group],
                    label=group,
                    edgecolors="none",
                    zorder=zorder,
                    rasterized=True,
                )

            ax.set_xlabel(f"PC{i} ({explained_pct[i-1]:.1f}%)")
            ax.set_ylabel(f"PC{j} ({explained_pct[j-1]:.1f}%)")

            # Add faint coordinate axes (x=0, y=0)
            ax.axhline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.5, zorder=0)
            ax.axvline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.5, zorder=0)
            
            # Add appropriate, slightly clearer background grid
            ax.grid(True, linestyle="--", alpha=0.75, linewidth=0.7, color="#C0C0C0", zorder=0)
            # Build a custom legend so that marker size is consistent
            legend_handles: list[Line2D] = []
            legend_labels: list[str] = []
            for group in order:
                n = group_counts.get(group, 0)
                if n == 0:
                    continue
                handle = Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    markersize=7.0,
                    markerfacecolor=color_map[group],
                    markeredgecolor="none",
                )
                legend_handles.append(handle)
                legend_labels.append(f"{group} (n={n:,})")

            # Place legend outside the main plotting area (right side)
            ax.legend(
                handles=legend_handles,
                labels=legend_labels,
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                frameon=False,
                borderaxespad=0.0,
            )

            # Matplotlib can save directly to PDF without intermediate PNGs.
            # We only save to the multi-page PDF here as requested.
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def main() -> None:
    args = parse_args()
    _configure_matplotlib()

    out_prefix = Path(args.out_prefix)

    eigenvals = _load_eigenval(args.bbj_eigenval, args.max_pcs)
    score_df = _load_projected_sscore(args.projected_sscore, len(eigenvals))
    sample_info_raw = _read_table_auto(args.sample_info)
    sample_info = _sanitize_sample_info(sample_info_raw, args.sample_id_col)

    group_map = _build_group_map(
        sample_info=sample_info,
        sample_id_col=args.sample_id_col,
        phenotype_col=args.phenotype_col,
        case_value=args.phenotype_case_value,
        ctrl_value=args.phenotype_ctrl_value,
        case_label=args.case_label,
        ctrl_label=args.ctrl_label,
    )

    # All samples are scored together in the merged .sscore.
    # BBJ labels are derived from IID prefix, while projected case/control
    # labels are derived from sample_info.
    iid_series = score_df["IID"].astype(str)
    is_bbj = iid_series.str.startswith(args.bbj_id_prefix)

    if args.keep_non_bbj_iids:
        keep_iids = _load_keep_iids(args.keep_non_bbj_iids)
        keep_mask = is_bbj | iid_series.isin(keep_iids)
        before_n = int(score_df.shape[0])
        score_df = score_df.loc[keep_mask].copy()
        after_n = int(score_df.shape[0])
        print(
            "[INFO] Applied --keep-non-bbj-iids: "
            f"kept {after_n:,}/{before_n:,} rows in projected sscore"
        )
        iid_series = score_df["IID"].astype(str)
        is_bbj = iid_series.str.startswith(args.bbj_id_prefix)

    score_df["SOURCE"] = np.where(is_bbj, "BBJ_MERGED", "PROJECTION")

    score_df["GROUP"] = iid_series.map(lambda x: group_map.get(str(x), "OTHER"))
    score_df.loc[is_bbj, "GROUP"] = args.bbj_label

    pc_cols = [f"PC{i}" for i in range(1, len(eigenvals) + 1)]
    plot_df = score_df[["IID", "SOURCE", "GROUP"] + pc_cols].copy()

    variance_png = Path(f"{out_prefix}.variance_summary.png")
    explained_pct, _ = plot_variance_summary(eigenvals, variance_png)

    plot_pairwise_scatter(
        plot_df=plot_df,
        explained_pct=explained_pct,
        out_prefix=out_prefix,
        bbj_label=args.bbj_label,
        case_label=args.case_label,
        ctrl_label=args.ctrl_label,
    )

    print(f"[OK] Variance summary PNG: {variance_png}")
    print(f"[OK] Pairwise PDF: {out_prefix}.pc_pairs.pdf")


if __name__ == "__main__":
    main()
