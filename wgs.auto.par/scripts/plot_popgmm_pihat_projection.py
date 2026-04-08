#!/usr/bin/env python3
"""Plot PCA projection results for PopGMM samples after relatedness-aware exclusion.

Outputs:
1) A PNG summarizing explained/cumulative variance from base PCA eigenvalues.
2) A PDF with pairwise PC scatter plots (PC1 vs PC2, PC3 vs PC4, ...).
3) A PNG with per-PC case/control KDE distributions.
4) A log file with per-PC summary statistics and group-comparison tests.

Compared with plot_bbj_projection.py:
- No BBJ group is present.
- Samples in the relatedness-flagged list are highlighted on top of the base case/control colors.
"""

from __future__ import annotations

import argparse
import logging
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, ttest_ind


CTRL_COLOR = "#1F78B4"
CASE_COLOR = "#E31A1C"
OTHER_COLOR = "#33A02C"
OVERLAY_FACE = "#BBBBBB"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot PopGMM PIHAT-relatedness projection (case/control with overlay highlights)."
    )
    parser.add_argument("--base-eigenval", required=True, help="Path to base PCA eigenval file")
    parser.add_argument("--projected-sscore", required=True, help="Path to all_samples_projection .sscore")
    parser.add_argument(
        "--relatedness-flagged-fid-iid",
        "--overlap-fid-iid",
        dest="relatedness_flagged_fid_iid",
        required=True,
        help="Path to relatedness-flagged FID/IID file (no header)",
    )
    parser.add_argument("--sample-info", required=True, help="Sample info table (xlsx/csv/tsv)")

    parser.add_argument("--sample-id-col", required=True, help="Sample ID column in sample info")
    parser.add_argument("--phenotype-col", required=True, help="Phenotype column in sample info")
    parser.add_argument("--phenotype-case-value", required=True, help="Case value in phenotype column")
    parser.add_argument("--phenotype-ctrl-value", required=True, help="Control value in phenotype column")

    parser.add_argument("--case-label", default="CTEPH", help="Label for case samples")
    parser.add_argument("--ctrl-label", default="AGP3K", help="Label for control samples")
    parser.add_argument("--max-pcs", type=int, default=20, help="Number of PCs to use")
    parser.add_argument("--kde-n-cols", type=int, default=5, help="Number of columns in the per-PC KDE grid")
    parser.add_argument("--out-prefix", required=True, help="Output prefix for figures")
    return parser.parse_args()


def _read_table_auto(path: str) -> pd.DataFrame:
    lower = path.lower()
    if lower.endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    if lower.endswith(".tsv"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _sanitize_sample_info(sample_info: pd.DataFrame, sample_id_col: str) -> pd.DataFrame:
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

    return out


def _load_overlap_iids(path: str) -> set[str]:
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


def _bh_fdr_adjust(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)
    mask = np.isfinite(pvals)
    if not mask.any():
        return out

    p = pvals[mask]
    n = p.size
    order = np.argsort(p)
    adj = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    tmp = np.empty_like(adj)
    tmp[order] = adj
    out[mask] = tmp
    return out


def _safe_stats(x: np.ndarray) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"mean": np.nan, "std": np.nan, "min": np.nan, "max": np.nan}
    return {
        "mean": float(np.mean(x)),
        "std": float(np.std(x, ddof=1)) if x.size >= 2 else float(np.std(x, ddof=0)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def _setup_file_logger(name: str, log_path: Path) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for handler in list(logger.handlers):
        logger.removeHandler(handler)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(sh)
    return logger


def _build_pc_stats_and_tests(
    plot_df: pd.DataFrame,
    pc_cols: list[str],
    case_label: str,
    ctrl_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    case_set = plot_df["GROUP"].eq(case_label)
    ctrl_set = plot_df["GROUP"].eq(ctrl_label)

    pc_stats: list[dict[str, object]] = []
    test_rows: list[dict[str, object]] = []

    for col in pc_cols:
        pc_num = int(str(col).replace("PC", ""))
        x_ctrl = pd.to_numeric(plot_df.loc[ctrl_set, col], errors="coerce").to_numpy(dtype=float, copy=False)
        x_case = pd.to_numeric(plot_df.loc[case_set, col], errors="coerce").to_numpy(dtype=float, copy=False)
        x_ctrl = x_ctrl[np.isfinite(x_ctrl)]
        x_case = x_case[np.isfinite(x_case)]
        s_ctrl = _safe_stats(x_ctrl)
        s_case = _safe_stats(x_case)

        pc_stats.append(
            {
                "pc": f"PC{pc_num}",
                "n_control": int(x_ctrl.size),
                "ctrl_mean": s_ctrl["mean"],
                "ctrl_std": s_ctrl["std"],
                "ctrl_min": s_ctrl["min"],
                "ctrl_max": s_ctrl["max"],
                "n_case": int(x_case.size),
                "case_mean": s_case["mean"],
                "case_std": s_case["std"],
                "case_min": s_case["min"],
                "case_max": s_case["max"],
                "diff": float(s_case["mean"] - s_ctrl["mean"])
                if np.isfinite(s_case["mean"]) and np.isfinite(s_ctrl["mean"])
                else np.nan,
            }
        )

        t_res = ttest_ind(x_case, x_ctrl, equal_var=False, nan_policy="omit")
        u_res = mannwhitneyu(x_case, x_ctrl, alternative="two-sided")
        test_rows.append(
            {
                "pc": f"PC{pc_num}",
                "t_stat": float(t_res[0]),
                "p_t": float(t_res[1]),
                "u_stat": float(u_res[0]),
                "p_u": float(u_res[1]),
            }
        )

    stats_df = pd.DataFrame(pc_stats).sort_values(
        "pc", key=lambda s: s.str.replace("PC", "", regex=False).astype(int)
    )
    tests_df = pd.DataFrame(test_rows).sort_values(
        "pc", key=lambda s: s.str.replace("PC", "", regex=False).astype(int)
    )
    tests_df["p_t_adj"] = _bh_fdr_adjust(tests_df["p_t"].to_numpy())
    tests_df["p_u_adj"] = _bh_fdr_adjust(tests_df["p_u"].to_numpy())
    tests_df["reject_t"] = tests_df["p_t_adj"] <= 0.05
    tests_df["reject_u"] = tests_df["p_u_adj"] <= 0.05
    return stats_df, tests_df


def _write_pc_distribution_log(
    *,
    log_path: Path,
    plot_df: pd.DataFrame,
    pc_stats: pd.DataFrame,
    tests_df: pd.DataFrame,
    case_label: str,
    ctrl_label: str,
) -> None:
    logger = _setup_file_logger("plot_popgmm_pihat_projection_stats", log_path)

    logger.info("=" * 88)
    logger.info(f"{'RELATEDNESS-AWARE ALL-PC DISTRIBUTION SUMMARY':^88}")
    logger.info("=" * 88)
    logger.info(f"  control_label         : {ctrl_label}")
    logger.info(f"  case_label            : {case_label}")
    logger.info(f"  total_samples         : {len(plot_df):,}")
    logger.info(f"  control_samples       : {int((plot_df['GROUP'] == ctrl_label).sum()):,}")
    logger.info(f"  case_samples          : {int((plot_df['GROUP'] == case_label).sum()):,}")
    logger.info(f"  other_samples         : {int((plot_df['GROUP'] == 'OTHER').sum()):,}")
    logger.info(f"  relatedness_flagged   : {int(plot_df['IS_RELATEDNESS_FLAGGED'].sum()):,}")
    logger.info("-" * 88)
    logger.info("Per-PC summary statistics:")
    logger.info("  PC     n_ctrl   mean_ctrl     std_ctrl     min_ctrl     max_ctrl    n_case   mean_case     std_case     min_case     max_case      diff")
    for _, row in pc_stats.iterrows():
        logger.info(
            f"  {row['pc']:<5}"
            f" {int(row['n_control']):>7,}"
            f" {float(row['ctrl_mean']):>11.4f}"
            f" {float(row['ctrl_std']):>11.4f}"
            f" {float(row['ctrl_min']):>11.4f}"
            f" {float(row['ctrl_max']):>11.4f}"
            f" {int(row['n_case']):>9,}"
            f" {float(row['case_mean']):>11.4f}"
            f" {float(row['case_std']):>11.4f}"
            f" {float(row['case_min']):>11.4f}"
            f" {float(row['case_max']):>11.4f}"
            f" {float(row['diff']):>11.4f}"
        )

    logger.info("-" * 88)
    logger.info("Statistical testing (case vs control):")
    logger.info("  PC      t_stat        p_t      p_t_adj  sig_t     u_stat        p_u      p_u_adj  sig_u")
    for _, row in tests_df.iterrows():
        logger.info(
            f"  {row['pc']:<5}"
            f" {float(row['t_stat']):>9.4f}"
            f" {float(row['p_t']):>11.4e}"
            f" {float(row['p_t_adj']):>11.4e}"
            f" {'Yes' if bool(row['reject_t']) else 'No ':>7}"
            f" {float(row['u_stat']):>11.4f}"
            f" {float(row['p_u']):>11.4e}"
            f" {float(row['p_u_adj']):>11.4e}"
            f" {'Yes' if bool(row['reject_u']) else 'No ':>7}"
        )

    logger.info("-" * 88)
    logger.info(f"  significant_t_tests   : {int(tests_df['reject_t'].sum())}")
    logger.info(f"  significant_u_tests   : {int(tests_df['reject_u'].sum())}")
    logger.info("=" * 88)


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


def _configure_kde_matplotlib() -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    sns.set_context("paper", font_scale=2.0)
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 18,
            "axes.titlesize": 24,
            "axes.labelsize": 20,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "legend.fontsize": 22,
            "figure.titlesize": 30,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.4,
            "xtick.major.width": 1.3,
            "ytick.major.width": 1.3,
            "grid.color": "#C3C3C3",
            "grid.linestyle": "--",
            "grid.alpha": 0.35,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _style_axes(ax: Axes) -> None:
    for spine in ["top", "right"]:
        if spine in ax.spines:
            ax.spines[spine].set_visible(False)
    ax.tick_params(axis="both", which="both", direction="out", top=False, right=False)
    ax.grid(False)


def plot_variance_summary(eigenvals: np.ndarray, out_png: Path) -> tuple[np.ndarray, np.ndarray]:
    pcs = np.arange(1, len(eigenvals) + 1)
    explained = eigenvals / np.sum(eigenvals)
    explained_pct = explained * 100.0
    cumulative_pct = np.cumsum(explained_pct)

    fig, ax1 = plt.subplots(figsize=(8, 4.5))
    xtick_labels = [f"PC{i}" for i in pcs]

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

    handles = [bars, line]
    labels = [h.get_label() for h in handles]
    ax1.legend(handles, labels, loc="center left", bbox_to_anchor=(1.20, 0.5), frameon=False, borderaxespad=0.0)

    fig.subplots_adjust(left=0.10, right=0.78, bottom=0.22, top=0.95)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

    return explained_pct, cumulative_pct


def plot_pairwise_scatter(
    plot_df: pd.DataFrame,
    explained_pct: np.ndarray,
    out_prefix: Path,
    case_label: str,
    ctrl_label: str,
) -> None:
    color_map = {
        ctrl_label: CTRL_COLOR,
        case_label: CASE_COLOR,
        "OTHER": OTHER_COLOR,
    }
    order = [ctrl_label, case_label, "OTHER"]

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

            base_size_map = {ctrl_label: 14.0, case_label: 18.0, "OTHER": 8.0}
            base_alpha_map = {ctrl_label: 0.80, case_label: 0.90, "OTHER": 0.60}
            base_zorder_map = {ctrl_label: 2, case_label: 3, "OTHER": 2}

            group_counts: dict[str, int] = {}
            overlap_counts: dict[str, int] = {}

            for group in order:
                sub = plot_df[(plot_df["GROUP"] == group) & (~plot_df["IS_RELATEDNESS_FLAGGED"])]
                if sub.empty:
                    continue
                group_counts[group] = int((plot_df["GROUP"] == group).sum())
                ax.scatter(
                    sub[x_col],
                    sub[y_col],
                    s=base_size_map.get(group, 8.0),
                    alpha=base_alpha_map.get(group, 0.6),
                    c=color_map[group],
                    edgecolors="none",
                    zorder=base_zorder_map.get(group, 2),
                    rasterized=True,
                )

            # Overlay relatedness-flagged samples (PIHAT vertex-cover derived)
            # with the same group color but a distinct marker and stroke.
            for group in order:
                sub = plot_df[(plot_df["GROUP"] == group) & (plot_df["IS_RELATEDNESS_FLAGGED"])]
                if sub.empty:
                    continue
                overlap_counts[group] = len(sub)
                ax.scatter(
                    sub[x_col],
                    sub[y_col],
                    s=36.0,
                    alpha=1.0,
                    c=color_map[group],
                    marker="^",
                    edgecolors="black",
                    linewidths=0.45,
                    zorder=6,
                    rasterized=True,
                )

            ax.set_xlabel(f"PC{i} ({explained_pct[i-1]:.1f}%)")
            ax.set_ylabel(f"PC{j} ({explained_pct[j-1]:.1f}%)")
            ax.axhline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.5, zorder=0)
            ax.axvline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.5, zorder=0)
            ax.grid(True, linestyle="--", alpha=0.75, linewidth=0.7, color="#C0C0C0", zorder=0)

            legend_handles: list[Line2D] = []
            legend_labels: list[str] = []
            for group in order:
                n_total = int((plot_df["GROUP"] == group).sum())
                if n_total == 0:
                    continue
                n_overlap = overlap_counts.get(group, 0)
                handle = Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    markersize=8.5,
                    markerfacecolor=color_map[group],
                    markeredgecolor="none",
                )
                legend_handles.append(handle)
                legend_labels.append(
                    f"{group} (n={n_total:,}; relatedness-flagged={n_overlap:,})"
                )

            overlap_handle = Line2D(
                [0],
                [0],
                marker="^",
                linestyle="",
                markersize=12.5,
                markerfacecolor=OVERLAY_FACE,
                markeredgecolor="black",
                markeredgewidth=1.1,
            )
            legend_handles.append(overlap_handle)
            legend_labels.append("Relatedness-flagged sample")

            ax.legend(
                handles=legend_handles,
                labels=legend_labels,
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                frameon=False,
                borderaxespad=0.0,
            )

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def _format_pc_label(pc_num: int, explained_pct: np.ndarray) -> str:
    return f"PC{pc_num} ({explained_pct[pc_num - 1]:.1f}%)"


def plot_grouped_pc_kde(
    plot_df: pd.DataFrame,
    explained_pct: np.ndarray,
    out_png: Path,
    case_label: str,
    ctrl_label: str,
    n_cols: int = 5,
) -> None:
    _configure_kde_matplotlib()

    pc_cols = [c for c in plot_df.columns if re.match(r"^PC\d+$", str(c))]
    pc_cols = sorted(pc_cols, key=lambda c: int(str(c).replace("PC", "")))
    if not pc_cols:
        raise ValueError("No PC columns available for KDE plotting.")

    n_cols = max(1, int(n_cols))
    n_rows = int(np.ceil(len(pc_cols) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.2 * n_cols, 4.8 * n_rows), squeeze=False)
    fig.subplots_adjust(wspace=0.30, hspace=0.45, left=0.06, right=0.98, top=0.82, bottom=0.07)

    ctrl_color = CTRL_COLOR
    case_color = CASE_COLOR

    for idx, col in enumerate(pc_cols, start=1):
        ax = axes[(idx - 1) // n_cols][(idx - 1) % n_cols]

        x_ctrl = pd.to_numeric(
            plot_df.loc[plot_df["GROUP"] == ctrl_label, col], errors="coerce"
        ).to_numpy(dtype=float, copy=False)
        x_case = pd.to_numeric(
            plot_df.loc[plot_df["GROUP"] == case_label, col], errors="coerce"
        ).to_numpy(dtype=float, copy=False)
        x_ctrl = x_ctrl[np.isfinite(x_ctrl)]
        x_case = x_case[np.isfinite(x_case)]

        combined = np.concatenate([x_ctrl, x_case]) if (x_ctrl.size or x_case.size) else np.array([], dtype=float)
        if combined.size > 0:
            x_min = float(np.min(combined))
            x_max = float(np.max(combined))
            span = x_max - x_min
            if span == 0:
                span = 1.0
            pad = span * 0.08
            ax.set_xlim(x_min - pad, x_max + pad)

        if x_ctrl.size > 1:
            sns.kdeplot(
                x=x_ctrl,
                color=ctrl_color,
                fill=True,
                alpha=0.60,
                linewidth=2.5,
                ax=ax,
            )
        if x_case.size > 1:
            sns.kdeplot(
                x=x_case,
                color=case_color,
                fill=True,
                alpha=0.60,
                linewidth=2.5,
                ax=ax,
            )

        pc_num = int(str(col).replace("PC", ""))
        ax.set_title(_format_pc_label(pc_num, explained_pct), pad=10, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Density")
        ax.grid(True, linestyle="--", alpha=0.35, color="#C3C3C3")
        ax.tick_params(axis="both", which="major", length=4.8, width=1.1)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.1)

        if ax.get_legend() is not None:
            ax.get_legend().remove()

    for j in range(len(pc_cols) + 1, n_rows * n_cols + 1):
        axes[(j - 1) // n_cols][(j - 1) % n_cols].axis("off")

    fig.suptitle(
        "Case/Control Distributions on Relatedness-Aware PC Axes",
        fontweight="bold",
        y=0.995,
    )
    fig.legend(
        handles=[
            Patch(color=ctrl_color, alpha=0.60, label=ctrl_label),
            Patch(color=case_color, alpha=0.60, label=case_label),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=2,
        frameon=False,
        fontsize=22,
    )
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    if bool(plt.isinteractive()):
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    args = parse_args()
    _configure_matplotlib()

    out_prefix = Path(args.out_prefix)

    eigenvals = _load_eigenval(args.base_eigenval, args.max_pcs)
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

    relatedness_flagged_iids = _load_overlap_iids(args.relatedness_flagged_fid_iid)

    score_df["GROUP"] = score_df["IID"].astype(str).map(lambda x: group_map.get(str(x), "OTHER"))
    score_df["IS_RELATEDNESS_FLAGGED"] = score_df["IID"].astype(str).isin(relatedness_flagged_iids)

    pc_cols = [f"PC{i}" for i in range(1, len(eigenvals) + 1)]
    plot_df = score_df[["IID", "GROUP", "IS_RELATEDNESS_FLAGGED"] + pc_cols].copy()

    variance_png = Path(f"{out_prefix}.variance_summary.png")
    explained_pct, _ = plot_variance_summary(eigenvals, variance_png)

    kde_png = Path(f"{out_prefix}.pc_group_distribution.png")
    log_path = Path(f"{out_prefix}.pc_group_distribution.log.txt")

    pc_stats, tests_df = _build_pc_stats_and_tests(
        plot_df=plot_df,
        pc_cols=pc_cols,
        case_label=args.case_label,
        ctrl_label=args.ctrl_label,
    )

    _write_pc_distribution_log(
        log_path=log_path,
        plot_df=plot_df,
        pc_stats=pc_stats,
        tests_df=tests_df,
        case_label=args.case_label,
        ctrl_label=args.ctrl_label,
    )

    plot_grouped_pc_kde(
        plot_df=plot_df,
        explained_pct=explained_pct,
        out_png=kde_png,
        case_label=args.case_label,
        ctrl_label=args.ctrl_label,
        n_cols=args.kde_n_cols,
    )

    plot_pairwise_scatter(
        plot_df=plot_df,
        explained_pct=explained_pct,
        out_prefix=out_prefix,
        case_label=args.case_label,
        ctrl_label=args.ctrl_label,
    )

    n_relatedness_flagged = int(plot_df["IS_RELATEDNESS_FLAGGED"].sum())
    print(f"[INFO] relatedness-flagged samples: {n_relatedness_flagged:,}")
    print(f"[OK] Variance summary PNG: {variance_png}")
    print(f"[OK] Group KDE PNG: {kde_png}")
    print(f"[OK] Group KDE log: {log_path}")
    print(f"[OK] Pairwise PDF: {out_prefix}.pc_pairs.pdf")


if __name__ == "__main__":
    main()
