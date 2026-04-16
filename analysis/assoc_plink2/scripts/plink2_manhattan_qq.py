#!/usr/bin/env python3
"""Generate Manhattan + QQ plots for PLINK2 GLM outputs.

The script is tuned for PLINK2 association files (e.g. *.glm.logistic.hybrid)
and can automatically split plots by TEST column when present.
"""

from __future__ import annotations

import argparse

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Manhattan and QQ plots for PLINK2 association results")
    parser.add_argument("--input", required=True, help="Path to PLINK2 GLM result file")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
    parser.add_argument("--title", default="PLINK2 Association Results", help="Base title for plot files")
    parser.add_argument("--gwas-threshold", type=float, default=5e-8, help="Genome-wide threshold (default: 5e-8)")
    return parser.parse_args()


def detect_col(cols: list[str], candidates: list[str]) -> str | None:
    lower_to_col = {c.lower(): c for c in cols}
    for candidate in candidates:
        if candidate.lower() in lower_to_col:
            return lower_to_col[candidate.lower()]
    return None


def map_chromosome(chrom: object) -> int:
    value = str(chrom).strip().lower().replace("chr", "")
    if value == "x":
        return 23
    if value == "y":
        return 24
    if value == "xy":
        return 25
    if value in {"m", "mt"}:
        return 26
    return int(value) if value.isdigit() else 99


def calculate_lambda_gc(pvals: np.ndarray) -> float:
    if pvals.size == 0:
        return float("nan")
    median_p = np.nanmedian(pvals)
    obs = stats.chi2.isf(median_p, df=1)
    exp = stats.chi2.ppf(0.5, df=1)
    return float(obs / exp)


def build_plot_frame(df: pd.DataFrame, chr_col: str, pos_col: str, p_col: str) -> pd.DataFrame:
    out = df[[chr_col, pos_col, p_col]].copy()
    out.columns = ["CHR", "POS", "P"]
    out["POS"] = pd.to_numeric(out["POS"], errors="coerce")
    out["P"] = pd.to_numeric(out["P"], errors="coerce")
    out = out.dropna(subset=["POS", "P"])
    out = out[(out["P"] > 0.0) & (out["P"] <= 1.0)]
    out["P"] = out["P"].clip(lower=1e-300, upper=1.0)
    out["CHR_NUM"] = out["CHR"].map(map_chromosome)
    out = out[out["CHR_NUM"] < 99]
    out = out.sort_values(["CHR_NUM", "POS"]).reset_index(drop=True)
    out["LOG10_P"] = -np.log10(out["P"])
    return out


def plot_single(df: pd.DataFrame, title: str, output_file: str, gwas_threshold: float) -> dict[str, object]:
    if df.empty:
        return {
            "output": output_file,
            "variants_used": 0,
            "lambda_gc": np.nan,
            "top_hit_chr": pd.NA,
            "top_hit_pos": pd.NA,
            "top_hit_p": pd.NA,
        }

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans", "sans-serif"],
            "font.size": 18,
            "axes.labelsize": 22,
            "axes.titlesize": 24,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "figure.dpi": 300,
            "axes.linewidth": 2.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    fig = plt.figure(figsize=(24, 9), facecolor="white")
    left, right, bottom, top = 0.08, 0.02, 0.16, 0.90
    gap = 0.06
    plot_h = top - bottom
    qq_w = plot_h * (9 / 24)
    man_w = 1.0 - left - right - gap - qq_w

    rect_man: tuple[float, float, float, float] = (left, bottom, man_w, plot_h)
    rect_qq: tuple[float, float, float, float] = (left + man_w + gap, bottom, qq_w, plot_h)
    ax_man = fig.add_axes(rect_man)
    ax_qq = fig.add_axes(rect_qq)

    n_tests = len(df)
    bonf = -np.log10(0.05 / n_tests) if n_tests > 0 else 8.0
    ymax = max(float(df["LOG10_P"].max()), bonf, 8.0) * 1.25

    chromosomes = sorted(df["CHR_NUM"].unique())
    chr_offsets: dict[int, tuple[float, float]] = {}
    xticks, xlabels = [], []
    current_offset = 0.0

    for chrom in chromosomes:
        c = df[df["CHR_NUM"] == chrom]
        cmin = float(c["POS"].min())
        cmax = float(c["POS"].max())
        span = max(cmax - cmin, 1.0)
        chr_offsets[chrom] = (current_offset, cmin)
        xticks.append(current_offset + span / 2.0)
        xlabels.append({23: "X", 24: "Y", 25: "XY", 26: "MT"}.get(chrom, str(chrom)))
        current_offset += span + 1.0

    colors = ["#1F4E79", "#8DB3E2"]
    keep_thresh = 2.0

    for i, chrom in enumerate(chromosomes):
        c = df[df["CHR_NUM"] == chrom]
        offset, cmin = chr_offsets[chrom]
        c = c.assign(X=offset + (c["POS"] - cmin))

        high = c[c["LOG10_P"] >= keep_thresh]
        low = c[c["LOG10_P"] < keep_thresh]

        if not low.empty:
            low_plot = low.sample(frac=0.1, random_state=42) if len(low) > 10000 else low
            ax_man.scatter(low_plot["X"], low_plot["LOG10_P"], s=18, color=colors[i % 2], alpha=0.9, linewidth=0, rasterized=True)

        if not high.empty:
            ax_man.scatter(high["X"], high["LOG10_P"], s=18, color=colors[i % 2], alpha=0.95, linewidth=0, rasterized=True)

    gws_y = -np.log10(gwas_threshold)
    ax_man.axhline(gws_y, color="#D32F2F", linestyle="--", linewidth=1.8, alpha=0.9, label=r"Genome-wide ($P<5\times10^{-8}$)")
    ax_man.axhline(bonf, color="#7B1FA2", linestyle=":", linewidth=1.8, alpha=0.9, label="Bonferroni")

    hits = df[df["LOG10_P"] >= gws_y]
    if not hits.empty:
        hx, hy = [], []
        for _, row in hits.iterrows():
            off, cmin = chr_offsets[int(row["CHR_NUM"])]
            hx.append(off + (float(row["POS"]) - cmin))
            hy.append(float(row["LOG10_P"]))
        ax_man.scatter(hx, hy, color="#B71C1C", s=24, edgecolor="black", linewidth=0.3, zorder=4)

    ax_man.set_xticks(xticks)
    ax_man.set_xticklabels([lab if i % 2 == 0 else f"\n{lab}" for i, lab in enumerate(xlabels)], fontweight="bold")
    ax_man.tick_params(axis="x", length=0, pad=8)
    ax_man.set_xlim(-current_offset * 0.015, current_offset * 1.015)
    ax_man.set_ylim(0, ymax)
    ax_man.set_xlabel("Chromosome", fontweight="bold")
    ax_man.set_ylabel(r"$-\log_{10}(P)$", fontweight="bold")
    ax_man.set_title(title, fontweight="bold", loc="left")
    ax_man.legend(loc="upper left", frameon=True, edgecolor="#BDBDBD", framealpha=0.9)

    p_sorted = np.sort(df["P"].to_numpy())
    observed = -np.log10(p_sorted)
    expected = -np.log10((np.arange(1, len(p_sorted) + 1) - 0.5) / len(p_sorted))

    try:
        if len(p_sorted) > 10000:
            head = np.arange(1, 1001)
            tail = np.unique(np.geomspace(1001, len(p_sorted), num=5000).astype(int))
            idx = np.concatenate([head, tail])
        else:
            idx = np.arange(1, len(p_sorted) + 1)
        ci_x = -np.log10((idx - 0.5) / len(p_sorted))
        low = stats.beta.ppf(0.025, idx, len(p_sorted) - idx + 1)
        high = stats.beta.ppf(0.975, idx, len(p_sorted) - idx + 1)
        ax_qq.fill_between(ci_x, -np.log10(high), -np.log10(low), color="#B0BEC5", alpha=0.4, zorder=1, label="95% CI")
    except Exception:
        pass

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "plink2_qq",
        list(zip([0.0, 0.3, 0.6, 1.0], ["#B0C4DE", "#4682B4", "#1F4E79", "#000080"])),
    )
    norm = mcolors.Normalize(vmin=0, vmax=max(8.0, float(observed.max())))
    ax_qq.scatter(expected, observed, c=observed, cmap=cmap, norm=norm, s=18, alpha=1.0, linewidth=0, rasterized=True)
    ax_qq.plot([0, ymax], [0, ymax], color="#D50000", linestyle="--", linewidth=1.8)

    lambda_gc = calculate_lambda_gc(p_sorted)
    ax_qq.text(
        0.05,
        0.95,
        f"$\\lambda_{{GC}}={lambda_gc:.4f}$",
        transform=ax_qq.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="left",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.9, "edgecolor": "#BDBDBD"},
    )
    ax_qq.set_aspect("equal")
    ax_qq.set_xlim(0, ymax)
    ax_qq.set_ylim(0, ymax)
    ax_qq.grid(True, linestyle="-", linewidth=0.5, color="#E0E0E0", alpha=1.0)
    ax_qq.set_xlabel(r"Expected $-\log_{10}(P)$", fontweight="bold")
    ax_qq.set_ylabel(r"Observed $-\log_{10}(P)$", fontweight="bold")
    ax_qq.set_title("Q-Q Plot", fontweight="bold", loc="left")

    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    top_idx = df["P"].idxmin()
    top_row = df.loc[top_idx]
    return {
        "output": output_file,
        "variants_used": int(n_tests),
        "lambda_gc": lambda_gc,
        "top_hit_chr": str(top_row["CHR"]),
        "top_hit_pos": int(top_row["POS"]),
        "top_hit_p": float(top_row["P"]),
    }


def main() -> None:
    args = parse_args()

    # Read only key columns first, then subset by TEST to reduce memory pressure.
    header = pd.read_csv(args.input, sep=r"\s+", nrows=0, engine="python")
    cols = list(header.columns)

    chr_col = detect_col(cols, ["#CHROM", "CHROM", "CHR"])
    pos_col = detect_col(cols, ["POS", "BP", "POSITION"])
    p_col = detect_col(cols, ["P", "PVAL", "P_VALUE", "p.value"])
    test_col = detect_col(cols, ["TEST"])

    missing = [name for name, col in [("CHR", chr_col), ("POS", pos_col), ("P", p_col)] if col is None]
    if missing:
        raise ValueError(f"Missing required columns in PLINK2 file: {', '.join(missing)}")

    usecols = [chr_col, pos_col, p_col]
    if test_col:
        usecols.append(test_col)

    df = pd.read_csv(args.input, sep=r"\s+", usecols=usecols, engine="python")

    jobs: list[tuple[str, pd.DataFrame, str]] = []
    if test_col:
        # Exclude non-association rows that can distort QQ/Manhattan interpretation.
        valid_tests = [t for t in sorted(df[test_col].dropna().astype(str).unique()) if t.upper() not in {"INTERCEPT", "FIRTH?"}]
        if not valid_tests:
            valid_tests = ["ALL"]

        for test_name in valid_tests:
            if test_name == "ALL":
                subset = df.copy()
            else:
                subset = df[df[test_col].astype(str) == test_name].copy()
            if subset.empty:
                continue
            safe_test = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in test_name)
            out_png = f"{args.output_prefix}.{safe_test}.manhattan_qq.png"
            plot_title = f"{args.title} [{test_name}]"
            jobs.append((test_name, subset, out_png))
    else:
        jobs.append(("ALL", df, f"{args.output_prefix}.manhattan_qq.png"))

    stats_rows = []
    for test_name, subset, out_png in jobs:
        frame = build_plot_frame(subset, chr_col, pos_col, p_col)
        summary = plot_single(frame, f"{args.title} [{test_name}]", out_png, args.gwas_threshold)
        summary["test"] = test_name
        stats_rows.append(summary)

    stats_df = pd.DataFrame(stats_rows)
    stats_path = f"{args.output_prefix}.plot_stats.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
