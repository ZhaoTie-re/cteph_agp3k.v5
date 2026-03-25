#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Run sample QC filtering from sample QC metrics table.")
	parser.add_argument("--metrics-tsv", required=True, help="Input sample QC metrics TSV")
	parser.add_argument("--config-json", required=True, help="Sample QC config JSON")
	parser.add_argument("--out-prefix", required=True, help="Output file prefix")
	parser.add_argument("--sample-info-xlsx", help="Sample info Excel file for phenotype summary")
	parser.add_argument("--sample-id-col", help="Column name containing sample IDs in sample info")
	parser.add_argument("--phenotype-col", help="Column name containing phenotype labels in sample info")
	parser.add_argument("--case-value", default="PH", help="Case phenotype value")
	parser.add_argument("--ctrl-value", default="AGP3K", help="Control phenotype value")
	parser.add_argument(
		"--case-label",
		help=(
			"Optional display label for case phenotype in plots; "
			"defaults to the value of --case-value if not provided."
		),
	)
	parser.add_argument(
		"--ctrl-label",
		help=(
			"Optional display label for control phenotype in plots; "
			"defaults to the value of --ctrl-value if not provided."
		),
	)
	return parser.parse_args()


def normalize_text(value) -> str | None:
	if value is None or pd.isna(value):
		return None
	text = str(value).strip()
	return text if text else None


def load_config(config_path: Path) -> dict:
	with config_path.open("r", encoding="utf-8") as handle:
		config = json.load(handle)

	condition1 = config.get("condition1", {}) if isinstance(config, dict) else {}
	condition2 = config.get("condition2", {}) if isinstance(config, dict) else {}

	mode = condition1.get("mode", config.get("condition1_mode", "dp_robustz"))
	if mode not in {"dp_robustz", "smiss"}:
		raise ValueError("condition1.mode must be 'dp_robustz' or 'smiss'")

	dp_robustz_threshold = condition1.get("dp_robustz", {}).get(
		"remove_if_less_than",
		config.get("dp_robustz_threshold", -3.0),
	)
	smiss_threshold = condition1.get("smiss", {}).get(
		"remove_if_greater_than",
		config.get("smiss_threshold", 0.05),
	)
	enable_het_f_sd_filter = condition2.get(
		"enabled",
		config.get("enable_het_f_sd_filter", True),
	)
	het_f_sd_multiplier = condition2.get(
		"sd_multiplier",
		config.get("het_f_sd_multiplier", 5.0),
	)

	return {
		"condition1_mode": mode,
		"dp_robustz_threshold": float(dp_robustz_threshold),
		"smiss_threshold": float(smiss_threshold),
		"enable_het_f_sd_filter": bool(enable_het_f_sd_filter),
		"het_f_sd_multiplier": float(het_f_sd_multiplier),
	}


def read_metrics(metrics_path: Path) -> pd.DataFrame:
	df = pd.read_csv(metrics_path, sep="\t", dtype=str)
	required_cols = ["#FID", "IID", "Het_F", "Target_DP", "DP", "SMISS", "DP_RobustZ_in_TargetDP"]
	missing_cols = [col for col in required_cols if col not in df.columns]
	if missing_cols:
		raise ValueError(f"Missing required columns in metrics table: {', '.join(missing_cols)}")

	df["Het_F"] = pd.to_numeric(df["Het_F"], errors="coerce")
	df["DP"] = pd.to_numeric(df["DP"], errors="coerce")
	df["SMISS"] = pd.to_numeric(df["SMISS"], errors="coerce")
	df["DP_RobustZ_in_TargetDP"] = pd.to_numeric(df["DP_RobustZ_in_TargetDP"], errors="coerce")

	return df


def attach_phenotype(
	df: pd.DataFrame,
	xlsx_path: Path,
	sample_id_col: str,
	phenotype_col: str,
	case_value: str,
	ctrl_value: str,
) -> pd.DataFrame:
	"""Attach phenotype labels (case/ctrl) to metrics table based on sample info.

	If any of the parameters are missing, the input dataframe is returned with
	PHENO_LABEL/PHENO_GROUP columns set to "unknown".
	"""
	result = df.copy()
	result["PHENO_LABEL"] = None
	result["PHENO_GROUP"] = "unknown"

	if not xlsx_path or not xlsx_path.exists() or not sample_id_col or not phenotype_col:
		return result

	info_df = pd.read_excel(xlsx_path, engine="openpyxl")
	required = [sample_id_col, phenotype_col]
	missing = [col for col in required if col not in info_df.columns]
	if missing:
		return result

	sub = info_df[[sample_id_col, phenotype_col]].copy()
	sub.columns = ["IID", "PHENO_RAW"]
	sub["IID"] = sub["IID"].map(normalize_text)
	sub["PHENO_RAW"] = sub["PHENO_RAW"].map(normalize_text)
	sub = sub[sub["IID"].notna()].copy()

	# Drop conflicting duplicate phenotype annotations
	duplicated = sub["IID"].duplicated(keep=False)
	if duplicated.any():
		dup = sub.loc[duplicated].drop_duplicates(subset=["IID", "PHENO_RAW"], keep=False)
		conflict_ids = dup["IID"].unique().tolist()
		if conflict_ids:
			preview = ", ".join(sorted(conflict_ids)[:5])
			print(
				f"WARNING: Conflicting phenotype labels for IIDs: {preview}; phenotype summary will ignore these samples.",
				file=sys.stderr,
			)
		sub = sub[~sub["IID"].isin(conflict_ids)]

	phenotype_map = sub.set_index("IID")["PHENO_RAW"]
	result["PHENO_LABEL"] = result["IID"].map(phenotype_map)

	case_norm = normalize_text(case_value)
	ctrl_norm = normalize_text(ctrl_value)

	def _group(label: str | None) -> str:
		if label is None or (isinstance(label, float) and pd.isna(label)):
			return "unknown"
		if label == case_norm:
			return "case"
		if label == ctrl_norm:
			return "ctrl"
		return "other"

	result["PHENO_GROUP"] = result["PHENO_LABEL"].map(_group)
	return result


def apply_filters(df: pd.DataFrame, config: dict) -> tuple[pd.DataFrame, dict]:
	mode = config["condition1_mode"]
	if mode == "dp_robustz":
		fail_cond1 = df["DP_RobustZ_in_TargetDP"] < config["dp_robustz_threshold"]
		cond1_desc = f"DP_RobustZ_in_TargetDP < {config['dp_robustz_threshold']}"
		cond1_col = "DP_RobustZ_in_TargetDP"
	else:
		fail_cond1 = df["SMISS"] > config["smiss_threshold"]
		cond1_desc = f"SMISS > {config['smiss_threshold']}"
		cond1_col = "SMISS"

	fail_cond1 = fail_cond1.fillna(False)

	if config["enable_het_f_sd_filter"]:
		het_mean = df["Het_F"].mean(skipna=True)
		het_sd = df["Het_F"].std(skipna=True, ddof=1)
		if pd.isna(het_mean) or pd.isna(het_sd) or het_sd == 0:
			fail_het = pd.Series([False] * len(df), index=df.index)
			het_low = float("nan")
			het_high = float("nan")
		else:
			k = config["het_f_sd_multiplier"]
			het_low = het_mean - k * het_sd
			het_high = het_mean + k * het_sd
			fail_het = (df["Het_F"] < het_low) | (df["Het_F"] > het_high)
			fail_het = fail_het.fillna(False)
	else:
		fail_het = pd.Series([False] * len(df), index=df.index)
		het_mean = float("nan")
		het_sd = float("nan")
		het_low = float("nan")
		het_high = float("nan")

	df_out = df.copy()
	df_out["QC_FAIL_COND1"] = fail_cond1
	df_out["QC_FAIL_HET_F_SD"] = fail_het
	df_out["QC_REMOVE"] = df_out["QC_FAIL_COND1"] | df_out["QC_FAIL_HET_F_SD"]

	def reason_row(row: pd.Series) -> str:
		reasons = []
		if row["QC_FAIL_COND1"]:
			reasons.append("COND1")
		if row["QC_FAIL_HET_F_SD"]:
			reasons.append("HET_F_SD")
		return ";".join(reasons)

	df_out["QC_REASON"] = df_out.apply(reason_row, axis=1)

	stats = {
		"n_total": int(len(df_out)),
		"n_remove": int(df_out["QC_REMOVE"].sum()),
		"n_keep": int((~df_out["QC_REMOVE"]).sum()),
		"n_fail_cond1": int(df_out["QC_FAIL_COND1"].sum()),
		"n_fail_het_f_sd": int(df_out["QC_FAIL_HET_F_SD"].sum()),
		"condition1_mode": mode,
		"condition1_description": cond1_desc,
		"condition1_col": cond1_col,
		"dp_robustz_threshold": config["dp_robustz_threshold"],
		"smiss_threshold": config["smiss_threshold"],
		"het_filter_enabled": bool(config["enable_het_f_sd_filter"]),
		"het_f_sd_multiplier": float(config["het_f_sd_multiplier"]),
		"het_f_mean": None if pd.isna(het_mean) else float(het_mean),
		"het_f_sd": None if pd.isna(het_sd) else float(het_sd),
		"het_f_low": None if pd.isna(het_low) else float(het_low),
		"het_f_high": None if pd.isna(het_high) else float(het_high),
	}

	return df_out, stats


def write_outputs(df: pd.DataFrame, stats: dict, out_prefix: Path) -> None:
	"""Write QC detail table, ID lists, and summary files."""
	detail_path = Path(f"{out_prefix}.sample_qc.detail.tsv")
	remove_path = Path(f"{out_prefix}.sample_qc.remove.id")
	keep_path = Path(f"{out_prefix}.sample_qc.keep.id")
	summary_json_path = Path(f"{out_prefix}.sample_qc.summary.json")
	summary_txt_path = Path(f"{out_prefix}.sample_qc.summary.txt")

	df.to_csv(detail_path, sep="\t", index=False)

	remove_df = df.loc[df["QC_REMOVE"], ["#FID", "IID"]].copy()
	keep_df = df.loc[~df["QC_REMOVE"], ["#FID", "IID"]].copy()
	remove_df.to_csv(remove_path, sep="\t", header=False, index=False)
	keep_df.to_csv(keep_path, sep="\t", header=False, index=False)

	with summary_json_path.open("w", encoding="utf-8") as handle:
		json.dump(stats, handle, ensure_ascii=False, indent=2)

	with summary_txt_path.open("w", encoding="utf-8") as handle:
		handle.write(f"n_total\t{stats['n_total']}\n")
		handle.write(f"n_remove\t{stats['n_remove']}\n")
		handle.write(f"n_keep\t{stats['n_keep']}\n")
		handle.write(f"n_fail_cond1\t{stats['n_fail_cond1']}\n")
		handle.write(f"n_fail_het_f_sd\t{stats['n_fail_het_f_sd']}\n")
		handle.write(f"condition1_mode\t{stats['condition1_mode']}\n")
		handle.write(f"condition1_description\t{stats['condition1_description']}\n")


def generate_visualization(
	df: pd.DataFrame,
	stats: dict,
	case_value: str,
	ctrl_value: str,
	out_prefix: Path,
	case_label: str | None = None,
	ctrl_label: str | None = None,
) -> None:
	"""Generate publication-quality QC visualization with joint scatter + KDE.

	Design目标：
	1. 散点图：x 轴为 mean DP，y 轴为 Het_F；
	   用两条虚线表示 Het_F outlier 阈值；
	   通过点形状标记 condition1（DP-robustZ / SMISS）outlier 样本。
	2. 针对 15x、30x 的 mean DP 分布（KDE），与散点图在 x 轴上对齐绘制。
	3. 右侧使用表格汇总 15x/30x + QC category（pass / F outlier /
	   condition1 outlier / both）的 case / ctrl 计数。
	"""

	# Make fonts reasonably large for publication-quality figures
	plt.rcParams.update(
		{
			"font.size": 13,
			"axes.titlesize": 18,
			"axes.labelsize": 15,
			"legend.fontsize": 12,
			"xtick.labelsize": 12,
			"ytick.labelsize": 12,
		}
	)

	cond1_mode = stats["condition1_mode"]
	het_f_low = stats["het_f_low"]
	het_f_high = stats["het_f_high"]
	het_k = stats["het_f_sd_multiplier"]
	# Display labels fall back to phenotype values when not explicitly provided
	case_label = str(case_label if case_label is not None else case_value)
	ctrl_label = str(ctrl_label if ctrl_label is not None else ctrl_value)
	cond1_name = "DP" if cond1_mode == "dp_robustz" else "SMISS"

	# Build mathematical QC title shared by the entire figure
	if cond1_mode == "dp_robustz":
		thr = stats.get("dp_robustz_threshold", -3.0)
		# Z^{robust}_DP: depth-normalised robust Z-score of coverage
		cond1_math = f"$Z_{{DP}}^{{\\mathrm{{robust}}}} < {thr:.1f}$"
	else:
		thr = stats.get("smiss_threshold", 0.05)
		# s_miss: sample-level missingness rate
		cond1_math = f"$s_{{\\mathrm{{miss}}}} > {thr:.3f}$"

	if stats.get("het_filter_enabled", False) and not pd.isna(het_f_low):
		# F_het: autosomal heterozygosity F; |F_het - mu_F| > k * sigma_F,
		# equivalent to F_het not in [mu_F - k * sigma_F, mu_F + k * sigma_F]
		het_math = f"$|F_{{\\mathrm{{het}}}} - \\mu_F| > {het_k:g}\\,\\sigma_F$"
		interval = f" ($F_{{\\mathrm{{het}}}} \\notin [{het_f_low:.3f}, {het_f_high:.3f}]$)"
		qc_title = f"QC thresholds: {cond1_math}; {het_math}{interval}"
	else:
		qc_title = f"QC thresholds: {cond1_math}; Het_F filter not applied"

	df_plot = df.copy()

	# Ensure phenotype groups are available (fall back to FID if needed)
	if "PHENO_GROUP" not in df_plot.columns:
		def _pheno_from_fid(fid: str | float) -> str:
			if isinstance(fid, float) and pd.isna(fid):
				return "unknown"
			if str(fid) == str(case_value):
				return "case"
			if str(fid) == str(ctrl_value):
				return "ctrl"
			return "other"

		df_plot["PHENO_GROUP"] = df_plot["#FID"].map(_pheno_from_fid)

	# Target_DP label: 原始列如果是 "15x" / "30x" 则直接保留
	def _depth_label(val) -> str:
		if val is None or (isinstance(val, float) and pd.isna(val)):
			return "Other"
		text = str(val).strip().lower()
		if text in {"15x", "15"}:
			return "15x"
		if text in {"30x", "30"}:
			return "30x"
		try:
			iv = int(round(float(text)))
			if iv == 15:
				return "15x"
			if iv == 30:
				return "30x"
		except Exception:  # noqa: BLE001
			pass
		return "Other"

	df_plot["TARGET_DP_LABEL"] = df_plot["Target_DP"].map(_depth_label)

	# QC base category flags（与原逻辑一致，供汇总表使用）
	def _qc_base(row: pd.Series) -> str:
		if not row["QC_REMOVE"]:
			return "Pass"
		fail_dp = bool(row["QC_FAIL_COND1"])
		fail_f = bool(row["QC_FAIL_HET_F_SD"])
		if fail_dp and fail_f:
			return "Both Outlier"
		if fail_dp:
			return "Only Cond1 Outlier"
		if fail_f:
			return "Only F Outlier"
		return "Pass"

	df_plot["QC_CATEGORY_BASE"] = df_plot.apply(_qc_base, axis=1)

	# 颜色方案：
	#  - 15x 与 30x 通过颜色深浅区分；
	#  - 每个 QC category（Pass / F / DP / Both）在 15x 和 30x 上使用同一色系的浅色 / 深色；
	#  - KDE：黑色曲线，15x 用虚线，30x 用实线。
	depth_linestyle = {"15x": "--", "30x": "-"}
	category_depth_colors: dict[tuple[str, str], str] = {
		("Pass", "15x"): "#c6dbef",   # light blue
		("Pass", "30x"): "#6baed6",   # darker blue
		("Only F Outlier", "15x"): "#d9f0d3",  # light green
		("Only F Outlier", "30x"): "#74c476",  # darker green
		("Only Cond1 Outlier", "15x"): "#fee6ce",  # light orange
		("Only Cond1 Outlier", "30x"): "#fdae6b",  # darker orange
		("Both Outlier", "15x"): "#fcbba1",  # light red
		("Both Outlier", "30x"): "#fb6a4a",  # darker red
	}
	depth_palette = {"15x": category_depth_colors[("Pass", "15x")], "30x": category_depth_colors[("Pass", "30x")], "Other": "0.6"}
	plotted_depths: list[str] = []

	# ---- Figure layout: 上方一行是 15x vs 30x KDE，下方是散点图，右侧是表格 ----
	# 更偏向左侧主图，同时整体尺寸和字体放大以满足发表级别需求
	fig = plt.figure(figsize=(17, 10.5))
	gs = fig.add_gridspec(
		2,
		2,
		width_ratios=[3.6, 2.0],
		height_ratios=[1.2, 2.8],
		hspace=0.1,
		wspace=0.35,
	)

	ax_scatter = fig.add_subplot(gs[1, 0])
	ax_dp = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
	ax_table = fig.add_subplot(gs[:, 1])
	ax_table.axis("off")

	# ---- 1. KDE of mean DP：只区分 15x / 30x，case-ctrl 合并 ----
	def _plot_dp_kde_combined(ax):
		plotted = False
		for depth in ["15x", "30x"]:
			sub = df_plot[df_plot["TARGET_DP_LABEL"] == depth]
			vals = sub["DP"].dropna()
			if len(vals) <= 1:
				continue
			sns.kdeplot(
				x=vals,
				ax=ax,
				label=depth,
				color="black",
				linestyle=depth_linestyle.get(depth, "-"),
				linewidth=2.0,
			)
			plotted = True
			if depth not in plotted_depths:
				plotted_depths.append(depth)

		if not plotted:
			ax.set_visible(False)
			return

		ax.set_ylabel("KDE")
		ax.grid(alpha=0.25)

	_plot_dp_kde_combined(ax_dp)
	# Figure-level title placed above the legend row
	fig.suptitle(qc_title, y=0.99, fontweight="bold")

	# 只在最下方散点图上显示 x 轴刻度；KDE 面板不再重复 x 轴标题
	plt.setp(ax_dp.get_xticklabels(), visible=False)
	ax_dp.set_xlabel("")

	# ---- 2. 散点图：x = Mean DP, y = Het_F
	# 形状：区分 case / ctrl；颜色：区分 QC category，与右侧表格底纹一致 ----
	scatter_df = df_plot.dropna(subset=["DP", "Het_F"]).copy()
	cat_colors = {
		# 用于 legend 的代表颜色：使用 30x 的较深色
		"Pass": category_depth_colors[("Pass", "30x")],
		"Only F Outlier": category_depth_colors[("Only F Outlier", "30x")],
		"Only Cond1 Outlier": category_depth_colors[("Only Cond1 Outlier", "30x")],
		"Both Outlier": category_depth_colors[("Both Outlier", "30x")],
	}
	pheno_markers = {
		"case": "o",      # circle
		"ctrl": "s",      # square
		"other": "D",     # diamond
		"unknown": "x",   # x marker
	}

	for pheno in ["case", "ctrl", "other", "unknown"]:
		for depth in ["15x", "30x", "Other"]:
			for base_key in ["Pass", "Only F Outlier", "Only Cond1 Outlier", "Both Outlier"]:
				mask = (
					(scatter_df["PHENO_GROUP"] == pheno)
					& (scatter_df["TARGET_DP_LABEL"] == depth)
					& (scatter_df["QC_CATEGORY_BASE"] == base_key)
				)
				sub = scatter_df[mask]
				if sub.empty:
					continue
				if depth in {"15x", "30x"}:
					color_val = category_depth_colors.get((base_key, depth), "0.6")
				else:
					color_val = "0.6"
				alpha_val = 0.9 if base_key == "Pass" else 1.0
				if depth == "30x":
					z_val = 3
				elif depth == "15x":
					z_val = 2
				else:
					z_val = 1
				ax_scatter.scatter(
					sub["DP"],
					sub["Het_F"],
					c=color_val,
					marker=pheno_markers.get(pheno, "o"),
					s=32,
					alpha=alpha_val,
					edgecolors="black",
					linewidths=0.3,
					zorder=z_val,
				)

	# 构造两个图例：一个解释形状（Phenotype），一个解释颜色（QC category）
	cond1_label = "Only DP outlier" if cond1_mode == "dp_robustz" else "Only SMISS outlier"
	legend_names = {
		"Only F Outlier": "Only F outlier",
		"Only Cond1 Outlier": cond1_label,
		"Both Outlier": "Both outlier",
	}

	shape_handles = [
		Line2D(
			[],
			[],
			marker=pheno_markers[ph],
			color="black",
			linestyle="None",
			markersize=7,
			label=(case_label if ph == "case" else ctrl_label if ph == "ctrl" else ph),
		)
		for ph in ["case", "ctrl"]
	]
	color_handles = [
		Line2D(
			[],
			[],
			marker="o",
			color="black",
			linestyle="None",
			markerfacecolor=cat_colors[key],
			markersize=8,
			label=legend_names[key],
		)
		for key in ["Only F Outlier", "Only Cond1 Outlier", "Both Outlier"]
	]

	dp_handles = [
		Line2D(
			[],
			[],
			color="black",
			linestyle=depth_linestyle.get(d, "-"),
			linewidth=2.0,
			label=f"{d}",
		)
		for d in plotted_depths
	]

	# Het_F outlier 阈值（两条虚线，水平线）
	if not pd.isna(het_f_low) and not pd.isna(het_f_high):
		ax_scatter.axhline(het_f_low, color="darkred", linestyle="--", linewidth=1.4)
		ax_scatter.axhline(het_f_high, color="darkred", linestyle="--", linewidth=1.4)

	ax_scatter.set_xlabel("Mean DP")
	ax_scatter.set_ylabel("Heterozygosity F")
	ax_scatter.grid(alpha=0.3)

	# 统一图例：在整张图上方分三组展示 Target DP / Phenotype / QC category
	if dp_handles:
		fig.legend(
			handles=dp_handles,
			title="Target DP",
			loc="upper center",
			bbox_to_anchor=(0.23, 0.955),
			ncol=len(dp_handles),
			fontsize=12,
			framealpha=0.95,
		)
	if shape_handles:
		fig.legend(
			handles=shape_handles,
			title="Phenotype",
			loc="upper center",
			bbox_to_anchor=(0.52, 0.955),
			ncol=len(shape_handles),
			fontsize=12,
			framealpha=0.95,
		)
	if color_handles:
		fig.legend(
			handles=color_handles,
			title="QC category",
			loc="upper center",
			bbox_to_anchor=(0.84, 0.955),
			ncol=len(color_handles),
			fontsize=12,
			framealpha=0.95,
		)

	# ---- 3. Summary table: counts by depth and QC category ----
	depths = ["15x", "30x"]
	base_keys = ["Pass", "Only F Outlier", "Only Cond1 Outlier", "Both Outlier"]
	display_names = {
		"Pass": "Pass",
		"Only F Outlier": "Only F outlier",
		"Only Cond1 Outlier": cond1_label,
		"Both Outlier": "Both outlier",
	}

	headers = ["Category", "Case", "Ctrl", "Total"]
	table_rows: list[list[str]] = []
	row_styles: list[str] = []  # one of "pass", "f", "cond1", "both"

	def _fmt(n: int) -> str:
		"""Format integer with comma as thousands separator."""
		return f"{n:,}"

	row_colors: list[str] = []
	for base_key in base_keys:
		for depth in depths:
			mask = (df_plot["TARGET_DP_LABEL"] == depth) & (df_plot["QC_CATEGORY_BASE"] == base_key)
			sub = df_plot[mask]
			case_count = int((sub["PHENO_GROUP"] == "case").sum()) if not sub.empty else 0
			ctrl_count = int((sub["PHENO_GROUP"] == "ctrl").sum()) if not sub.empty else 0
			total_count = int(len(sub))
			row_label = f"{display_names[base_key]} ({depth})"
			table_rows.append(
				[row_label, _fmt(case_count), _fmt(ctrl_count), _fmt(total_count)],
			)
			if (base_key, depth) in category_depth_colors:
				row_colors.append(category_depth_colors[(base_key, depth)])
			else:
				row_colors.append("#ffffff")

	cell_text = [headers] + table_rows
	table = ax_table.table(
		cellText=cell_text,
		cellLoc="center",
		loc="center",
		colWidths=[0.6, 0.13, 0.13, 0.14],
	)
	table.auto_set_font_size(False)
	table.set_fontsize(13)
	table.scale(1.1, 1.8)

	for i in range(len(cell_text)):
		for j in range(len(headers)):
			cell = table[(i, j)]
			if i == 0:
				cell.set_facecolor("#40466e")
				cell.set_text_props(weight="bold", color="white")
			else:
				cell.set_facecolor(row_colors[i - 1])

	plot_path = Path(f"{out_prefix}.sample_qc.png")
	plt.savefig(plot_path, dpi=300, bbox_inches="tight")
	plt.close()


def main() -> int:
	args = parse_args()
	try:
		metrics = read_metrics(Path(args.metrics_tsv))
		if args.sample_info_xlsx and args.sample_id_col and args.phenotype_col:
			metrics = attach_phenotype(
				metrics,
				Path(args.sample_info_xlsx),
				args.sample_id_col,
				args.phenotype_col,
				args.case_value,
				args.ctrl_value,
			)
		config = load_config(Path(args.config_json))
		filtered, stats = apply_filters(metrics, config)
		write_outputs(filtered, stats, Path(args.out_prefix))
		case_label = args.case_label if hasattr(args, "case_label") and args.case_label else args.case_value
		ctrl_label = args.ctrl_label if hasattr(args, "ctrl_label") and args.ctrl_label else args.ctrl_value
		generate_visualization(
			filtered,
			stats,
			args.case_value,
			args.ctrl_value,
			Path(args.out_prefix),
			case_label,
			ctrl_label,
		)

		print(
			(
				f"Sample QC finished: total={stats['n_total']}, "
				f"remove={stats['n_remove']}, keep={stats['n_keep']}"
			),
			file=sys.stderr,
		)
		return 0
	except Exception as exc:  # noqa: BLE001
		print(f"ERROR: {exc}", file=sys.stderr)
		return 1


if __name__ == "__main__":
	sys.exit(main())