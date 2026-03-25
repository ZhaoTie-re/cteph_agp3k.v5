#!/usr/bin/env python3

"""PI_HAT-based relatedness QC using weighted vertex cover.

This script:

1. Reads a PLINK .genome file with pairwise PI_HAT statistics.
2. Filters pairs with PI_HAT > threshold.
3. Builds an undirected graph where nodes are samples (IID) and edges are
   high-PI_HAT pairs.
4. Assigns node weights using sample-level missingness (SMISS) and
   case/control status:

   - Case samples: very large weights (hard to remove), scaled by 1/SMISS.
   - Control samples: smaller weights, also scaled by 1/SMISS.

5. Runs a minimum weighted vertex cover approximation to find a minimal-
   weight set of samples whose removal would break all high-PI_HAT edges.
6. Writes:

   - "*.pi_hat.pairs.tsv": all high-PI_HAT pairs with PI_HAT and pair type.
   - "*.pi_hat.vertex_cover_samples.tsv": per-sample table including
     status, SMISS, degree, weight, and whether selected by vertex cover.
   - "*.pi_hat.log.txt": human-readable description of the algorithm and
     summary statistics of the relatedness pruning strategy.

Note: this script DOES NOT modify any genotype files; it only produces
annotations and documentation of which samples would be removed by the
strategy.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Circle
from networkx.algorithms.approximation import min_weighted_vertex_cover


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run PI_HAT-based relatedness QC (weighted vertex cover)."
    )
    parser.add_argument("--genome", required=True, help="PLINK .genome file path")
    parser.add_argument("--metrics-tsv", required=True, help="Sample QC metrics TSV (with SMISS)")
    parser.add_argument(
        "--sample-info-xlsx",
        required=True,
        help="Sample info Excel, used for case/control labels",
    )
    parser.add_argument(
        "--sample-id-col",
        required=True,
        help="Column in sample info corresponding to IID (e.g. 'ID JHRPv5')",
    )
    parser.add_argument(
        "--phenotype-col",
        required=True,
        help="Column in sample info containing phenotype labels (e.g. 'Outcome')",
    )
    parser.add_argument(
        "--case-value",
        default="PH",
        help="Phenotype value representing cases (default: PH)",
    )
    parser.add_argument(
        "--ctrl-value",
        default="AGP3K",
        help="Phenotype value representing controls (default: AGP3K)",
    )
    parser.add_argument(
        "--pi-hat-threshold",
        type=float,
        default=0.20,
        help="PI_HAT threshold for defining related pairs (default: 0.20)",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix for result tables and log",
    )
    return parser.parse_args()


def normalize_text(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    return text


def load_genome(genome_path: Path, threshold: float) -> pd.DataFrame:
    """Load PLINK .genome and filter by PI_HAT > threshold."""
    df = pd.read_csv(
        genome_path,
        sep=r"\s+",
        dtype={0: str, 1: str},
        low_memory=False,
    )

    required = ["IID1", "IID2", "PI_HAT"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s) in genome file: {', '.join(missing)}")

    df["PI_HAT"] = pd.to_numeric(df["PI_HAT"], errors="coerce")
    df = df[df["PI_HAT"] > threshold].copy()
    df = df.dropna(subset=["IID1", "IID2", "PI_HAT"])
    return df.reset_index(drop=True)


def load_metrics(metrics_path: Path) -> pd.DataFrame:
    df = pd.read_csv(metrics_path, sep="\t", dtype=str)
    if "IID" not in df.columns or "SMISS" not in df.columns:
        raise ValueError("metrics TSV must contain 'IID' and 'SMISS' columns")
    df["IID"] = df["IID"].map(normalize_text)
    df["SMISS"] = pd.to_numeric(df["SMISS"], errors="coerce")
    df = df[df["IID"].notna()].copy()
    return df


def load_phenotype(
    xlsx_path: Path,
    sample_id_col: str,
    phenotype_col: str,
    case_value: str,
    ctrl_value: str,
) -> pd.DataFrame:
    df = pd.read_excel(xlsx_path, engine="openpyxl")
    missing = [c for c in [sample_id_col, phenotype_col] if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required column(s) in sample info: {', '.join(missing)}"
        )

    sub = df[[sample_id_col, phenotype_col]].copy()
    sub.columns = ["IID", "PHENO_RAW"]
    sub["IID"] = sub["IID"].map(normalize_text)
    sub["PHENO_RAW"] = sub["PHENO_RAW"].map(normalize_text)
    sub = sub[sub["IID"].notna()].copy()

    # Drop conflicting duplicate phenotype annotations
    duplicated = sub["IID"].duplicated(keep=False)
    if duplicated.any():
        dup_rows = sub.loc[duplicated].drop_duplicates(
            subset=["IID", "PHENO_RAW"], keep=False
        )
        conflict_ids = dup_rows["IID"].unique().tolist()
        if conflict_ids:
            preview = ", ".join(sorted(conflict_ids)[:5])
            print(
                f"WARNING: Conflicting phenotype labels for IIDs: {preview}; "
                "phenotype labels for these samples will be treated as 'other'.",
                file=sys.stderr,
            )
        sub = sub[~sub["IID"].isin(conflict_ids)]

    case_norm = normalize_text(case_value)
    ctrl_norm = normalize_text(ctrl_value)

    def _group(label: str | None) -> str:
        if label is None or (isinstance(label, float) and pd.isna(label)):
            return "other"
        if label == case_norm:
            return "case"
        if label == ctrl_norm:
            return "ctrl"
        return "other"

    sub["PHENO_GROUP"] = sub["PHENO_RAW"].map(_group)
    return sub[["IID", "PHENO_GROUP"]].copy()


def build_graph(
    pairs: pd.DataFrame,
    metrics_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
) -> tuple[nx.Graph, pd.DataFrame]:
    """Build weighted graph from high-PI_HAT pairs and annotate node info.

    Node weight rule:
      - case  : weight = 1e6 / SMISS
      - ctrl  : weight = 1.0 / SMISS
      - other : treat as control (1.0 / SMISS)
    """

    node_ids = sorted(set(pairs["IID1"]).union(set(pairs["IID2"])))

    # Merge metrics (SMISS) and phenotype labels
    node_df = pd.DataFrame({"IID": node_ids})
    node_df = node_df.merge(
        metrics_df[["IID", "SMISS"]], on="IID", how="left", validate="one_to_one"
    )
    node_df = node_df.merge(pheno_df, on="IID", how="left")

    # Default SMISS for missing values to a small positive number to avoid div-by-zero
    node_df["SMISS"] = node_df["SMISS"].fillna(0.001)

    def _status(group: str | None) -> str:
        if group == "case":
            return "case"
        if group == "ctrl":
            return "ctrl"
        return "other"

    node_df["STATUS"] = node_df["PHENO_GROUP"].map(_status)

    def _weight(row: pd.Series) -> float:
        smiss = float(row["SMISS"]) if not pd.isna(row["SMISS"]) else 0.001
        smiss = max(smiss, 1e-4)  # safety floor
        if row["STATUS"] == "case":
            return 1e6 / smiss
        else:  # ctrl or other
            return 1.0 / smiss

    node_df["WEIGHT"] = node_df.apply(_weight, axis=1)

    G = nx.Graph()
    for _, r in node_df.iterrows():
        G.add_node(
            r["IID"],
            status=r["STATUS"],
            smiss=float(r["SMISS"]),
            weight=float(r["WEIGHT"]),
        )

    # Add edges with PI_HAT as an edge attribute so that
    # visualization can scale line width and add labels.
    for _, row in pairs.iterrows():
        a = row["IID1"]
        b = row["IID2"]
        try:
            pi_hat_val = float(row["PI_HAT"])
        except Exception:  # noqa: BLE001
            continue
        u, v = (a, b) if a <= b else (b, a)
        if G.has_edge(u, v):
            existing = G[u][v].get("pi_hat", pi_hat_val)
            if pi_hat_val > existing:
                G[u][v]["pi_hat"] = pi_hat_val
        else:
            G.add_edge(u, v, pi_hat=pi_hat_val)

    # degree
    node_df["DEGREE"] = [G.degree(iid) for iid in node_df["IID"]]

    return G, node_df


def classify_pair_type(row: pd.Series, status_map: dict[str, str]) -> str:
    s1 = status_map.get(row["IID1"], "other")
    s2 = status_map.get(row["IID2"], "other")
    if s1 == "case" and s2 == "case":
        return "case-case"
    if s1 == "ctrl" and s2 == "ctrl":
        return "ctrl-ctrl"
    if (s1 == "case" and s2 == "ctrl") or (s1 == "ctrl" and s2 == "case"):
        return "case-ctrl"
    return "other"


def run_vertex_cover(G: nx.Graph) -> set[str]:
    """Compute (near-)exact minimum-weight vertex cover.

    For each connected component we run a small branch-and-bound search on a
    bitmask representation of the component to obtain the true minimum-weight
    vertex cover. This is exact for components up to a moderate size. For very
    large components (more than 25 nodes), we fall back to the NetworkX
    approximation for safety.
    """

    if G.number_of_edges() == 0:
        return set()

    selected: set[str] = set()

    for comp_nodes in nx.connected_components(G):
        H = G.subgraph(comp_nodes).copy()
        if H.number_of_edges() == 0:
            continue

        nodes = list(H.nodes())
        n = len(nodes)

        # For unusually large components, fall back to approximation to avoid
        # exponential blow-up in worst-case graphs.
        if n > 25:
            approx_vc = min_weighted_vertex_cover(H, weight="weight")
            selected.update(approx_vc)
            continue

        index = {node: i for i, node in enumerate(nodes)}
        weights = [float(H.nodes[node].get("weight", 1.0)) for node in nodes]
        edges = [(index[u], index[v]) for u, v in H.edges()]

        best_weight = float("inf")
        best_mask = 0

        def backtrack(mask: int, current_weight: float) -> None:
            nonlocal best_weight, best_mask

            # Prune if already worse than the best known solution
            if current_weight >= best_weight:
                return

            # Find the first uncovered edge
            uncovered_edge = None
            for i, j in edges:
                if not ((mask >> i) & 1 or (mask >> j) & 1):
                    uncovered_edge = (i, j)
                    break

            # All edges are covered: update best solution
            if uncovered_edge is None:
                best_weight = current_weight
                best_mask = mask
                return

            u, v = uncovered_edge

            # Branch 1: include u in the cover
            if not ((mask >> u) & 1):
                backtrack(mask | (1 << u), current_weight + weights[u])

            # Branch 2: include v in the cover
            if not ((mask >> v) & 1):
                backtrack(mask | (1 << v), current_weight + weights[v])

        backtrack(0, 0.0)

        # Decode best_mask to node IDs
        for i, node in enumerate(nodes):
            if (best_mask >> i) & 1:
                selected.add(node)

    return selected


def write_outputs(
    pairs: pd.DataFrame,
    node_df: pd.DataFrame,
    to_remove: set[str],
    threshold: float,
    out_prefix: Path,
) -> None:
    # Annotate pair type for pairs table
    status_map = dict(zip(node_df["IID"], node_df["STATUS"]))
    pairs = pairs.copy()
    pairs["PAIR_TYPE"] = pairs.apply(
        lambda r: classify_pair_type(r, status_map), axis=1
    )

    # Use explicit suffix appending so that if out_prefix already
    # ends with ".pi_hat" we get files like "*.pi_hat.pairs.tsv",
    # which matches the Nextflow patterns.
    pairs_out = Path(f"{out_prefix}.pairs.tsv")
    pairs.to_csv(pairs_out, sep="\t", index=False)

    node_df = node_df.copy()
    node_df["SELECTED_FOR_REMOVAL"] = node_df["IID"].isin(to_remove)
    node_out = Path(f"{out_prefix}.vertex_cover_samples.tsv")
    node_df[[
        "IID",
        "STATUS",
        "SMISS",
        "DEGREE",
        "WEIGHT",
        "SELECTED_FOR_REMOVAL",
    ]].sort_values(["SELECTED_FOR_REMOVAL", "STATUS", "DEGREE"], ascending=[False, True, False]).to_csv(
        node_out, sep="\t", index=False
    )

    # Summary log
    n_pairs = len(pairs)
    n_nodes = len(node_df)
    n_to_remove = int(node_df["SELECTED_FOR_REMOVAL"].sum())
    n_case = int((node_df["STATUS"] == "case").sum())
    n_ctrl = int((node_df["STATUS"] == "ctrl").sum())
    n_other = n_nodes - n_case - n_ctrl

    log_out = Path(f"{out_prefix}.log.txt")
    with log_out.open("w", encoding="utf-8") as handle:
        handle.write("PI_HAT-based relatedness QC (weighted vertex cover)\n")
        handle.write("=" * 70 + "\n\n")

        handle.write("1. Input and thresholds\n")
        handle.write(f"  - PI_HAT source  : PLINK .genome\n")
        handle.write(f"  - PI_HAT cutoff  : PI_HAT > {threshold:.3f}\n")
        handle.write(
            "  - Graph model    : nodes = samples (IID); edges = pairs with PI_HAT above cutoff\n"
        )
        handle.write(
            "  - Node weights   :\n"
            "      * case  : weight = 1e6 / SMISS\n"
            "      * ctrl/other : weight = 1.0 / SMISS\n"
        )
        handle.write(
            "    Samples with higher missingness (SMISS) thus receive smaller weights\n"
        )
        handle.write(
            "    and are more likely to be selected for removal, especially controls.\n\n"
        )

        handle.write("2. Graph statistics (after PI_HAT filtering)\n")
        handle.write(f"  - Number of high-PI_HAT pairs  : {n_pairs}\n")
        handle.write(f"  - Number of unique samples     : {n_nodes}\n")
        handle.write(f"      * case   : {n_case}\n")
        handle.write(f"      * ctrl   : {n_ctrl}\n")
        handle.write(f"      * other  : {n_other}\n\n")

        handle.write("3. Minimum weighted vertex cover\n")
        handle.write(
            "  - Objective: find a minimum-total-weight set of vertices whose\n"
        )
        handle.write(
            "    removal covers all high-PI_HAT edges (no remaining pair has\n"
        )
        handle.write(
            "    PI_HAT above the threshold).\n"
        )
        handle.write(
            "  - Algorithm: the graph is decomposed into connected components;\n"
        )
        handle.write(
            "    components with up to 25 vertices are solved exactly via a\n"
        )
        handle.write(
            "    branch-and-bound search on a bitmask representation, while\n"
        )
        handle.write(
            "    larger components use NetworkX's min_weighted_vertex_cover\n"
        )
        handle.write(
            "    approximation to avoid exponential worst-case time.\n"
        )
        handle.write(
            "  - Interpretation: samples with SELECTED_FOR_REMOVAL=TRUE in\n"
        )
        handle.write(
            "    vertex_cover_samples.tsv form a minimum-weight set under the\n"
        )
        handle.write(
            "    specified case/control and SMISS-based weights, and are\n"
        )
        handle.write(
            "    highlighted in the network figure with a red outline and\n"
        )
        handle.write(
            "    '[REMOVED]' label.\n"
        )
        handle.write(
            "\n"
        )
        handle.write(
            "4. Important note\n"
        )
        handle.write(
            "  This step does NOT directly remove samples from the genotype.\n"
        )
        handle.write(
            "  Instead, it provides annotations (vertex_cover_samples.tsv) that\n"
        )
        handle.write(
            "  document which samples would be removed under this strategy.\n"
        )
        handle.write(
            "  Downstream analyses may choose to exclude these samples based on\n"
        )
        handle.write(
            "  study-specific requirements.\n"
        )
        handle.write("\n5. Summary of removal counts\n")
        handle.write(f"  - Selected samples for removal : {n_to_remove}\n")

        if n_to_remove > 0:
            sub = node_df[node_df["SELECTED_FOR_REMOVAL"]].copy()
            n_case_del = int((sub["STATUS"] == "case").sum())
            n_ctrl_del = int((sub["STATUS"] == "ctrl").sum())
            n_other_del = n_to_remove - n_case_del - n_ctrl_del
            handle.write(f"      * case   : {n_case_del}\n")
            handle.write(f"      * ctrl   : {n_ctrl_del}\n")
            handle.write(f"      * other  : {n_other_del}\n")


def plot_network_figure(
    G: nx.Graph,
    node_df: pd.DataFrame,
    to_remove: set[str],
    threshold: float,
    out_prefix: Path,
) -> None:
    """Visualize the kinship graph by connected components in a grid.

    Each connected component (cluster of related samples) is drawn in its own
    subplot, which avoids the overplotting that occurs when all components
    are embedded in a single panel. Nodes are colored by case/control/other
    status, and nodes selected by the weighted vertex cover are highlighted
    with a red outline.

    The figure is saved as a PNG at "{out_prefix}.network.png" and is intended
    for academic/paper-quality inspection of the relatedness pruning.
    """
    # Global style for publication-quality visualization
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    if G.number_of_nodes() == 0:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(
            0.5,
            0.5,
            "No PI_HAT pairs above threshold",
            ha="center",
            va="center",
            fontsize=10,
        )
        ax.axis("off")
        out_png = Path(f"{out_prefix}.network.png")
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    # Colorblind-friendly, publication-style palette
    status_color = {"case": "#D55E00", "ctrl": "#0072B2", "other": "#999999"}
    status_map = dict(zip(node_df["IID"], node_df["STATUS"]))
    smiss_map = dict(zip(node_df["IID"], node_df["SMISS"]))
    to_remove = set(to_remove)

    # Connected components sorted by size (largest first). Show ALL components.
    components = [list(c) for c in nx.connected_components(G)]
    components.sort(key=len, reverse=True)
    n_comp = len(components)

    # Determine grid layout (up to 3 columns)
    n_cols = min(3, n_comp)
    n_rows = (n_comp + n_cols - 1) // n_cols

    # Use relatively large per-panel size for clarity
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.0 * n_cols, 4.5 * n_rows))
    if n_rows * n_cols == 1:
        axes = [[axes]]
    elif n_rows == 1:
        axes = [axes]

    for idx, comp_nodes in enumerate(components):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row][col]

        H = G.subgraph(comp_nodes).copy()

        # Layout per component to avoid overlap between clusters
        try:
            from networkx.drawing.nx_agraph import graphviz_layout  # type: ignore

            pos = graphviz_layout(H, prog="neato")
        except Exception:  # noqa: BLE001
            pos = nx.spring_layout(H, seed=42)

        node_colors = [status_color.get(status_map.get(iid, "other"), "#999999") for iid in H.nodes]

        # Scale edge width by PI_HAT value and annotate PI_HAT on edges.
        edges = list(H.edges())
        if edges:
            edge_pihats = [
                float(H[u][v].get("pi_hat", threshold)) for u, v in edges
            ]
            max_pihat = max(edge_pihats)

            edge_widths: list[float] = []
            for val in edge_pihats:
                if max_pihat <= threshold:
                    w = 0.6
                else:
                    norm = (val - threshold) / (max_pihat - threshold)
                    norm = max(0.0, min(1.0, norm))
                    w = 0.6 + 2.9 * norm  # roughly 0.6–3.5
                edge_widths.append(w)

            nx.draw_networkx_edges(
                H,
                pos,
                ax=ax,
                width=edge_widths,
                alpha=0.7,
                edge_color="#B0B0B0",
            )

            # Add PI_HAT labels at edge midpoints (mathematical-style text)
            for (u, v), val in zip(edges, edge_pihats):
                if u not in pos or v not in pos:
                    continue
                x1, y1 = pos[u]
                x2, y2 = pos[v]
                xm, ym = (x1 + x2) / 2.0, (y1 + y2) / 2.0
                # Use a unicode pi with combining hat to avoid LaTeX parsing
                # issues while still looking mathematical.
                label = f"π̂={val:.2f}"
                ax.text(
                    xm,
                    ym,
                    label,
                    fontsize=5,
                    ha="center",
                    va="center",
                    color="#555555",
                    bbox={
                        "boxstyle": "round,pad=0.1",
                        "fc": "white",
                        "ec": "none",
                        "alpha": 0.7,
                    },
                    clip_on=True,
                )

        # Node sizes: emphasise removed samples.
        node_sizes = [120.0 if iid in to_remove else 70.0 for iid in H.nodes]
        nx.draw_networkx_nodes(
            H,
            pos,
            ax=ax,
            node_color=node_colors,
            node_size=node_sizes,
            linewidths=0.0,
        )

        # Highlight nodes selected for removal in this component
        for iid in H.nodes:
            if iid in to_remove and iid in pos:
                x, y = pos[iid]
                ax.add_patch(
                    Circle(
                        (x, y),
                        radius=0.08,
                        fill=False,
                        edgecolor="red",
                        linewidth=1.5,
                    )
                )

        # Annotate each node with IID and SMISS
        for iid in H.nodes:
            if iid not in pos:
                continue
            x, y = pos[iid]
            smiss_val = smiss_map.get(iid)
            if pd.isna(smiss_val):
                label = f"{iid}\nSMISS=NA"
            else:
                label = f"{iid}\nSMISS={float(smiss_val):.3f}"

            # Explicitly mark removed samples in the label.
            if iid in to_remove:
                label = f"{label}\n[REMOVED]"
                text_color = "red"
            else:
                text_color = "black"
            ax.text(
                x,
                y,
                label,
                fontsize=6,
                ha="center",
                va="center",
                color=text_color,
                clip_on=True,
            )

        removed_here = sum(1 for iid in H.nodes if iid in to_remove)
        ax.set_title(
            f"Component {idx + 1} (n={H.number_of_nodes()}, removed={removed_here})",
            fontsize=9,
        )
        ax.axis("off")

    # Turn off any unused subplots
    for k in range(n_comp, n_rows * n_cols):
        row = k // n_cols
        col = k % n_cols
        axes[row][col].axis("off")

    legend_elements = [
        Patch(facecolor=status_color["case"], edgecolor="none", label="Case"),
        Patch(facecolor=status_color["ctrl"], edgecolor="none", label="Control"),
        Patch(facecolor=status_color["other"], edgecolor="none", label="Other"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=3,
        frameon=False,
        fontsize=9,
    )

    fig.suptitle(
        f"Kinship components by π̂ (π̂ > {threshold:.2f})",
        fontsize=12,
        fontweight="bold",
        y=1.04,
    )
    fig.tight_layout(rect=[0, 0.0, 1, 0.98])

    out_png = Path(f"{out_prefix}.network.png")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()

    genome_path = Path(args.genome)
    metrics_path = Path(args.metrics_tsv)
    xlsx_path = Path(args.sample_info_xlsx)
    out_prefix = Path(args.out_prefix)

    try:
        pairs = load_genome(genome_path, args.pi_hat_threshold)
        metrics_df = load_metrics(metrics_path)
        pheno_df = load_phenotype(
            xlsx_path,
            args.sample_id_col,
            args.phenotype_col,
            args.case_value,
            args.ctrl_value,
        )

        G, node_df = build_graph(pairs, metrics_df, pheno_df)
        to_remove = run_vertex_cover(G)
        write_outputs(pairs, node_df, to_remove, args.pi_hat_threshold, out_prefix)
        # Academic-style visualization of the kinship network and vertex cover
        plot_network_figure(G, node_df, to_remove, args.pi_hat_threshold, out_prefix)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
