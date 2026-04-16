#!/usr/bin/env python3
"""Generate a production-ready GWAS model summary report.

This script integrates additive/dominant/recessive PLINK2 association outputs,
annotates significant loci with cohort and ToMMo VCF evidence, computes
case/control genotype metrics (counts, AAF, HWE), and exports a stable,
schema-controlled CSV for downstream statistical and reporting workflows.
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

LOGGER = logging.getLogger("gwas_model_summary")
DEFAULT_P_THRESHOLD = 5e-8
DEFAULT_PLINK2_PATH = "/home/b/b37974/plink2_alpha6/plink2"
OUTPUT_COLUMNS = [
    "#CHROM", "POS", "ID", "rsID", "Gene", "REF", "ALT",
    "TEST", "OR", "CI95", "P",
    "Multi-Allelic TAG", "Allelic NUM", "Allelic Records",
    "Case HomRef Count", "Case Het Count", "Case HomAlt Count", "Case Miss Count", "Case Total Count",
    "Ctrl HomRef Count", "Ctrl Het Count", "Ctrl HomAlt Count", "Ctrl Miss Count", "Ctrl Total Count",
    "Case HomRef Freq", "Case Het Freq", "Case HomAlt Freq", "Case Miss Freq",
    "Ctrl HomRef Freq", "Ctrl Het Freq", "Ctrl HomAlt Freq", "Ctrl Miss Freq",
    "Case-Ctrl Total Miss Freq",
    "Case AAF", "Ctrl AAF", "ToMMo AAF", "ToMMo Filter",
    "Case HWE P-value", "Ctrl HWE P-value",
    "ANN Summary", "ANN Records",
]


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments.

    Hyphenated options are primary. Legacy underscore options are retained as
    aliases for backward compatibility.
    """
    parser = argparse.ArgumentParser(description="Integrate GWAS model results and variant annotations.")
    parser.add_argument("--bed-prefix", "--bed_prefix", dest="bed_prefix", required=True,
                        help="PLINK bed/bim/fam prefix.")
    parser.add_argument("--tommo-vcf-file", "--tommo_vcf_file", dest="tommo_vcf_file", required=True,
                        help="Path to ToMMo VCF/VCF.GZ file.")
    parser.add_argument("--cohort-vcf-path", "--cohort_vcf_path", "--jhrp4_vcf_path", dest="cohort_vcf_path", required=True,
                        help="Path to cohort VCF source (directory of split VCFs or single VCF file).")
    parser.add_argument("--p-threshold", "--p_threshold", dest="p_threshold", type=float, default=DEFAULT_P_THRESHOLD,
                        help=f"P-value threshold for significant loci (default: {DEFAULT_P_THRESHOLD}).")
    parser.add_argument("--add-path", "--add_path", dest="add_path", required=True,
                        help="Path to additive model PLINK2 result.")
    parser.add_argument("--dom-path", "--dom_path", dest="dom_path", required=True,
                        help="Path to dominant model PLINK2 result.")
    parser.add_argument("--rec-path", "--rec_path", dest="rec_path", required=True,
                        help="Path to recessive model PLINK2 result.")
    parser.add_argument("--plink2-path", "--plink2_path", dest="plink2_path", default=DEFAULT_PLINK2_PATH,
                        help=f"PLINK2 executable path (default: {DEFAULT_PLINK2_PATH}).")
    parser.add_argument("--plink2-threads", "--plink2_threads", dest="plink2_threads", type=int, default=16,
                        help="Thread count used for all PLINK2 sub-commands (default: 16).")
    parser.add_argument("--output", default="gwas_summary.plink2.csv",
                        help="Output CSV file path (default: gwas_summary.plink2.csv).")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO",
                        help="Logging level (default: INFO).")
    return parser.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s [%(levelname)s] [%(name)s] [pid=%(process)d] %(message)s",
    )


def validate_runtime(args: argparse.Namespace) -> None:
    """Validate runtime prerequisites before heavy computation starts."""
    if args.plink2_threads < 1:
        raise ValueError(f"plink2-threads must be >= 1, got: {args.plink2_threads}")

    required_files = [
        f"{args.bed_prefix}.bed",
        f"{args.bed_prefix}.bim",
        f"{args.bed_prefix}.fam",
        args.add_path,
        args.dom_path,
        args.rec_path,
        args.tommo_vcf_file,
    ]
    missing = [path for path in required_files if not Path(path).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required input files: {', '.join(missing)}")

    cohort_path = Path(args.cohort_vcf_path)
    if not (cohort_path.is_file() or cohort_path.is_dir()):
        raise FileNotFoundError(
            "cohort-vcf-path must be a VCF file or a directory of split VCFs: "
            f"{args.cohort_vcf_path}"
        )

    plink2_path = Path(args.plink2_path)
    if plink2_path.exists():
        if not os.access(str(plink2_path), os.X_OK):
            raise PermissionError(f"PLINK2 exists but is not executable: {args.plink2_path}")
    elif shutil.which(args.plink2_path) is None:
        raise FileNotFoundError(f"PLINK2 executable not found: {args.plink2_path}")

    if shutil.which("bcftools") is None:
        raise FileNotFoundError("Required executable not found in PATH: bcftools")


def ensure_output_parent(output_path: str) -> None:
    """Create output directory when writing to nested locations."""
    parent = Path(output_path).resolve().parent
    parent.mkdir(parents=True, exist_ok=True)


def enforce_columns(df: pd.DataFrame, expected: list[str]) -> pd.DataFrame:
    """Guarantee stable schema and fail early if required columns are missing."""
    missing = [column for column in expected if column not in df.columns]
    if missing:
        raise KeyError(f"Final summary is missing expected columns: {missing}")
    return df[expected]


def normalize_contig(contig: str) -> str:
    """Normalize contig labels for robust comparison."""
    normalized = str(contig).strip()
    if normalized.lower().startswith("chr"):
        normalized = normalized[3:]
    if normalized.upper() in {"M", "MT"}:
        return "MT"
    return normalized


def contig_query_candidates(contig: str) -> list[str]:
    """Generate query aliases for region-based VCF lookup."""
    raw = str(contig).strip()
    if not raw:
        return []

    base = normalize_contig(raw)
    candidates: list[str] = []

    for candidate in (raw, base):
        if candidate and candidate not in candidates:
            candidates.append(candidate)

    chr_form = base if base.lower().startswith("chr") else f"chr{base}"
    if chr_form not in candidates:
        candidates.append(chr_form)

    if base.upper() == "MT":
        for candidate in ("MT", "M", "chrMT", "chrM"):
            if candidate not in candidates:
                candidates.append(candidate)

    return candidates


def fetch_vcf_records(vcf_file: str, chrom: str, pos: str | int) -> list[str]:
    """Fetch records at one locus from a single VCF file using alias queries."""
    for contig in contig_query_candidates(chrom):
        region = f"{contig}:{pos}-{pos}"
        try:
            result = subprocess.run(
                ["bcftools", "view", "-r", region, vcf_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError:
            continue

        records = [line for line in result.stdout.splitlines() if line and not line.startswith("#")]
        if records:
            return records

    return []


def candidate_split_vcf_files(cohort_vcf_path: str, chrom: str) -> list[str]:
    """Return split VCF files that likely contain the requested contig.

    Matching is prefix-agnostic and based on filename tokens + contig aliases.
    """
    path = Path(cohort_vcf_path)
    if path.is_file():
        return [str(path)]
    if not path.is_dir():
        return []

    contig_aliases = {normalize_contig(alias) for alias in contig_query_candidates(chrom)}
    matched_files: list[str] = []

    for entry in sorted(path.iterdir()):
        name = entry.name
        if not (name.endswith(".vcf.gz") or name.endswith(".vcf")):
            continue

        stem = name
        for suffix in (".vcf.gz", ".vcf"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break

        tokens = [normalize_contig(t) for t in stem.replace("-", ".").replace("_", ".").split(".")]
        if any(token in contig_aliases for token in tokens):
            matched_files.append(str(entry))

    return matched_files


def fetch_cohort_records(cohort_vcf_path: str, chrom: str, pos: str | int) -> list[str]:
    """Fetch records at one locus from cohort VCF source (single or split)."""
    source = Path(cohort_vcf_path)
    if source.is_file():
        return fetch_vcf_records(str(source), chrom, pos)

    combined: list[str] = []
    for vcf_file in candidate_split_vcf_files(cohort_vcf_path, chrom):
        records = fetch_vcf_records(vcf_file, chrom, pos)
        if records:
            combined.extend(records)
    return combined


def read_significant_snps(plink2_result_path: str, p_threshold: float) -> pd.DataFrame:
    """Load significant SNPs from one PLINK2 GLM result file."""
    desired_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "TEST", "OR", "CI95", "P"]
    significant_rows: list[list[object]] = []

    with open(plink2_result_path, "r", encoding="utf-8") as handle:
        header = handle.readline().strip().split()
        idx = {name: header.index(name) for name in ["#CHROM", "POS", "ID", "REF", "ALT", "TEST", "OR", "L95", "U95", "P"]}

        for line in handle:
            fields = line.strip().split()
            try:
                if float(fields[idx["P"]]) >= p_threshold:
                    continue
            except ValueError:
                continue

            significant_rows.append([
                fields[idx["#CHROM"]],
                fields[idx["POS"]],
                fields[idx["ID"]],
                fields[idx["REF"]],
                fields[idx["ALT"]],
                fields[idx["TEST"]],
                fields[idx["OR"]],
                (fields[idx["L95"]], fields[idx["U95"]]),
                fields[idx["P"]],
            ])

    if not significant_rows:
        return pd.DataFrame(columns=desired_cols)
    return pd.DataFrame(significant_rows, columns=desired_cols)


def merge_gwas_models(df_add: pd.DataFrame, df_dom: pd.DataFrame, df_rec: pd.DataFrame) -> pd.DataFrame:
    """Merge ADD/DOM/REC results by locus."""
    model_frames = []
    for frame, label in ((df_add, "ADD"), (df_dom, "DOM"), (df_rec, "REC")):
        if frame.empty:
            continue
        updated = frame.copy()
        updated["MODEL"] = label
        model_frames.append(updated)

    if not model_frames:
        return pd.DataFrame(columns=["#CHROM", "POS", "ID", "REF", "ALT", "TEST", "OR", "CI95", "P"])

    merged_rows = []
    all_df = pd.concat(model_frames, ignore_index=True)
    for key, group in all_df.groupby(["#CHROM", "POS", "ID", "REF", "ALT"]):
        by_model = {row["MODEL"]: row for _, row in group.iterrows()}
        tests, ors, cis, pvals = [], [], [], []
        for model in ("ADD", "DOM", "REC"):
            if model not in by_model:
                continue
            row = by_model[model]
            tests.append(row["TEST"])
            ors.append(row["OR"])
            cis.append(row["CI95"])
            pvals.append(row["P"])
        merged_rows.append(list(key) + [tests, ors, cis, pvals])

    merged_df = pd.DataFrame(merged_rows, columns=["#CHROM", "POS", "ID", "REF", "ALT", "TEST", "OR", "CI95", "P"])

    def chrom_sort_key(chrom: str) -> tuple[object, ...]:
        value = normalize_contig(chrom)
        try:
            return (0, int(value))
        except ValueError:
            return (1, value)

    merged_df = merged_df.sort_values(
        by=["#CHROM", "POS"],
        key=lambda col: col.map(chrom_sort_key) if col.name == "#CHROM" else pd.to_numeric(col),
    ).reset_index(drop=True)
    return merged_df


def check_multiallelic_sites(focus_loci: list[str], cohort_vcf_path: str) -> pd.DataFrame:
    """Check multi-allelic status of all loci in cohort VCF source."""
    tags, counts, records_map = {}, {}, {}

    for locus in focus_loci:
        chrom, pos, ref, _ = locus.split(":")
        try:
            alt_values: list[str] = []
            locus_records: list[str] = []
            for line in fetch_cohort_records(cohort_vcf_path, chrom, pos):
                fields = line.split("\t")
                if normalize_contig(fields[0]) != normalize_contig(chrom) or fields[1] != pos or fields[3] != ref:
                    continue
                alts = fields[4].split(",")
                alt_values.extend(alts)
                locus_records.extend([f"{fields[0]}:{fields[1]}:{fields[3]}:{alt}" for alt in alts])

            if not alt_values:
                tags[locus] = None
                counts[locus] = None
                records_map[locus] = None
                continue

            unique_alts = list(dict.fromkeys(alt_values))
            tags[locus] = len(unique_alts) > 1
            counts[locus] = len(unique_alts) + 1
            records_map[locus] = list(dict.fromkeys(locus_records))
        except Exception as exc:  # pylint: disable=broad-except
            LOGGER.warning("Failed multi-allelic check for %s: %s", locus, exc)
            tags[locus] = None
            counts[locus] = None
            records_map[locus] = None

    return pd.DataFrame(
        {
            "ID": focus_loci,
            "Multi-Allelic TAG": [tags.get(x) for x in focus_loci],
            "Allelic NUM": [counts.get(x) for x in focus_loci],
            "Allelic Records": [records_map.get(x) for x in focus_loci],
        }
    )


def run_plink_subset(
    plink2_path: str,
    plink2_threads: int,
    bed_prefix: str,
    sample_df: pd.DataFrame,
    out_suffix: str,
    plink_args: list[str],
    focus_loci: list[str] | None,
) -> pd.DataFrame:
    """Run PLINK2 on a sample subset and return the output table."""
    keep_file = tempfile.NamedTemporaryFile(mode="w+", delete=False)
    extract_file = None
    out_file = tempfile.NamedTemporaryFile(suffix=out_suffix, delete=False)

    try:
        sample_df.to_csv(keep_file.name, sep="\t", index=False, header=False)
        keep_file.close()

        if focus_loci:
            extract_file = tempfile.NamedTemporaryFile(mode="w+", delete=False)
            extract_file.write("\n".join(focus_loci) + "\n")
            extract_file.close()

        out_prefix = out_file.name.rsplit(".", 1)[0]
        out_file.close()

        cmd = [
            plink2_path,
            "--bfile", bed_prefix,
            "--keep", keep_file.name,
            "--threads", str(plink2_threads),
        ] + plink_args + ["--out", out_prefix]
        if extract_file is not None:
            cmd += ["--extract", extract_file.name]

        LOGGER.debug("Running PLINK2: %s", " ".join(cmd))
        subprocess.run(cmd, check=True)

        result_path = f"{out_prefix}{out_suffix}"
        return pd.read_csv(result_path, sep=r"\s+")
    finally:
        for temp_path in (keep_file.name, getattr(extract_file, "name", None), out_file.name):
            if temp_path and os.path.exists(temp_path):
                try:
                    os.remove(temp_path)
                except OSError:
                    pass
        result_file = f"{out_file.name.rsplit('.', 1)[0]}{out_suffix}"
        if os.path.exists(result_file):
            try:
                os.remove(result_file)
            except OSError:
                pass


def get_case_ctrl_iids(bed_prefix: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    fam = pd.read_csv(f"{bed_prefix}.fam", sep=r"\s+", header=None)
    fam.columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
    return fam[fam["PHENO"] == 2][["FID", "IID"]], fam[fam["PHENO"] == 1][["FID", "IID"]]


def count_genotypes_by_group(
    plink2_path: str,
    plink2_threads: int,
    bed_prefix: str,
    focus_loci: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute genotype count and frequency metrics for case/control."""
    case_iids, ctrl_iids = get_case_ctrl_iids(bed_prefix)
    df_case = run_plink_subset(plink2_path, plink2_threads, bed_prefix, case_iids, ".gcount", ["--geno-counts"], focus_loci)
    df_ctrl = run_plink_subset(plink2_path, plink2_threads, bed_prefix, ctrl_iids, ".gcount", ["--geno-counts"], focus_loci)

    for frame in (df_case, df_ctrl):
        frame["TOTAL_CT"] = frame[["HOM_REF_CT", "HET_REF_ALT_CTS", "TWO_ALT_GENO_CTS", "MISSING_CT"]].sum(axis=1)
        frame["HOM_REF_FREQ"] = np.round(np.where(frame["TOTAL_CT"] != 0, frame["HOM_REF_CT"] / frame["TOTAL_CT"], np.nan), 4)
        frame["HET_REF_ALT_FREQ"] = np.round(np.where(frame["TOTAL_CT"] != 0, frame["HET_REF_ALT_CTS"] / frame["TOTAL_CT"], np.nan), 4)
        frame["TWO_ALT_GENO_FREQ"] = np.round(np.where(frame["TOTAL_CT"] != 0, frame["TWO_ALT_GENO_CTS"] / frame["TOTAL_CT"], np.nan), 4)
        frame["MISSING_FREQ"] = np.round(np.where(frame["TOTAL_CT"] != 0, frame["MISSING_CT"] / frame["TOTAL_CT"], np.nan), 4)

    return df_case, df_ctrl


def compute_hwe_by_group(
    plink2_path: str,
    plink2_threads: int,
    bed_prefix: str,
    focus_loci: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute HWE p-values for case/control."""
    case_iids, ctrl_iids = get_case_ctrl_iids(bed_prefix)
    case_hwe = run_plink_subset(plink2_path, plink2_threads, bed_prefix, case_iids, ".hardy", ["--hardy"], focus_loci)[["ID", "P"]]
    ctrl_hwe = run_plink_subset(plink2_path, plink2_threads, bed_prefix, ctrl_iids, ".hardy", ["--hardy"], focus_loci)[["ID", "P"]]
    return case_hwe, ctrl_hwe


def compute_aaf_by_group(
    plink2_path: str,
    plink2_threads: int,
    bed_prefix: str,
    focus_loci: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute ALT allele frequency (AAF) for case/control."""
    case_iids, ctrl_iids = get_case_ctrl_iids(bed_prefix)
    case_aaf = run_plink_subset(plink2_path, plink2_threads, bed_prefix, case_iids, ".afreq", ["--freq"], focus_loci)[["ID", "ALT_FREQS"]]
    ctrl_aaf = run_plink_subset(plink2_path, plink2_threads, bed_prefix, ctrl_iids, ".afreq", ["--freq"], focus_loci)[["ID", "ALT_FREQS"]]
    return case_aaf, ctrl_aaf


def query_variant_from_tommo(tommo_vcf_file: str, variant_key: str) -> tuple[str | None, set[str], str | None, list[str] | None]:
    """Query rsID/gene/filter/AF from ToMMo VCF for one variant key."""
    chrom, pos, ref, alt = variant_key.strip().split(":")
    records = fetch_vcf_records(tommo_vcf_file, chrom, pos)
    if not records:
        return None, set(), None, None

    rsid = None
    genes: set[str] = set()
    filter_status = None
    allele_freq = None

    for line in records:
        fields = line.split("\t")
        if normalize_contig(fields[0]) != normalize_contig(chrom) or fields[1] != pos or fields[3] != ref:
            continue
        alts = fields[4].split(",")
        if alt not in alts:
            continue

        rsid = fields[2]
        filter_status = fields[6]
        info = fields[7]

        for item in info.split(";"):
            if item.startswith("AF="):
                af = item.split("=", 1)[1]
                allele_freq = af.split(",") if "," in af else [af]
            if item.startswith("ANN="):
                for ann in item[4:].split(","):
                    parts = ann.split("|")
                    if len(parts) > 3 and parts[3]:
                        genes.add(parts[3])
        break

    return rsid, genes, filter_status, allele_freq


def batch_query_variants(tommo_vcf_file: str, variant_list: list[str]) -> pd.DataFrame:
    """Batch query ToMMo annotations for all loci."""
    rows = []
    for variant in variant_list:
        rsid, genes, filter_status, af = query_variant_from_tommo(tommo_vcf_file, variant)
        rows.append(
            {
                "ID": variant,
                "rsID": rsid,
                "Gene": genes if genes else pd.NA,
                "ToMMo Filter": filter_status,
                "ToMMo AAF": ",".join(af) if af else None,
            }
        )
    return pd.DataFrame(rows)


def parse_ann_entries(info_field: str, alt_filter: str | None = None) -> list[dict[str, str]]:
    """Parse ANN entries from one INFO field."""
    ann_blob = None
    for item in info_field.split(";"):
        if item.startswith("ANN="):
            ann_blob = item[4:]
            break
    if not ann_blob:
        return []

    ann_fields = [
        "Allele", "Consequence", "Impact", "Gene_symbol", "Gene_ID",
        "Feature_type", "Transcript_ID", "Biotype", "Rank",
        "HGVS_c", "HGVS_p", "cDNA_pos", "CDS_pos", "AA_pos",
        "Distance", "Errors",
    ]

    annotations: list[dict[str, str]] = []
    for entry in ann_blob.split(","):
        values = entry.split("|")
        parsed = {ann_fields[i]: values[i] if i < len(values) else "" for i in range(len(ann_fields))}
        if alt_filter and parsed["Allele"] != alt_filter:
            continue
        annotations.append(parsed)
    return annotations


def extract_ann_records_from_variant(tommo_vcf_file: str, variant_key: str) -> list[dict[str, str]]:
    """Extract ANN records for one variant from ToMMo VCF."""
    chrom, pos, ref, alt = variant_key.strip().split(":")
    records = fetch_vcf_records(tommo_vcf_file, chrom, pos)
    for line in records:
        fields = line.split("\t")
        if normalize_contig(fields[0]) != normalize_contig(chrom) or fields[1] != str(pos) or fields[3] != ref:
            continue
        if alt not in fields[4].split(","):
            continue
        return parse_ann_entries(fields[7], alt_filter=alt)
    return []


def rename_columns(df: pd.DataFrame, group: str, col_map: dict[str, str]) -> pd.DataFrame:
    return df.rename(columns={k: f"{group} {v}" for k, v in col_map.items()})


def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)
    validate_runtime(args)
    ensure_output_parent(args.output)

    LOGGER.info("Loading significant loci from model results")
    df_add = read_significant_snps(args.add_path, args.p_threshold)
    df_dom = read_significant_snps(args.dom_path, args.p_threshold)
    df_rec = read_significant_snps(args.rec_path, args.p_threshold)

    merged_df = merge_gwas_models(df_add, df_dom, df_rec)
    focus_loci = merged_df["ID"].tolist()
    LOGGER.info("Significant merged loci: %d", len(focus_loci))

    if not focus_loci:
        LOGGER.warning("No loci passed the p-threshold (%.3e); writing an empty summary table.", args.p_threshold)
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(args.output, index=False)
        LOGGER.info("Summary written: %s", args.output)
        return

    LOGGER.info("Annotating multi-allelic status from cohort VCF source")
    multiallelic_df = check_multiallelic_sites(focus_loci, args.cohort_vcf_path)

    LOGGER.info("Computing genotype metrics with PLINK2")
    df_geno_case, df_geno_ctrl = count_genotypes_by_group(args.plink2_path, args.plink2_threads, args.bed_prefix, focus_loci)
    df_hwe_case, df_hwe_ctrl = compute_hwe_by_group(args.plink2_path, args.plink2_threads, args.bed_prefix, focus_loci)
    df_aaf_case, df_aaf_ctrl = compute_aaf_by_group(args.plink2_path, args.plink2_threads, args.bed_prefix, focus_loci)

    LOGGER.info("Annotating loci from ToMMo VCF")
    df_tommo_info = batch_query_variants(args.tommo_vcf_file, focus_loci)

    LOGGER.info("Extracting ANN summaries")
    ann_rows = []
    for variant in focus_loci:
        ann_records = extract_ann_records_from_variant(args.tommo_vcf_file, variant)
        if ann_records:
            summary_lines = [
                f"{rec['Gene_symbol']} ({rec['Transcript_ID']}): {rec['Consequence']} [{rec['Biotype']}]"
                for rec in ann_records
            ]
            summary = "\n".join(f"[{idx + 1}] {line}" for idx, line in enumerate(summary_lines))
            ann_payload = ann_records
        else:
            summary = pd.NA
            ann_payload = pd.NA

        ann_rows.append({"ID": variant, "ANN Summary": summary, "ANN Records": ann_payload})
    df_ann = pd.DataFrame(ann_rows)

    col_maps = {
        "geno": {
            "HOM_REF_CT": "HomRef Count",
            "HET_REF_ALT_CTS": "Het Count",
            "TWO_ALT_GENO_CTS": "HomAlt Count",
            "MISSING_CT": "Miss Count",
            "TOTAL_CT": "Total Count",
            "HOM_REF_FREQ": "HomRef Freq",
            "HET_REF_ALT_FREQ": "Het Freq",
            "TWO_ALT_GENO_FREQ": "HomAlt Freq",
            "MISSING_FREQ": "Miss Freq",
        },
        "aaf": {"ALT_FREQS": "AAF"},
        "hwe": {"P": "HWE P-value"},
    }

    LOGGER.info("Merging all annotation tables")
    merge_sources = [
        (multiallelic_df, "ID"),
        (rename_columns(df_geno_case, "Case", col_maps["geno"]), "ID"),
        (rename_columns(df_geno_ctrl, "Ctrl", col_maps["geno"]), "ID"),
        (rename_columns(df_aaf_case, "Case", col_maps["aaf"]), "ID"),
        (rename_columns(df_aaf_ctrl, "Ctrl", col_maps["aaf"]), "ID"),
        (rename_columns(df_hwe_case, "Case", col_maps["hwe"]), "ID"),
        (rename_columns(df_hwe_ctrl, "Ctrl", col_maps["hwe"]), "ID"),
        (df_tommo_info, "ID"),
        (df_ann, "ID"),
    ]

    for df_source, key in merge_sources:
        merged_df = merged_df.merge(df_source, on=key, how="left")

    merged_df["Case-Ctrl Total Miss Freq"] = np.round(
        (merged_df["Case Miss Count"] + merged_df["Ctrl Miss Count"])
        / (merged_df["Case Total Count"] + merged_df["Ctrl Total Count"]),
        4,
    )

    merged_df = enforce_columns(merged_df, OUTPUT_COLUMNS)
    merged_df.to_csv(args.output, index=False)
    LOGGER.info("Summary written: %s", args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # pylint: disable=broad-except
        LOGGER.exception("Pipeline failed: %s", exc)
        raise SystemExit(1) from exc
