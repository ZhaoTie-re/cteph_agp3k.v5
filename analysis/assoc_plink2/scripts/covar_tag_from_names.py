#!/usr/bin/env python3
"""Build covariate tag string from Nextflow covariate list.

Example:
- SEX,PC1_AVG-PC10_AVG -> sex.10pc
- PC1_AVG,PC2_AVG -> nosex.2pc
"""

from __future__ import annotations

import argparse
import re


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert covariate names to compact output tag")
    parser.add_argument("--covar-names", required=True, help="Comma-separated covariate names")
    return parser.parse_args()


def build_covar_tag(covar_names: str) -> str:
    tokens = [item.strip() for item in covar_names.split(",") if item.strip()]
    has_sex = any(token.upper() == "SEX" for token in tokens)

    pc_indices: set[int] = set()
    range_pattern = re.compile(r"^PC(\d+)(?:_AVG)?-PC(\d+)(?:_AVG)?$", re.IGNORECASE)
    single_pattern = re.compile(r"^PC(\d+)(?:_AVG)?$", re.IGNORECASE)

    for token in tokens:
        range_match = range_pattern.match(token)
        if range_match:
            start = int(range_match.group(1))
            end = int(range_match.group(2))
            lo, hi = sorted((start, end))
            pc_indices.update(range(lo, hi + 1))
            continue

        single_match = single_pattern.match(token)
        if single_match:
            pc_indices.add(int(single_match.group(1)))

    sex_tag = "sex" if has_sex else "nosex"
    pc_tag = f"{len(pc_indices)}pc" if pc_indices else "nopc"
    return f"{sex_tag}.{pc_tag}"


def main() -> None:
    args = parse_args()
    print(build_covar_tag(args.covar_names))


if __name__ == "__main__":
    main()
