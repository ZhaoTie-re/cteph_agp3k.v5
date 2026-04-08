#!/usr/bin/env python3
"""Compare BBJ preprocess signatures with numeric equivalence.

Exit code:
- 0: signatures match
- 1: signatures do not match
"""

from __future__ import annotations

import math
import sys
from pathlib import Path


def parse_signature(path: str) -> dict[str, str]:
    text = Path(path).read_text(encoding="utf-8").strip()
    out: dict[str, str] = {}
    for item in text.split("|"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        out[key.strip()] = value.strip()
    return out


def main() -> int:
    if len(sys.argv) != 3:
        return 1

    expected = parse_signature(sys.argv[1])
    found = parse_signature(sys.argv[2])

    for key in ["mind", "geno", "maf", "hwe"]:
        if key not in expected or key not in found:
            return 1
        try:
            if not math.isclose(float(expected[key]), float(found[key]), rel_tol=0.0, abs_tol=1e-12):
                return 1
        except ValueError:
            if expected[key] != found[key]:
                return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
