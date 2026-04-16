#!/usr/bin/env zsh
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <bim_file> <high_ld_file> <output_file>" >&2
    exit 1
fi

bim_file="$1"
high_ld_file="$2"
out_file="$3"

# Detect chromosome labels from BIM and High-LD files, then normalize High-LD to match BIM.
bim_chr=$(head -n 1 "$bim_file" | awk '{print $1}')
ld_chr=$(head -n 1 "$high_ld_file" | awk '{print $1}')

cp "$high_ld_file" "$out_file"

if [[ "$bim_chr" == chr* ]] && [[ "$ld_chr" != chr* ]]; then
    echo "BIM has chr prefix, High-LD does not. Adding chr prefix to High-LD file."
    awk '{print "chr"$0}' "$high_ld_file" > "$out_file"
elif [[ "$bim_chr" != chr* ]] && [[ "$ld_chr" == chr* ]]; then
    echo "BIM does not have chr prefix, High-LD has it. Removing chr prefix from High-LD file."
    sed 's/^chr//' "$high_ld_file" > "$out_file"
else
    echo "Chromosome formats already match."
fi
