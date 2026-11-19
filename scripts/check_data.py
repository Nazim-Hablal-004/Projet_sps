#!/usr/bin/env python3
import argparse
import gzip
import os
import sys

REQUIRED = ["GeneSymbol", "Type", "ClinicalSignificance", "PhenotypeList", "Assembly"]


def read_header(path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        line = f.readline()
    return line.rstrip("\n").lstrip("#").split("\t")


def main():
    parser = argparse.ArgumentParser(
        description="Check that the ClinVar file has the expected columns."
    )
    parser.add_argument("--input", required=True, help="path to variant_summary.txt.gz")
    parser.add_argument("--output", default="data/clean/colonnes.txt",
                        help="where to write the column numbers")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"ERROR: file not found: {args.input}", file=sys.stderr)
        sys.exit(3)

    header = read_header(args.input)

    missing = [c for c in REQUIRED if c not in header]
    if missing:
        print(f"ERROR: missing columns: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as out:
        for col in REQUIRED:
            number = header.index(col) + 1
            print(f"{col}: {number}")
            out.write(f"{col}\t{number}\n")

    print("OK: all required columns found")


if __name__ == "__main__":
    main()