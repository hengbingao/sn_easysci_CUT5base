#!/usr/bin/env python3
"""
Generate a bcl-convert SampleSheet.csv from an I7 barcode list.

Usage:
  python make_samplesheet.py \
    --i7-list barcodes/I7_index.txt \
    --out SampleSheet.csv \
    [--software-version 3.9.3]
"""

import argparse
import os


def make_samplesheet(args):
    with open(args.i7_list) as f:
        barcodes = [l.strip() for l in f if l.strip()]

    if not barcodes:
        raise ValueError(f"No barcodes found in {args.i7_list}")

    lines = [
        "[Header]",
        "FileFormatVersion,2",
        "",
        "[BCLConvert_Settings]",
        f"SoftwareVersion,{args.software_version}",
        "FastqCompressionFormat,gzip",
        "",
        "[BCLConvert_Data]",
        "Sample_ID,index",
    ]
    for i, bc in enumerate(barcodes, 1):
        lines.append(f"sample_i7_{i:03d}_{bc},{bc}")

    with open(args.out, "w") as fout:
        fout.write("\n".join(lines) + "\n")

    print(f"Written: {args.out}  ({len(barcodes)} I7 samples)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--i7-list",          required=True, help="I7 barcode list, one per line")
    p.add_argument("--out",              default="SampleSheet.csv")
    p.add_argument("--software-version", default="3.9.3")
    make_samplesheet(p.parse_args())


if __name__ == "__main__":
    main()
