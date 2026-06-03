#!/usr/bin/env python3
"""
BAM → fragments file, preserving all index information.

Reads a BAM produced from labeled FASTQs (read names contain barcode suffixes,
and CB/RG/XM/XI tags are present). Outputs:

1. fragments.tsv.gz  — standard 5-column fragments file, but column 4 is a
   composite barcode in the form:
       {i5}+{n5}+{i7}+{tn5}
   e.g.  TGAGCCGCGG+GTCGCCAACC+TAAGGCGA+AACACC
   This is directly usable by SnapATAC2, ArchR, Signac (they treat col4 as
   an opaque string for grouping). The '+' separator lets you split it back
   at any time.

2. barcode_metadata.tsv  — sidecar mapping table:
   composite_bc  i5  n5  i7  tn5
   One row per unique composite barcode observed. Join this onto any
   per-barcode result table to recover individual index identities.

Fragment definition (same as 10x Cell Ranger):
  - Proper paired-end, both mates mapped, not duplicate
  - start = min(pos_r1, pos_r2),  end = max(end_r1, end_r2)
  - Both mates on same chromosome
  - Fragment length filter: 1 bp – 2000 bp (configurable)

Usage:
  python bam_to_fragments.py \\
    --bam  merged.bam \\
    --out  fragments.tsv.gz \\
    --meta barcode_metadata.tsv \\
    [--min-mapq 30] \\
    [--max-fraglen 2000] \\
    [--shift-plus 4] [--shift-minus 5]   # Tn5 cut-site shift (ATAC style)
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam is required.  pip install pysam --break-system-packages")


# ── Parse composite barcode from read ────────────────────────────────────────

_BC_RE = re.compile(r'_i5=([^_]+)_n5=([^_]+)_i7=([^_]+)_tn5=([^_\s/]+)')

def parse_barcodes(read):
    """
    Extract i5, n5, i7, tn5 from a labeled read.
    Tries BAM tags first (XI:Z), then falls back to qname regex.
    Returns (i5, n5, i7, tn5) or None.
    """
    # Try XI tag (full identity string: i5=...,n5=...,i7=...,tn5=...)
    try:
        xi = read.get_tag("XI")
        parts = dict(kv.split("=") for kv in xi.split(","))
        return parts["i5"], parts["n5"], parts["i7"], parts["tn5"]
    except (KeyError, ValueError):
        pass

    # Fallback: parse read name
    m = _BC_RE.search(read.query_name or "")
    if m:
        return m.group(1), m.group(2), m.group(3), m.group(4)

    # Last resort: try CB + RG + XM tags separately
    try:
        cb  = read.get_tag("CB")  # i5+n5 concatenated
        rg  = read.get_tag("RG")
        xm  = read.get_tag("XM")
        # CB is i5+n5 but we don't know the split point without the original lengths
        # Store as cb_unknown+rg+xm and flag
        return cb, "", rg, xm
    except KeyError:
        pass

    return None


def composite_bc(i5, n5, i7, tn5):
    """Build composite barcode string for fragments col4."""
    return f"{i5}+{n5}+{i7}+{tn5}"


# ── BAM → fragments ───────────────────────────────────────────────────────────

def bam_to_fragments(args):
    bam = pysam.AlignmentFile(args.bam, "rb")

    # Collect proper pairs: key = read name (without /1 /2 suffix)
    # Value = list of reads (expect exactly 2 per pair)
    pairs = defaultdict(list)
    n_total = n_unmapped = n_unpaired = n_nobc = n_lowmapq = n_dup = 0

    print(f"Reading BAM: {args.bam}", file=sys.stderr)
    for read in bam:
        n_total += 1
        if n_total % 5_000_000 == 0:
            print(f"  {n_total:,} reads processed...", file=sys.stderr)

        # Basic filters
        if read.is_unmapped or read.mate_is_unmapped:
            n_unmapped += 1; continue
        if not read.is_proper_pair:
            n_unpaired += 1; continue
        if read.mapping_quality < args.min_mapq:
            n_lowmapq += 1; continue
        if read.is_duplicate:
            n_dup += 1; continue
        if read.is_secondary or read.is_supplementary:
            continue

        # Strip /1 /2 from read name for pairing
        rname = re.sub(r'[/\s][12]$', '', read.query_name or "")
        pairs[rname].append(read)

    bam.close()

    # Process pairs
    fragments = []
    meta = {}   # composite_bc -> (i5, n5, i7, tn5)
    n_written = n_nobc = n_diffchrom = n_fraglen = 0

    print(f"Processing {len(pairs):,} read name groups...", file=sys.stderr)
    for rname, reads in pairs.items():
        if len(reads) != 2:
            continue

        r1, r2 = reads if reads[0].is_read1 else reads[::-1]

        # Chromosome filter
        if r1.reference_name != r2.reference_name:
            n_diffchrom += 1; continue

        # Barcode
        bc_tuple = parse_barcodes(r1) or parse_barcodes(r2)
        if bc_tuple is None:
            n_nobc += 1; continue
        i5, n5, i7, tn5 = bc_tuple
        cb = composite_bc(i5, n5, i7, tn5)
        if cb not in meta:
            meta[cb] = (i5, n5, i7, tn5)

        # Fragment coordinates (0-based, half-open like BED)
        # Apply Tn5 cut-site shift for ATAC-seq (optional, default 0 for CUT&TAG)
        pos1 = r1.reference_start + (args.shift_plus if not r1.is_reverse else -args.shift_minus)
        pos2 = r2.reference_end   + (-args.shift_minus if r2.is_reverse else args.shift_plus)
        frag_start = min(pos1, pos2)
        frag_end   = max(pos1, pos2)

        # Fragment length filter
        flen = frag_end - frag_start
        if flen < 1 or flen > args.max_fraglen:
            n_fraglen += 1; continue

        fragments.append((r1.reference_name, frag_start, frag_end, cb))
        n_written += 1

    # Sort fragments: chr, start, end
    print("Sorting fragments...", file=sys.stderr)
    fragments.sort(key=lambda x: (x[0], x[1], x[2]))

    # Write fragments.tsv.gz
    print(f"Writing {args.out} ...", file=sys.stderr)
    opener = gzip.open if args.out.endswith(".gz") else open
    with opener(args.out, "wt") as fout:
        for chrom, start, end, cb in fragments:
            # score column = 1 (can be changed to read count if desired)
            fout.write(f"{chrom}\t{start}\t{end}\t{cb}\t1\n")

    # Write sidecar metadata table
    print(f"Writing metadata: {args.meta} ...", file=sys.stderr)
    with open(args.meta, "w") as fmeta:
        fmeta.write("composite_bc\ti5\tn5\ti7\ttn5\n")
        for cb, (i5, n5, i7, tn5) in sorted(meta.items()):
            fmeta.write(f"{cb}\t{i5}\t{n5}\t{i7}\t{tn5}\n")

    # Stats
    print(f"\n=== BAM → fragments statistics ===", file=sys.stderr)
    print(f"Total reads          : {n_total:>12,}", file=sys.stderr)
    print(f"  unmapped / unpaired: {n_unmapped+n_unpaired:>12,}", file=sys.stderr)
    print(f"  low MAPQ (<{args.min_mapq})     : {n_lowmapq:>12,}", file=sys.stderr)
    print(f"  duplicates         : {n_dup:>12,}", file=sys.stderr)
    print(f"Pairs processed      : {len(pairs):>12,}", file=sys.stderr)
    print(f"  diff chromosome    : {n_diffchrom:>12,}", file=sys.stderr)
    print(f"  no barcode found   : {n_nobc:>12,}", file=sys.stderr)
    print(f"  fragment len filter: {n_fraglen:>12,}", file=sys.stderr)
    print(f"Fragments written    : {n_written:>12,}", file=sys.stderr)
    print(f"Unique barcodes      : {len(meta):>12,}", file=sys.stderr)
    print(f"\nOutput fragments : {args.out}", file=sys.stderr)
    print(f"Output metadata  : {args.meta}", file=sys.stderr)
    print(f"\nComposite barcode format in col4: {{i5}}+{{n5}}+{{i7}}+{{tn5}}", file=sys.stderr)
    print(f"Split back with: cb.split('+')", file=sys.stderr)


# ── CLI ───────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam",         required=True,  help="Input BAM (sorted, indexed)")
    p.add_argument("--out",         default="fragments.tsv.gz", help="Output fragments file")
    p.add_argument("--meta",        default="barcode_metadata.tsv",
                   help="Output sidecar barcode metadata TSV")
    p.add_argument("--min-mapq",    type=int, default=30,  help="Min mapping quality (default 30)")
    p.add_argument("--max-fraglen", type=int, default=2000, help="Max fragment length bp (default 2000)")
    p.add_argument("--shift-plus",  type=int, default=0,
                   help="Tn5 shift on + strand in bp (default 0; use 4 for ATAC-style)")
    p.add_argument("--shift-minus", type=int, default=0,
                   help="Tn5 shift on - strand in bp (default 0; use 5 for ATAC-style)")
    bam_to_fragments(p.parse_args())

if __name__ == "__main__":
    main()
