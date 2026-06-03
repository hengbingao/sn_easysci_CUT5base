#!/usr/bin/env python3
"""
snCUT&TAG read labeling — no splitting, one merged FASTQ per I7 sample.

For each read pair, parses i5 ligation / N5 well / Tn5 modality barcodes
from R1, then:
  1. Rewrites the read name to embed all four barcodes
  2. Adds SAM-style tags to the FASTQ comment field (compatible with
     STARsolo --soloFeatures, samtools, and downstream single-cell tools)
  3. Trims Tn5 ME from R1 (immediately after N5 well) and R2 (5' end)

R1 layout:
  [i5 ligation 10bp] [TruSeq Read1 33bp] [Tn5 modality 6-7bp] [N5 well 10bp] [Tn5 ME 19bp] [gDNA ...]

Output read name format:
  @<original_id>_i5=TGAGCC_n5=GTCGCC_i7=TAAGGCGA_tn5=AACACC CB:Z:TGAGCCGTCGCC CR:Z:TGAGCCGTCGCC RG:Z:TAAGGCGA

  - CB:Z  = cell barcode (i5+N5 concatenated, canonical for sc tools)
  - CR:Z  = raw cell barcode (same here; could differ if corrected)
  - RG:Z  = read group / I7 sample barcode
  - XM:Z  = modality tag (Tn5 modality barcode)
  - XI:Z  = full composite identity string (all four barcodes)

This allows:
  - Unified mapping with any aligner (barcodes in qname survive everything)
  - samtools view -d RG:TAAGGCGA to pull one I7 sample from merged BAM
  - STARsolo or custom CB-based grouping using CB tag
  - Full traceability: grep on read name or BAM tag

Usage:
  python label_snCUTTAG.py \
    --r1 sample_i7_001_R1.fastq.gz \
    --r2 sample_i7_001_R2.fastq.gz \
    --i7 TAAGGCGA \
    --i5-barcodes I5_ligation_index.txt \
    --n5-barcodes N5_well_index.txt \
    --tn5-barcodes Tn5_modility.txt \
    --r1-out labeled_R1.fastq.gz \
    --r2-out labeled_R2.fastq.gz \
    [--unmatched-r1 unmatched_R1.fastq.gz]  # optional: save unmatched reads
    [--max-mismatch 1]
    [--min-length 20]
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

# ── Constants ─────────────────────────────────────────────────────────────────
TN5_ME     = "AGATGTGTATAAGAGACAG"   # 19 bp
TN5_ME_RC  = "CTGTCTCTTATACACATCT"   # RC of Tn5 ME (can appear at 3' of R1)
TRUSEQ_R1  = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"  # 33 bp anchor
I5_LIG_LEN = 10


# ── Helpers ───────────────────────────────────────────────────────────────────
def revcomp(seq):
    t = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(t)[::-1]


def hamming(a, b):
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


def best_match(query, barcodes, max_mm):
    best, best_mm = None, max_mm + 1
    for bc in barcodes:
        blen = len(bc)
        if len(query) < blen:
            continue
        mm = hamming(query[:blen], bc)
        if mm < best_mm:
            best_mm, best = mm, bc
    return (best, best_mm) if best is not None and best_mm <= max_mm else (None, None)


def find_anchor(seq, anchor=TRUSEQ_R1, window=20, max_mm=3):
    alen = len(anchor)
    best_pos, best_mm = None, max_mm + 1
    for i in range(I5_LIG_LEN, min(I5_LIG_LEN + window, len(seq) - alen + 1)):
        mm = hamming(seq[i:i + alen], anchor)
        if mm < best_mm:
            best_mm, best_pos = mm, i
    return best_pos if best_pos is not None and best_mm <= max_mm else None


def trim_me(seq, qual, end="left", max_mm=2):
    me_len = len(TN5_ME)
    ltrim = rtrim = 0
    if end in ("left", "both") and hamming(seq[:me_len], TN5_ME) <= max_mm:
        ltrim = me_len
    if end in ("right", "both") and len(seq) >= me_len and \
            hamming(seq[-me_len:], TN5_ME_RC) <= max_mm:
        rtrim = me_len
    s = seq[ltrim: len(seq) - rtrim if rtrim else None]
    q = qual[ltrim: len(qual) - rtrim if rtrim else None]
    return s, q


def open_fastq(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def fastq_iter(fh):
    while True:
        name = fh.readline().rstrip()
        if not name:
            break
        yield name, fh.readline().rstrip(), fh.readline().rstrip(), fh.readline().rstrip()


def load_barcodes(path):
    with open(path) as f:
        return [l.strip() for l in f if l.strip()]


# ── Core labeling logic ───────────────────────────────────────────────────────
def label_reads(args):
    i5_bcs  = load_barcodes(args.i5_barcodes)
    n5_bcs  = load_barcodes(args.n5_barcodes)
    tn5_bcs = load_barcodes(args.tn5_barcodes)
    i7      = args.i7.upper()
    mm      = args.max_mismatch

    os.makedirs(os.path.dirname(os.path.abspath(args.r1_out)), exist_ok=True)

    fout1 = open_fastq(args.r1_out, "wt")
    fout2 = open_fastq(args.r2_out, "wt")
    fout_u1 = open_fastq(args.unmatched_r1, "wt") if args.unmatched_r1 else None
    fout_u2 = open_fastq(args.unmatched_r2, "wt") if args.unmatched_r2 else None

    stats = defaultdict(int)

    with open_fastq(args.r1) as f1, open_fastq(args.r2) as f2:
        for (n1, s1, p1, q1), (n2, s2, p2, q2) in zip(fastq_iter(f1), fastq_iter(f2)):
            stats["total"] += 1

            # Strip existing comment from name (everything after first space)
            read_id = n1[1:].split()[0]  # drop '@', take first token

            # ── Parse barcodes ────────────────────────────────────────────
            i5_bc, _ = best_match(s1[:I5_LIG_LEN], i5_bcs, mm)
            if i5_bc is None:
                stats["fail_i5"] += 1
                if fout_u1:
                    fout_u1.write(f"{n1}\n{s1}\n{p1}\n{q1}\n")
                    fout_u2.write(f"{n2}\n{s2}\n{p2}\n{q2}\n")
                continue

            # R1 layout after i5:
            #   [TruSeq Read1 33bp] [Tn5 modality 6-7bp] [N5 well 10bp] [Tn5 ME 19bp] [gDNA]
            anchor_pos = find_anchor(s1) or I5_LIG_LEN
            tn5_start  = anchor_pos + len(TRUSEQ_R1)
            tn5_len    = max(len(b) for b in tn5_bcs) if tn5_bcs else 9

            tn5_bc, _ = best_match(s1[tn5_start: tn5_start + tn5_len], tn5_bcs, mm)
            if tn5_bc is None:
                stats["fail_tn5"] += 1
                if fout_u1:
                    fout_u1.write(f"{n1}\n{s1}\n{p1}\n{q1}\n")
                    fout_u2.write(f"{n2}\n{s2}\n{p2}\n{q2}\n")
                continue

            n5_start = tn5_start + len(tn5_bc)
            n5_len   = len(n5_bcs[0]) if n5_bcs else 10

            n5_bc, _ = best_match(s1[n5_start: n5_start + n5_len], n5_bcs, mm)
            if n5_bc is None:
                stats["fail_n5"] += 1
                if fout_u1:
                    fout_u1.write(f"{n1}\n{s1}\n{p1}\n{q1}\n")
                    fout_u2.write(f"{n2}\n{s2}\n{p2}\n{q2}\n")
                continue

            # ── Trim Tn5 ME ───────────────────────────────────────────────
            # ME sits immediately after N5 well
            me_start    = n5_start + n5_len
            me_len      = len(TN5_ME)
            gdna_start  = me_start + me_len \
                          if hamming(s1[me_start:me_start+me_len], TN5_ME) <= 2 \
                          else me_start
            r1_seq  = s1[gdna_start:]
            r1_qual = q1[gdna_start:]
            # also trim any 3' adapter contamination on R1
            r1_seq, r1_qual = trim_me(r1_seq, r1_qual, end="right")
            # R2: Tn5 ME at 5' end
            r2_seq, r2_qual = trim_me(s2, q2, end="left")

            if len(r1_seq) < args.min_length or len(r2_seq) < args.min_length:
                stats["too_short"] += 1
                continue

            # ── Build labeled read name ───────────────────────────────────
            # Composite cell barcode: i5 + N5 (concatenated, like 10x CB)
            cb = i5_bc + n5_bc

            # SAM-style tags in FASTQ comment field:
            #   CB:Z = corrected cell barcode
            #   CR:Z = raw cell barcode (same as CB here)
            #   RG:Z = read group (I7 index — identifies the Illumina sample)
            #   XM:Z = modality (Tn5 modality barcode)
            #   XI:Z = full composite identity
            xi  = f"i5={i5_bc},n5={n5_bc},i7={i7},tn5={tn5_bc}"
            comment = f"CB:Z:{cb} CR:Z:{cb} RG:Z:{i7} XM:Z:{tn5_bc} XI:Z:{xi}"

            # New read name: original_id + barcode suffix + space + comment
            # The suffix makes every combination greppable even if comment is lost
            new_id = f"{read_id}_i5={i5_bc}_n5={n5_bc}_i7={i7}_tn5={tn5_bc}"
            new_name1 = f"@{new_id} {comment}"
            new_name2 = f"@{new_id} {comment}"  # same name for R2 (paired)

            fout1.write(f"{new_name1}\n{r1_seq}\n{p1}\n{r1_qual}\n")
            fout2.write(f"{new_name2}\n{r2_seq}\n{p2}\n{r2_qual}\n")
            stats["pass"] += 1

    fout1.close(); fout2.close()
    if fout_u1: fout_u1.close(); fout_u2.close()

    total = stats["total"]
    print(f"\n=== Labeling statistics for I7={i7} ===")
    print(f"Total read pairs   : {total:>10,}")
    print(f"Labeled (pass)     : {stats['pass']:>10,}  ({100*stats['pass']/total:.1f}%)")
    print(f"Failed i5          : {stats['fail_i5']:>10,}  ({100*stats['fail_i5']/total:.1f}%)")
    print(f"Failed N5          : {stats['fail_n5']:>10,}  ({100*stats['fail_n5']/total:.1f}%)")
    print(f"Failed Tn5 modality: {stats['fail_tn5']:>10,}  ({100*stats['fail_tn5']/total:.1f}%)")
    print(f"Too short          : {stats['too_short']:>10,}  ({100*stats['too_short']/total:.1f}%)")
    print(f"\nOutputs: {args.r1_out}  {args.r2_out}")


# ── CLI ───────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--r1",            required=True)
    p.add_argument("--r2",            required=True)
    p.add_argument("--i7",            required=True,  help="I7 barcode for this file")
    p.add_argument("--i5-barcodes",   required=True)
    p.add_argument("--n5-barcodes",   required=True)
    p.add_argument("--tn5-barcodes",  required=True)
    p.add_argument("--r1-out",        required=True,  help="Output labeled R1")
    p.add_argument("--r2-out",        required=True,  help="Output labeled R2")
    p.add_argument("--unmatched-r1",  default=None,   help="Optional: save unmatched R1")
    p.add_argument("--unmatched-r2",  default=None,   help="Optional: save unmatched R2")
    p.add_argument("--max-mismatch",  type=int, default=1)
    p.add_argument("--min-length",    type=int, default=20)
    label_reads(p.parse_args())


if __name__ == "__main__":
    main()
