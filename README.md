# easysci_CUT5base-label

A lightweight pipeline for processing snCUT&TAG libraries that **labels every read with its full barcode identity instead of splitting into per-cell FASTQ files**. All four index layers — i5 ligation, N5 well, I7 sample, and Tn5 modality — are embedded in the read name and FASTQ comment field, enabling a single mapping job followed by flexible BAM-level grouping.

---

## Table of contents

- [Overview](#overview)
- [Library structure](#library-structure)
- [How it works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Step-by-step usage](#step-by-step-usage)
  - [Step 1 — Generate SampleSheet](#step-1--generate-samplesheet)
  - [Step 2 — bcl-convert](#step-2--bcl-convert)
  - [Step 3 — Label reads](#step-3--label-reads)
  - [Step 4 — Merge and map](#step-4--merge-and-map)
  - [Step 5 — BAM to fragments](#step-5--bam-to-fragments)
- [Read name format](#read-name-format)
- [Downstream analysis](#downstream-analysis)
- [Script reference](#script-reference)
- [Repository structure](#repository-structure)
- [FAQ](#faq)

---

## Overview

```
BCL run folder
      │
      ▼
make_samplesheet.py ──► SampleSheet.csv
      │
      ▼
bcl-convert                          I7 index (8 bp)
      │                              → one FASTQ pair per I7 barcode
      ▼
label_snCUTTAG.py  ◄── I5 ligation, N5 well, Tn5 modality barcode lists
      │
      │  For every read pair in R1:
      │    1. Match i5 ligation barcode   (R1 bp  0–9)
      │    2. Locate TruSeq Read1 anchor  (R1 bp 10–42, 33 bp)
      │    3. Match N5 well barcode       (R1 bp 43–52)
      │    4. Match Tn5 modality barcode  (R1 bp 53–58/61)
      │    5. Trim Tn5 ME from R1 and R2
      │    6. Write labeled FASTQ with barcodes in read name + comment tags
      │
      ▼
labeled/per_i7/  ──► cat ──► labeled/merged_R1.fastq.gz
                              labeled/merged_R2.fastq.gz
                                      │
                                      ▼
                              bowtie2 / bwa mem -C
                                      │
                                      ▼
                              merged.bam  (all reads, all barcodes preserved)
                                      │
                        ┌─────────────┴──────────────┐
                        ▼                            ▼
                samtools view                bam_to_fragments.py
                filter by tag               fragments.tsv.gz
                (RG, XM, CB, XI)            barcode_metadata.tsv
```

---

## Library structure

```
Read 1 (5′ → 3′):
┌──────────────┬───────────────┬───────────────┬──────────────┬─────────┬─────────────┐
│ i5 ligation  │ TruSeq Read 1 │ Tn5 modality  │ N5 well idx  │ Tn5 ME  │ genomic DNA │
│    10 bp     │    33 bp      │    6–7 bp     │    10 bp     │  19 bp  │      →      │
└──────────────┴───────────────┴───────────────┴──────────────┴─────────┴─────────────┘

I7 index read:
┌──────────────┐
│  I7 barcode  │  8 bp  — handled by bcl-convert
└──────────────┘

Read 2 (5′ → 3′):
┌─────────┬─────────────┐
│ Tn5 ME  │ genomic DNA │
│  19 bp  │      →      │
└─────────┴─────────────┘
```

The TruSeq Read 1 sequence (`ACACTCTTTCCCTACACGACGCTCTTCCGATCT`) serves as an anchor to precisely locate the N5 well index and Tn5 modality barcode within R1, tolerating up to 3 mismatches. The Tn5 ME sequence (`AGATGTGTATAAGAGACAG`) is trimmed from both R1 (between the barcode region and genomic DNA) and from the 5′ end of R2.

---

## How it works

Traditional snCUT&TAG demultiplexing splits reads into thousands of per-cell FASTQ files — one per unique barcode combination. This creates storage overhead, requires running the aligner once per file, and makes traceability depend entirely on the filename.

This pipeline takes a different approach: instead of splitting, it **labels every read** with its full barcode identity and then maps everything together.

| | Split approach | **Label approach (this pipeline)** |
|---|---|---|
| FASTQ files | thousands | **2 merged files** |
| Mapping jobs | one per split | **one** |
| Barcode traceability | filename only | **read name + BAM tags** |
| Downstream grouping | re-split required | **filter BAM by any tag** |
| Storage | redundant headers everywhere | minimal overhead |
| Tool compatibility | universal | universal |

---

## Requirements

| Tool | Version | Notes |
|---|---|---|
| Python | ≥ 3.8 | standard library only for core scripts |
| bcl-convert | ≥ 3.9 | Illumina; for Step 2 |
| bowtie2 | ≥ 2.4 | or `bwa mem`; for mapping |
| samtools | ≥ 1.15 | for BAM processing |
| pysam | ≥ 0.22 | Python package; for `bam_to_fragments.py` only |

---

## Installation

```bash
git clone https://github.com/yourlab/snCUTTAG-label.git
cd snCUTTAG-label

# Install pysam for the optional bam_to_fragments.py step
pip install -r requirements.txt
```

No build step is required. The core labeling scripts run with Python ≥ 3.8 and no third-party dependencies.

---

## Quick start

```bash
# One command runs the full pipeline (BCL → labeled merged FASTQs → BAM → fragments)
BCL_DIR=/path/to/RunFolder \
GENOME_INDEX=/path/to/bowtie2/index \
bash scripts/run_pipeline.sh
```

Set `SKIP_MAPPING=1` to stop after producing merged FASTQs without mapping:

```bash
BCL_DIR=/path/to/RunFolder SKIP_MAPPING=1 bash scripts/run_pipeline.sh
```

---

## Step-by-step usage

### Step 1 — Generate SampleSheet

```bash
python3 scripts/make_samplesheet.py \
    --i7-list barcodes/I7_index.txt \
    --out     SampleSheet.csv
```

This creates a bcl-convert–compatible `SampleSheet.csv` with one entry per I7 barcode, naming samples as `sample_i7_NNN_BARCODE`.

---

### Step 2 — bcl-convert

```bash
bcl-convert \
    --bcl-input-directory  /path/to/RunFolder \
    --output-directory     ./fastq_bcl \
    --sample-sheet         SampleSheet.csv \
    --no-lane-splitting    true \
    --bcl-num-conversion-threads  8 \
    --bcl-num-compression-threads 8
```

This produces one R1/R2 FASTQ pair per I7 barcode. These files are not yet cell-level demultiplexed — they still contain reads from all i5 ligation, N5 well, and Tn5 modality combinations.

---

### Step 3 — Label reads

Run `label_snCUTTAG.py` on each I7 FASTQ pair. The script parses the remaining three barcode layers from R1, trims Tn5 ME from both reads, and writes labeled FASTQs.

**Single sample:**

```bash
python3 scripts/label_snCUTTAG.py \
    --r1           fastq_bcl/sample_i7_001_TAAGGCGA_R1_001.fastq.gz \
    --r2           fastq_bcl/sample_i7_001_TAAGGCGA_R2_001.fastq.gz \
    --i7           TAAGGCGA \
    --i5-barcodes  barcodes/I5_ligation_index.txt \
    --n5-barcodes  barcodes/N5_well_index.txt \
    --tn5-barcodes barcodes/Tn5_modility.txt \
    --r1-out       labeled/per_i7/sample_i7_001_TAAGGCGA_labeled_R1.fastq.gz \
    --r2-out       labeled/per_i7/sample_i7_001_TAAGGCGA_labeled_R2.fastq.gz \
    --unmatched-r1 labeled/unmatched/sample_i7_001_TAAGGCGA_unmatched_R1.fastq.gz \
    --unmatched-r2 labeled/unmatched/sample_i7_001_TAAGGCGA_unmatched_R2.fastq.gz \
    --max-mismatch 1 \
    --min-length   20
```

The script prints a per-sample statistics summary on completion:

```
=== Labeling statistics for I7=TAAGGCGA ===
Total read pairs   :  5,000,000
Labeled (pass)     :  4,312,441  (86.2%)
Failed i5          :    412,300  ( 8.2%)
Failed N5          :    198,100  ( 4.0%)
Failed Tn5 modality:     77,159  ( 1.5%)
Too short          :          0  ( 0.0%)
```

---

### Step 4 — Merge and map

```bash
# Merge all labeled FASTQs
cat labeled/per_i7/*_labeled_R1.fastq.gz > labeled/merged_R1.fastq.gz
cat labeled/per_i7/*_labeled_R2.fastq.gz > labeled/merged_R2.fastq.gz

# Map with bowtie2
# Note: bowtie2 automatically passes FASTQ comment fields to SAM optional tags
bowtie2 \
    -x /path/to/genome/index \
    -1 labeled/merged_R1.fastq.gz \
    -2 labeled/merged_R2.fastq.gz \
    --no-mixed --no-discordant \
    -X 2000 \
    -p 16 \
  | samtools sort -@ 8 -o aligned/merged.bam

samtools index aligned/merged.bam
```

> **Using BWA instead of bowtie2?** Add the `-C` flag to pass FASTQ comment fields through to the BAM:
> ```bash
> bwa mem -C -t 16 genome.fa labeled/merged_R1.fastq.gz labeled/merged_R2.fastq.gz \
>   | samtools sort -o aligned/merged.bam
> ```

---

### Step 5 — BAM to fragments

```bash
python3 scripts/bam_to_fragments.py \
    --bam  aligned/merged.bam \
    --out  aligned/fragments.tsv.gz \
    --meta aligned/barcode_metadata.tsv
```

This produces a standard 5-column fragments file and a sidecar barcode metadata table (see [Downstream analysis](#downstream-analysis)).

---

## Read name format

After labeling, each read name carries the full barcode identity in two complementary forms:

```
@{original_id}_{barcodes}  {SAM-style comment tags}
```

**Example:**

```
@NS500123:1:HYFK3BGX3:1:11101:9071:1048_i5=TGAGCCGCGG_n5=GTCGCCAACC_i7=TAAGGCGA_tn5=AACACC CB:Z:TGAGCCGCGGGTCGCCAACC CR:Z:TGAGCCGCGGGTCGCCAACC RG:Z:TAAGGCGA XM:Z:AACACC XI:Z:i5=TGAGCCGCGG,n5=GTCGCCAACC,i7=TAAGGCGA,tn5=AACACC
```

| Field | Content | Purpose |
|---|---|---|
| `_i5=..._n5=..._i7=..._tn5=...` | barcode suffix in qname | survives any aligner; always greppable |
| `CB:Z:` | i5 + N5 concatenated | cell barcode for SnapATAC2, STARsolo, Signac |
| `CR:Z:` | same as CB (raw barcode) | standard single-cell tools convention |
| `RG:Z:` | I7 barcode | `samtools view -r` by Illumina sample |
| `XM:Z:` | Tn5 modality barcode | identifies antibody / modality |
| `XI:Z:` | full identity string | single-field audit trail |

The barcode suffix in the read name is redundant with the comment tags by design: the qname survives even if an aligner strips optional fields, making every read independently traceable.

---

## Downstream analysis

### Filter BAM by modality

```bash
# Extract all reads for a specific Tn5 modality (e.g. tn5=AACACC = H3K27ac)
samtools view -h aligned/merged.bam \
  | awk '/^@/ || $1 ~ /_tn5=AACACC/' \
  | samtools view -bS -o aligned/H3K27ac.bam
```

### Filter BAM by I7 sample

```bash
# Extract one Illumina sample using the RG tag
samtools view -r TAAGGCGA -b aligned/merged.bam -o aligned/sample_TAAGGCGA.bam
```

### Count reads per barcode combination

```bash
samtools view aligned/merged.bam \
  | grep -oP 'XI:Z:\S+' \
  | sort | uniq -c | sort -rn \
  | head -30
```

### Load fragments in SnapATAC2

```python
import snapatac2 as snap
import pandas as pd

# Load fragments — composite barcode in col4, fully compatible
data = snap.pp.import_fragments(
    "aligned/fragments.tsv.gz",
    chrom_sizes=snap.genome.hg38,
)

# Join sidecar metadata to recover individual barcode identities
meta = pd.read_csv("aligned/barcode_metadata.tsv", sep="\t", index_col="composite_bc")
data.obs = data.obs.join(meta)

# Now data.obs has columns: i5, n5, i7, tn5
# Subset by modality
h3k27ac = data[data.obs["tn5"] == "AACACC"]
```

### Fragments file format

The output `fragments.tsv.gz` follows the standard 5-column BED-like format. Column 4 contains the composite barcode with `+` as a separator, which can be split back into individual indices at any time:

```
chr1  1000  1200  TGAGCCGCGG+GTCGCCAACC+TAAGGCGA+AACACC  1
chr1  1500  1700  TGAGCCGCGG+GTCGCCAACC+TAAGGCGA+AACACC  1
chr2  3000  3180  ACGAGGAGCC+AACGATCTAC+CGTACTAG+AGTCCA   1
```

```python
# Split composite barcode back into components
cb = "TGAGCCGCGG+GTCGCCAACC+TAAGGCGA+AACACC"
i5, n5, i7, tn5 = cb.split("+")
```

The accompanying `barcode_metadata.tsv` maps every observed composite barcode to its component indices:

```
composite_bc                              i5          n5          i7        tn5
TGAGCCGCGG+GTCGCCAACC+TAAGGCGA+AACACC   TGAGCCGCGG  GTCGCCAACC  TAAGGCGA  AACACC
```

---

## Script reference

### `make_samplesheet.py`

| Parameter | Default | Description |
|---|---|---|
| `--i7-list` | required | I7 barcode list, one per line |
| `--out` | `SampleSheet.csv` | Output path |
| `--software-version` | `3.9.3` | bcl-convert software version string |

---

### `label_snCUTTAG.py`

| Parameter | Default | Description |
|---|---|---|
| `--r1` | required | Input R1 FASTQ (.fastq or .fastq.gz) |
| `--r2` | required | Input R2 FASTQ |
| `--i7` | required | I7 barcode for this file (8 bp) |
| `--i5-barcodes` | required | i5 ligation barcode list |
| `--n5-barcodes` | required | N5 well barcode list |
| `--tn5-barcodes` | required | Tn5 modality barcode list |
| `--r1-out` | required | Output labeled R1 |
| `--r2-out` | required | Output labeled R2 |
| `--unmatched-r1` | — | Optional: write unmatched R1 reads here |
| `--unmatched-r2` | — | Optional: write unmatched R2 reads here |
| `--max-mismatch` | `1` | Max mismatches allowed per barcode |
| `--min-length` | `20` | Min read length after trimming (bp) |

---

### `bam_to_fragments.py`

| Parameter | Default | Description |
|---|---|---|
| `--bam` | required | Input BAM (coordinate-sorted, indexed) |
| `--out` | `fragments.tsv.gz` | Output fragments file |
| `--meta` | `barcode_metadata.tsv` | Output barcode metadata TSV |
| `--min-mapq` | `30` | Minimum mapping quality |
| `--max-fraglen` | `2000` | Maximum fragment length (bp) |
| `--shift-plus` | `0` | Tn5 cut-site shift on + strand (use `4` for ATAC-style) |
| `--shift-minus` | `0` | Tn5 cut-site shift on − strand (use `5` for ATAC-style) |

---

### `run_pipeline.sh`

Controlled via environment variables:

| Variable | Default | Description |
|---|---|---|
| `BCL_DIR` | required | Path to Illumina run folder |
| `GENOME_INDEX` | — | bowtie2 index prefix (required for mapping) |
| `FASTQ_DIR` | `./fastq_bcl` | bcl-convert output directory |
| `LABELED_DIR` | `./labeled` | Labeled FASTQ output directory |
| `ALIGN_DIR` | `./aligned` | Mapping output directory |
| `THREADS` | `8` | CPU threads |
| `MAX_MISMATCH` | `1` | Barcode mismatch tolerance |
| `MIN_LENGTH` | `20` | Min read length after trimming |
| `SKIP_MAPPING` | `0` | Set to `1` to stop after producing merged FASTQs |

---

## Repository structure

```
snCUTTAG-label/
├── scripts/
│   ├── make_samplesheet.py     # Generate bcl-convert SampleSheet from I7 list
│   ├── label_snCUTTAG.py       # Core: embed barcodes, trim Tn5 ME
│   ├── bam_to_fragments.py     # BAM → fragments.tsv.gz + metadata
│   └── run_pipeline.sh         # End-to-end pipeline wrapper
├── barcodes/
│   ├── I7_index.txt            # I7 sample barcodes (8 bp, 192 barcodes)
│   ├── I5_ligation_index.txt   # i5 ligation barcodes (10 bp)
│   ├── N5_well_index.txt       # N5 well barcodes (10 bp)
│   └── Tn5_modility.txt        # Tn5 modality barcodes (6–7 bp)
├── requirements.txt
├── CHANGELOG.md
└── README.md
```

---

## FAQ

**Q: Does labeling increase file size?**
Each read name grows by approximately 50 characters. For a typical experiment this adds roughly 5–8% to uncompressed size, which compresses back to ~1–2% overhead in `.fastq.gz`.

**Q: What happens to reads that fail barcode matching?**
Reads that fail any barcode match (i5, N5, or Tn5) are excluded from the labeled output and optionally written to separate unmatched FASTQ files via `--unmatched-r1` / `--unmatched-r2`. Per-sample statistics are printed at the end of each run.

**Q: What is the Tn5 ME sequence?**
`AGATGTGTATAAGAGACAG` (19 bp). Its reverse complement `CTGTCTCTTATACACATCT` is additionally trimmed from the 3′ end of R1 when present as adapter contamination.

**Q: Can I use BWA instead of bowtie2?**
Yes. Use `bwa mem -C` to pass FASTQ comment fields through to the BAM as optional tags:
```bash
bwa mem -C -t 16 genome.fa labeled/merged_R1.fastq.gz labeled/merged_R2.fastq.gz \
  | samtools sort -o aligned/merged.bam
```

**Q: Can I use this with STARsolo?**
The `CB:Z:` tag in the FASTQ comment is directly compatible with STARsolo. Provide `--soloCBwhitelist` with all valid i5+N5 combinations and point `--readFilesIn` at the merged FASTQs.

**Q: How do I reconstruct per-modality fragments files after mapping?**
Filter the BAM by the `XM:Z:` tag (Tn5 modality) and then run `bam_to_fragments.py` on the filtered BAM, or use the `barcode_metadata.tsv` to subset an existing fragments file in Python or R.
