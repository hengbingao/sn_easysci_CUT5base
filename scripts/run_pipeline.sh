#!/usr/bin/env bash
# =============================================================================
# snCUTTAG-label  —  full pipeline from BCL to labeled + merged FASTQs
#
# Steps:
#   1. Generate SampleSheet.csv from I7 barcode list
#   2. Run bcl-convert (I7-level demultiplexing)
#   3. Label each I7 FASTQ pair (embed i5/N5/Tn5 barcodes, trim Tn5 ME)
#   4. Merge all labeled FASTQs into a single R1/R2 pair
#   5. (Optional) Map with bowtie2 and convert to fragments
#
# Usage:
#   bash run_pipeline.sh [options]
#
# Options (override via environment variables or edit below):
#   BCL_DIR         Path to Illumina run folder         [required]
#   GENOME_INDEX    bowtie2 genome index prefix          [required for mapping]
#   FASTQ_DIR       bcl-convert output directory         [./fastq_bcl]
#   LABELED_DIR     Labeled FASTQ output directory       [./labeled]
#   ALIGN_DIR       Mapping output directory             [./aligned]
#   THREADS         CPU threads                          [8]
#   MAX_MISMATCH    Barcode mismatch tolerance           [1]
#   MIN_LENGTH      Min read length after trimming       [20]
#   SKIP_MAPPING    Set to 1 to stop after labeling      [0]
# =============================================================================
set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS="${REPO_DIR}/scripts"
BARCODES="${REPO_DIR}/barcodes"

BCL_DIR="${BCL_DIR:-}"
GENOME_INDEX="${GENOME_INDEX:-}"
FASTQ_DIR="${FASTQ_DIR:-./fastq_bcl}"
LABELED_DIR="${LABELED_DIR:-./labeled}"
ALIGN_DIR="${ALIGN_DIR:-./aligned}"
SAMPLESHEET="${FASTQ_DIR}/SampleSheet.csv"
THREADS="${THREADS:-8}"
MAX_MISMATCH="${MAX_MISMATCH:-1}"
MIN_LENGTH="${MIN_LENGTH:-20}"
SKIP_MAPPING="${SKIP_MAPPING:-0}"

# ── Validate required inputs ──────────────────────────────────────────────────
if [[ -z "${BCL_DIR}" ]]; then
    echo "ERROR: BCL_DIR is not set."
    echo "  Usage: BCL_DIR=/path/to/RunFolder bash run_pipeline.sh"
    exit 1
fi

if [[ ! -d "${BCL_DIR}" ]]; then
    echo "ERROR: BCL_DIR does not exist: ${BCL_DIR}"
    exit 1
fi

if [[ "${SKIP_MAPPING}" == "0" && -z "${GENOME_INDEX}" ]]; then
    echo "WARNING: GENOME_INDEX is not set — skipping mapping step."
    echo "  Set GENOME_INDEX=/path/to/bowtie2/index to enable mapping."
    SKIP_MAPPING=1
fi

# ── Step 1: Generate SampleSheet ─────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════"
echo " Step 1: Generate SampleSheet.csv"
echo "════════════════════════════════════════════════"
mkdir -p "${FASTQ_DIR}"
python3 "${SCRIPTS}/make_samplesheet.py" \
    --i7-list "${BARCODES}/I7_index.txt" \
    --out     "${SAMPLESHEET}"

# ── Step 2: bcl-convert ───────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════"
echo " Step 2: bcl-convert (I7 demultiplexing)"
echo "════════════════════════════════════════════════"
bcl-convert \
    --bcl-input-directory  "${BCL_DIR}" \
    --output-directory     "${FASTQ_DIR}" \
    --sample-sheet         "${SAMPLESHEET}" \
    --no-lane-splitting    true \
    --bcl-num-conversion-threads  "${THREADS}" \
    --bcl-num-compression-threads "${THREADS}"

echo "bcl-convert complete. FASTQs in: ${FASTQ_DIR}"

# ── Step 3: Label each I7 sample ─────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════"
echo " Step 3: Label reads (embed barcodes, trim ME)"
echo "════════════════════════════════════════════════"
mkdir -p "${LABELED_DIR}/per_i7"
mkdir -p "${LABELED_DIR}/unmatched"

N_LABELED=0
for R1 in "${FASTQ_DIR}"/*_R1_001.fastq.gz; do
    [[ -f "${R1}" ]] || continue
    SAMPLE=$(basename "${R1}" _R1_001.fastq.gz)
    R2="${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz"

    # Extract I7 barcode from sample name (format: sample_i7_NNN_BARCODE)
    I7=$(echo "${SAMPLE}" | grep -oP '[ACGT]{8}$' || true)
    if [[ -z "${I7}" ]]; then
        echo "  WARNING: Cannot extract I7 from '${SAMPLE}', skipping."
        continue
    fi

    echo "  [${SAMPLE}]  I7=${I7}"
    python3 "${SCRIPTS}/label_snCUTTAG.py" \
        --r1           "${R1}" \
        --r2           "${R2}" \
        --i7           "${I7}" \
        --i5-barcodes  "${BARCODES}/I5_ligation_index.txt" \
        --n5-barcodes  "${BARCODES}/N5_well_index.txt" \
        --tn5-barcodes "${BARCODES}/Tn5_modility.txt" \
        --r1-out       "${LABELED_DIR}/per_i7/${SAMPLE}_labeled_R1.fastq.gz" \
        --r2-out       "${LABELED_DIR}/per_i7/${SAMPLE}_labeled_R2.fastq.gz" \
        --unmatched-r1 "${LABELED_DIR}/unmatched/${SAMPLE}_unmatched_R1.fastq.gz" \
        --unmatched-r2 "${LABELED_DIR}/unmatched/${SAMPLE}_unmatched_R2.fastq.gz" \
        --max-mismatch "${MAX_MISMATCH}" \
        --min-length   "${MIN_LENGTH}"

    N_LABELED=$((N_LABELED + 1))
done

if [[ "${N_LABELED}" -eq 0 ]]; then
    echo "ERROR: No I7 FASTQ files were labeled. Check FASTQ_DIR and sample naming."
    exit 1
fi
echo "Labeled ${N_LABELED} I7 sample(s)."

# ── Step 4: Merge ─────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════"
echo " Step 4: Merge labeled FASTQs"
echo "════════════════════════════════════════════════"
MERGED_R1="${LABELED_DIR}/merged_R1.fastq.gz"
MERGED_R2="${LABELED_DIR}/merged_R2.fastq.gz"

cat "${LABELED_DIR}"/per_i7/*_labeled_R1.fastq.gz > "${MERGED_R1}"
cat "${LABELED_DIR}"/per_i7/*_labeled_R2.fastq.gz > "${MERGED_R2}"

echo "Merged R1: ${MERGED_R1}  ($(du -sh "${MERGED_R1}" | cut -f1))"
echo "Merged R2: ${MERGED_R2}  ($(du -sh "${MERGED_R2}" | cut -f1))"

# ── Step 5: Map + fragments (optional) ───────────────────────────────────────
if [[ "${SKIP_MAPPING}" == "1" ]]; then
    echo ""
    echo "Skipping mapping (SKIP_MAPPING=1 or no GENOME_INDEX set)."
    echo "Run manually:"
    echo ""
    echo "  bowtie2 -x \${GENOME_INDEX} \\"
    echo "    -1 ${MERGED_R1} \\"
    echo "    -2 ${MERGED_R2} \\"
    echo "    --no-mixed --no-discordant -X 2000 -p ${THREADS} \\"
    echo "  | samtools sort -@ ${THREADS} -o aligned/merged.bam"
    echo ""
    echo "  python3 scripts/bam_to_fragments.py \\"
    echo "    --bam aligned/merged.bam \\"
    echo "    --out aligned/fragments.tsv.gz \\"
    echo "    --meta aligned/barcode_metadata.tsv"
else
    echo ""
    echo "════════════════════════════════════════════════"
    echo " Step 5: Map with bowtie2"
    echo "════════════════════════════════════════════════"
    mkdir -p "${ALIGN_DIR}"
    BAM="${ALIGN_DIR}/merged.bam"

    bowtie2 \
        -x "${GENOME_INDEX}" \
        -1 "${MERGED_R1}" \
        -2 "${MERGED_R2}" \
        --no-mixed --no-discordant \
        -X 2000 \
        -p "${THREADS}" \
      | samtools sort -@ "${THREADS}" -o "${BAM}"

    samtools index "${BAM}"
    echo "BAM: ${BAM}"

    echo ""
    echo "════════════════════════════════════════════════"
    echo " Step 6: BAM → fragments"
    echo "════════════════════════════════════════════════"
    python3 "${SCRIPTS}/bam_to_fragments.py" \
        --bam  "${BAM}" \
        --out  "${ALIGN_DIR}/fragments.tsv.gz" \
        --meta "${ALIGN_DIR}/barcode_metadata.tsv"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════"
echo " Pipeline complete"
echo "════════════════════════════════════════════════"
echo ""
echo "Key outputs:"
echo "  ${MERGED_R1}"
echo "  ${MERGED_R2}"
[[ "${SKIP_MAPPING}" == "0" ]] && echo "  ${ALIGN_DIR}/merged.bam"
[[ "${SKIP_MAPPING}" == "0" ]] && echo "  ${ALIGN_DIR}/fragments.tsv.gz"
[[ "${SKIP_MAPPING}" == "0" ]] && echo "  ${ALIGN_DIR}/barcode_metadata.tsv"
echo ""
