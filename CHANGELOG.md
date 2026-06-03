# Changelog

## [0.1.0] — initial release

- `make_samplesheet.py` — generate bcl-convert SampleSheet from I7 barcode list
- `label_snCUTTAG.py` — embed i5/N5/I7/Tn5 barcodes into read names and FASTQ comment tags; trim Tn5 ME
- `bam_to_fragments.py` — convert labeled BAM to fragments format with composite barcode and sidecar metadata
- `run_pipeline.sh` — end-to-end pipeline from BCL to merged FASTQ (and optionally to fragments)
- Bundled barcode lists for I7, I5 ligation, N5 well, and Tn5 modality
