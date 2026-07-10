# SPACE-Tag-HD Snakemake pipeline

FASTQ → 40 µm spatial-chromatin **fragment files** for SPACE-Tag-HD (spatial CUT&Tag on
10x Visium HD). A thin Snakemake wrapper around **Space Ranger** plus per-spot deduplication,
fragment-file generation, library-complexity, and a fragments-per-spot QC map.

Works for both **Illumina** and **Ultima V120** data (both delivered as R1 = spatial barcode,
R2 = genomic read). Optional in-pipeline trimming picks the correct adapter set per platform.

> ⚠️ Paths in `Snakefile` and `run_pipeline.sh` are hard-coded for the **NYGC cluster**
> (micromamba at `/gpfs/commons/home/cyan/...`, Space Ranger via `module load`, mm10 GEX
> reference). Edit those for another environment.

---

## What it produces

Per sample, under `<fastq_folder>/<spaceranger_id>_<binstr>/` (e.g. `.../2604il_27ac_040um/`):

| Output | Description |
|---|---|
| `outs/possorted_genome_bam.bam` | Space Ranger alignment (STAR) + `barcode_mappings.parquet` |
| `possorted_genome_bam_deduplicated.bam` + `markdup_stats.txt` | per-spot dedup (`samtools markdup --barcode-tag CB`) |
| `possorted_genome_bam_deduplicated_fragments2.bed.gz` (+ `.tbi`) | **2 µm** fragment file |
| `possorted_genome_bam_deduplicated_fragments_040um.bed.gz` (+ `.tbi`) | **40 µm binned** fragment file *(main downstream input)* |
| `possorted_genome_bam_deduplicated_fragments_CELLS.bed.gz` | cell-level fragments *(only if `run_cells: true`)* |
| `50M_complexity_curve.txt` | preseq `lc_extrap` on a 50 M-read subsample |
| `qc_fragspot_stats.txt` + `qc_fragspot_040um.pdf/.png` | fragments-per-spot summary + spatial map |

## Rules (DAG)

```
[trim*] → spaceranger → deduplicate → fragments_2um → convert_bin → qc_fragspot
                     ↘ complexity                   ↘ convert_to_cells*
```
`*` optional. `trim` runs only when `trim: true`; `convert_to_cells` only when `run_cells: true`.

- **trim** *(optional)* — cutadapt on every raw lane FASTQ (`R2` 3′ adapters; handles lane-split
  input), writing trimmed FASTQs into `fastq_folder` with lane names preserved. Adapter set is
  chosen by `platform` (see below).
- **spaceranger** — `spaceranger count` (Visium HD, `--custom-bin-size`, `--create-bam`).
- **complexity** — subsample BAM to 50 M reads → `preseq lc_extrap`.
- **deduplicate** — name-sort → `fixmate -m` → coord-sort → `markdup -r --barcode-tag CB`
  (same coordinate + different spot barcode is kept; = per-spot dedup).
- **fragments_2um** — `samtools view -q 30 -d CB | awk` → BED (chrom, start, start+readlen, CB, 1).
- **convert_bin** — remap 2 µm barcodes → 40 µm square bins via `barcode_mappings.parquet`.
- **qc_fragspot** — fragments per in-tissue 40 µm spot: stats + spatial map (`qc_fragspot.py`).

## Requirements (NYGC)

- `module load spaceranger/4.0.1`
- micromamba env **`spacetag_hd_analysis`**: snakemake, samtools (≥1.17), preseq, bgzip/tabix,
  python ≥3.10 + pandas + pyarrow, bc.
- cutadapt + fastqc from env **`spatial_multiome`** (used by the `trim` rule).
- matplotlib/pandas/pyarrow from env **`sctm`** (used by `qc_fragspot`).
- mm10 GEX reference: `/gpfs/commons/home/cyan/genome/refdata-gex-mm10-2020-A`.

## Configure

One `config.yaml` per sample (copy `config_template.yaml`):

```yaml
sample: "2604il_27ac"          # FASTQ sample prefix (spaceranger --sample)
spaceranger_id: "2604il_27ac"  # run id / output-dir stem
slide: "H1-Y3Z2ZNH"
area:  "A1"
# --- input FASTQs ---
fastq_folder: ".../2604il_27ac/fastqs/trimmed"   # SpaceRanger reads here (trim writes here)
# --- optional in-pipeline trimming ---
trim: true                     # false = FASTQs are already trimmed (skip cutadapt)
platform: "illumina"           # illumina | ultima  (selects the cutadapt adapters)
raw_fastq_folder: ".../2604il_27ac/fastqs/raw"   # raw lane FASTQs; only needed when trim: true
# --- images / reference ---
cytaimage: ".../..._A1_27ac.tif"
image: ".../...27ac-Spot000002.jpg"
transcriptome: "/gpfs/commons/home/cyan/genome/refdata-gex-mm10-2020-A"
segmentation_file: ""          # leave blank unless using a cell-seg mask
run_cells: false               # true also emits cell-level fragments
bin_size: 40
```

**Trimming / platform.** Set `trim: true` to trim inside the pipeline; `platform` picks the R2 3′
adapter set:
- `illumina` → `CTGTCTCTTATACACATCT`, `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`, poly-A
- `ultima`   → `AGATGTGTATAAGAGACAG`, `ACACTCTTTCCCTACACGACGCTCTTCCGATCT`, poly-T

Common cutadapt params: `--minimum-length 28:50 --nextseq-trim=20 -q 0,15`.
If your FASTQs are already trimmed, set `trim: false` (default) and point `fastq_folder` at them.
Lane-split FASTQs (e.g. `<sample>_S1_L001..L016_R{1,2}_001.fastq.gz`) are handled automatically.

## Run

```bash
sbatch --chdir=<sample_dir> /path/to/run_pipeline.sh <ABS/path/to/config.yaml>
```

`run_pipeline.sh` activates the `spacetag_hd_analysis` env and runs
`snakemake --cores 16 --rerun-incomplete` in one 16-core / 256 GB / 48 h SLURM allocation.
Dry-run first if you like:

```bash
micromamba run -n spacetag_hd_analysis \
  snakemake --snakefile Snakefile --configfile <config.yaml> --cores 16 -n
```

## Notes

- The pipeline does **not** branch on single- vs paired-end: Space Ranger always ingests
  R1 = barcode / R2 = genomic. Only R2 is genomic, so fragment/insert size is **not** measurable
  (NucleosomeSignal / FragmentHistogram N/A) regardless of instrument.
- `deduplicate` runs `fixmate -m` unconditionally (needed for PE BAMs; a no-op otherwise).
