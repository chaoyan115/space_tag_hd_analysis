import os

# ── Setup ──────────────────────────────────────────────────────────────────────
BIN_SIZE = config["bin_size"]
BIN_STR  = f"{int(BIN_SIZE):03d}um"
FULL_ID  = f"{config['spaceranger_id']}_{BIN_STR}"
SPACERANGER_OUTS = os.path.join(config["fastq_folder"], FULL_ID)

seg_path  = config.get("segmentation_file")
SEG_PARAMS = f"--custom-segmentation-file={seg_path} --nucleus-expansion-distance-micron=20" if seg_path else ""

# ── Optional in-pipeline trimming (Illumina now / Ultima option retained) ───────
# Backward-compatible: configs without these keys default to no-trim, so existing
# (already-trimmed) Ultima/Illumina sample dirs are unaffected.
PLATFORM         = config.get("platform", "illumina")
DO_TRIM          = bool(config.get("trim", False))
RAW_FASTQ_FOLDER = config.get("raw_fastq_folder", "")
TRIM_SENTINEL    = os.path.join(config["fastq_folder"], ".trim_complete")

# cutadapt 3' adapters are applied to R2 (genomic read) only; R1 is the spatial barcode.
# poly-A/T written as literal 9-mers (== cutadapt 'A{9}') to avoid any brace ambiguity in the shell.
ILLUMINA_ADAPTERS = "-A 'CTGTCTCTTATACACATCT;o=6' -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;o=6' -A 'AAAAAAAAA'"
ULTIMA_ADAPTERS   = "-A 'AGATGTGTATAAGAGACAG;o=6' -A 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT;o=6' -A 'TTTTTTTTT'"
ADAPTER_ARGS = ULTIMA_ADAPTERS if PLATFORM == "ultima" else ILLUMINA_ADAPTERS
CUTADAPT = "/gpfs/commons/home/cyan/miniconda3/envs/spatial_multiome/bin/cutadapt"
FASTQC   = "/gpfs/commons/home/cyan/miniconda3/envs/spatial_multiome/bin/fastqc"

# ── Targets ────────────────────────────────────────────────────────────────────
TARGETS = [
    os.path.join(SPACERANGER_OUTS, "outs", "possorted_genome_bam.bam"),
    os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated.bam"),
    os.path.join(SPACERANGER_OUTS, "markdup_stats.txt"),
    os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments2.bed.gz"),
    os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments2.bed.gz.tbi"),
    os.path.join(SPACERANGER_OUTS, f"possorted_genome_bam_deduplicated_fragments_{BIN_STR}.bed.gz"),
    os.path.join(SPACERANGER_OUTS, f"possorted_genome_bam_deduplicated_fragments_{BIN_STR}.bed.gz.tbi"),
    os.path.join(SPACERANGER_OUTS, "50M_complexity_curve.txt"),
    os.path.join(SPACERANGER_OUTS, "qc_fragspot_stats.txt"),
    os.path.join(SPACERANGER_OUTS, f"qc_fragspot_{BIN_STR}.pdf"),
]
if config.get("run_cells", False):
    TARGETS += [
        os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments_CELLS.bed.gz"),
        os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments_CELLS.bed.gz.tbi"),
    ]

rule all:
    input: TARGETS


# ── Rule 0: Trim (optional, in-pipeline) ───────────────────────────────────────
# Runs cutadapt on every raw lane FASTQ of this sample (handles lane-split input),
# writing trimmed FASTQs into fastq_folder with names preserved so SpaceRanger
# globs all lanes. Only wired in when config `trim: true`. Adapter set chosen by
# `platform` (illumina|ultima). Emits a sentinel so downstream rules can depend on it.
rule trim:
    output:
        sentinel = TRIM_SENTINEL
    threads: 16
    resources:
        mem_mb = 32000
    shell:
        r"""
        set -euo pipefail
        mkdir -p {config[fastq_folder]}
        cd {RAW_FASTQ_FOLDER}
        shopt -s nullglob
        for r1 in {config[sample]}_*_R1_001.fastq.gz; do
            [[ $r1 == Undetermined* ]] && continue
            r2="${{r1/_R1_001.fastq.gz/_R2_001.fastq.gz}}"
            base=$(basename "$r1" _R1_001.fastq.gz)
            echo "[trim] $base  (platform={PLATFORM})"
            {CUTADAPT} --cores={threads} \
                --minimum-length 28:50 \
                --nextseq-trim=20 \
                -q 0,15 \
                {ADAPTER_ARGS} \
                -o "{config[fastq_folder]}/${{base}}_R1_001.fastq.gz" \
                -p "{config[fastq_folder]}/${{base}}_R2_001.fastq.gz" \
                "$r1" "$r2" > "{config[fastq_folder]}/${{base}}_cutadapt.log" 2>&1
        done
        mkdir -p {config[fastq_folder]}/fastqc_trimmed
        {FASTQC} -t {threads} -o {config[fastq_folder]}/fastqc_trimmed {config[fastq_folder]}/*_R[12]_001.fastq.gz || true
        touch {output.sentinel}
        """


# ── Rule 1: SpaceRanger ────────────────────────────────────────────────────────
rule spaceranger:
    input:
        ([TRIM_SENTINEL] if DO_TRIM else [])
    output:
        bam     = os.path.join(SPACERANGER_OUTS, "outs", "possorted_genome_bam.bam"),
        bai     = os.path.join(SPACERANGER_OUTS, "outs", "possorted_genome_bam.bam.bai"),
        parquet = os.path.join(SPACERANGER_OUTS, "outs", "barcode_mappings.parquet")
    threads: 16
    shell:
        """
        source /etc/profile.d/modules.sh
        module load spaceranger/4.0.1

        # Restart guard: skip if already complete
        if [ -f {output.bam} ]; then
            echo "SpaceRanger output already exists — skipping."
            exit 0
        fi

        rm -rf {SPACERANGER_OUTS}

        spaceranger count --id={FULL_ID} \
            --transcriptome={config[transcriptome]} \
            --fastqs={config[fastq_folder]} \
            --sample={config[sample]} \
            --image={config[image]} \
            --cytaimage={config[cytaimage]} \
            --slide={config[slide]} \
            --area={config[area]} \
            --custom-bin-size={BIN_SIZE} \
            --output-dir={SPACERANGER_OUTS} \
            --create-bam=true \
            --localcores={threads} \
            --localmem=114 \
            {SEG_PARAMS}
        """


# ── Rule 2: Complexity curve ───────────────────────────────────────────────────
rule complexity:
    input:
        bam = rules.spaceranger.output.bam,
        bai = rules.spaceranger.output.bai
    output:
        subsampled_bam = temp(os.path.join(SPACERANGER_OUTS, "subsampled_50M.bam")),
        sorted_bam     = temp(os.path.join(SPACERANGER_OUTS, "subsampled_50M.sorted.bam")),
        complexity     = os.path.join(SPACERANGER_OUTS, "50M_complexity_curve.txt")
    threads: 16
    resources:
        mem_mb = 64000
    shell:
        """

        TOTAL=$(samtools idxstats {input.bam} | awk '{{s+=$3+$4}} END {{print s}}')
        FRACTION=$(echo "scale=6; 50000000 / $TOTAL" | bc)

        samtools view -@ {threads} -b -s 42$FRACTION {input.bam} > {output.subsampled_bam}
        samtools sort  -@ {threads} {output.subsampled_bam} -o {output.sorted_bam}
        preseq lc_extrap -B -o {output.complexity} {output.sorted_bam}
        """


# ── Rule 3: Deduplicate ────────────────────────────────────────────────────────
# Uses samtools markdup ≥1.17 --barcode-tag CB for per-spot deduplication.
# Reads at the same coordinate but different CB barcodes are NOT flagged as
# duplicates — equivalent to umi_tools dedup --per-cell but ~10x faster and
# without the chr-splitting loop.
#
# Paired-end BAMs require fixmate -m (adds mate-score tag) before markdup.
# Pipeline: name-sort → fixmate → coord-sort → markdup --barcode-tag CB
rule deduplicate:
    input:  rules.spaceranger.output.bam
    output:
        bam   = os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated.bam"),
        stats = os.path.join(SPACERANGER_OUTS, "markdup_stats.txt")
    threads: 16
    resources:
        mem_mb = 200000
    shell:
        """
        TMPDIR=$(dirname {output.bam})
        TMP_NSORT=$TMPDIR/tmp_namesorted.bam
        TMP_FIX=$TMPDIR/tmp_fixmate.bam
        TMP_CSORT=$TMPDIR/tmp_coordsorted.bam

        echo "[dedup] Step 1: name-sort"
        samtools sort -n -m 1G -@ {threads} {input} -o $TMP_NSORT

        echo "[dedup] Step 2: fixmate"
        samtools fixmate -m -@ {threads} $TMP_NSORT $TMP_FIX
        rm -f $TMP_NSORT

        echo "[dedup] Step 3: coord-sort"
        samtools sort -m 1G -@ {threads} $TMP_FIX -o $TMP_CSORT
        rm -f $TMP_FIX

        echo "[dedup] Step 4: markdup"
        samtools markdup -r --barcode-tag CB -f {output.stats} -@ {threads} $TMP_CSORT {output.bam}
        rm -f $TMP_CSORT

        echo "[dedup] Indexing"
        samtools index -@ {threads} {output.bam}
        """


# ── Rule 4: Generate 2 µm base fragments ──────────────────────────────────────
# Key improvements vs original:
#   -d CB  : pre-filters to reads with CB tag (avoids awk guard check)
#   match(): single regex on the full line is faster than looping fields 12..NF
#   bgzip -@ : multi-threaded compression
rule fragments_2um:
    input:  rules.deduplicate.output.bam
    output:
        bed = os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments2.bed.gz"),
        tbi = os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments2.bed.gz.tbi")
    threads: 8
    resources:
        mem_mb = 40000
    shell:
        """

        samtools view -@ {threads} -q 30 -d CB {input} | \
        awk 'BEGIN{{OFS="\\t"}} {{
            if (match($0, /CB:Z:[^\\t]+/)) {{
                bc = substr($0, RSTART+5, RLENGTH-5)
                print $3, $4, $4+length($10), bc, "1"
            }}
        }}' | \
        sort -k1,1 -k2,2n -k3,3n -S 30G --parallel={threads} | \
        bgzip -@ {threads} > {output.bed}

        tabix -p bed {output.bed}
        """


# ── Rule 5: Bin-level fragments ────────────────────────────────────────────────
rule convert_bin:
    input:
        frag    = rules.fragments_2um.output.bed,
        mapping = rules.spaceranger.output.parquet
    output:
        frag_bin = os.path.join(SPACERANGER_OUTS, f"possorted_genome_bam_deduplicated_fragments_{BIN_STR}.bed.gz"),
        tbi      = os.path.join(SPACERANGER_OUTS, f"possorted_genome_bam_deduplicated_fragments_{BIN_STR}.bed.gz.tbi")
    params:
        bin_col = f"square_{BIN_STR}"
    threads: 5
    resources:
        mem_mb = 32000
    shell:
        """
python3 - <<'EOF' | sort -k1,1 -k2,2n -S 50G --parallel={threads} | bgzip -@ {threads} > {output.frag_bin}
import pandas as pd
import gzip, sys, io

mappings = pd.read_parquet("{input.mapping}", columns=["square_002um", "{params.bin_col}"])
mappings = mappings.dropna(subset=["{params.bin_col}"])
translate = dict(zip(mappings["square_002um"], mappings["{params.bin_col}"]))

out = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", write_through=True)
with gzip.open("{input.frag}", "rt") as f:
    for line in f:
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) >= 4:
            new_bc = translate.get(parts[3])
            if new_bc:
                parts[3] = new_bc
                out.write("\\t".join(parts) + "\\n")
EOF
        tabix -p bed {output.frag_bin}
        """


# ── Rule 6: Cell-level fragments ───────────────────────────────────────────────
rule convert_to_cells:
    input:
        frag    = rules.fragments_2um.output.bed,
        mapping = rules.spaceranger.output.parquet
    output:
        frag_cells = os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments_CELLS.bed.gz"),
        tbi        = os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments_CELLS.bed.gz.tbi")
    threads: 5
    resources:
        mem_mb = 32000
    shell:
        """
python3 - <<'EOF' | sort -k1,1 -k2,2n -S 50G --parallel={threads} | bgzip -@ {threads} > {output.frag_cells}
import pandas as pd
import gzip, sys, io

mappings = pd.read_parquet("{input.mapping}", columns=["square_002um", "cell_id", "in_cell"])
mappings = mappings[mappings["in_cell"] == True]
translate = dict(zip(mappings["square_002um"], mappings["cell_id"]))

out = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", write_through=True)
with gzip.open("{input.frag}", "rt") as f:
    for line in f:
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) >= 4:
            new_bc = translate.get(parts[3])
            if new_bc:
                parts[3] = new_bc
                out.write("\\t".join(parts) + "\\n")
EOF
        tabix -p bed {output.frag_cells}
        """


# ── Rule 7: Fragments-per-spot QC (40 µm) ──────────────────────────────────────
# Reads the binned fragments + the SpaceRanger tissue_positions.parquet, keeps
# in_tissue spots, and writes summary stats + a spatial map coloured by fragments
# per capturing spot. Depends on convert_bin so the SpaceRanger outs (parquet) exist.
rule qc_fragspot:
    input:
        frag = rules.convert_bin.output.frag_bin
    output:
        stats = os.path.join(SPACERANGER_OUTS, "qc_fragspot_stats.txt"),
        pdf   = os.path.join(SPACERANGER_OUTS, f"qc_fragspot_{BIN_STR}.pdf")
    params:
        positions = os.path.join(SPACERANGER_OUTS, "outs", "binned_outputs", f"square_{BIN_STR}", "spatial", "tissue_positions.parquet"),
        outprefix = os.path.join(SPACERANGER_OUTS, f"qc_fragspot_{BIN_STR}"),
        sample    = config["spaceranger_id"]
    threads: 4
    resources:
        mem_mb = 32000
    shell:
        r"""
        export MAMBA_ROOT_PREFIX=/gpfs/commons/home/cyan/miniconda3
        /gpfs/commons/home/cyan/.local/bin/micromamba run -n sctm python \
            /gpfs/commons/home/cyan/data/7_space_tag_hd/pipeline/qc_fragspot.py \
            --frag {input.frag} \
            --positions {params.positions} \
            --outprefix {params.outprefix} \
            --stats {output.stats} \
            --sample {params.sample}
        """
