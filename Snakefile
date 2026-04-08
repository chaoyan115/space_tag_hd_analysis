import os

configfile: "config.yaml"

# ... (your existing ID and Path setup) ...
# 1. Setup Dynamic Strings
BIN_SIZE = config["bin_size"]
BIN_STR = f"{int(BIN_SIZE):03d}um" # Results in '024um'

# 2. Combine for the specific ID you requested
# Results in e.g., "27ac_06x_24um"
FULL_ID = f"{config['spaceranger_id']}_{BIN_STR}"

# 3. Define the exact path where 'outs' will live
# /gpfs/.../27ac_0.6x/27ac_06x_24um
SPACERANGER_OUTS = os.path.join(config["fastq_folder"], FULL_ID)

# Check if segmentation_file is present and not empty
seg_path = config.get("segmentation_file")

if seg_path:
    # This string will be inserted into the shell command
    SEG_PARAMS = f"--custom-segmentation-file={seg_path} --nucleus-expansion-distance-micron=20"
else:
    # This leaves the command clean if no file is provided
    SEG_PARAMS = ""

rule all:
    input:
        # 1. Basic SpaceRanger output
        os.path.join(SPACERANGER_OUTS, "outs", "possorted_genome_bam.bam"),
        
        # 2. Complexity/Quality metrics
        os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_complexity.txt"),
        
        # 3. The Cleaned/Deduplicated BAM
        os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated.bam"),
        
        # 4. The 2um Base Fragments (Base level)
        os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments2.bed.gz"),
        
        # 5. The Binned Fragments (The file you'll use in Signac/Seurat)
        os.path.join(SPACERANGER_OUTS, f"possorted_genome_bam_deduplicated_fragments_{BIN_STR}.bed.gz"),
        
        os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments_CELLS.bed.gz")

        os.path.join(SPACERANGER_OUTS, "50M_complexity_curve.txt")

rule spaceranger:
    output:
        bam = os.path.join(SPACERANGER_OUTS, "outs", "possorted_genome_bam.bam"),
        parquet = os.path.join(SPACERANGER_OUTS, "outs", "barcode_mappings.parquet")
    threads: 16
    shell:
        """
        module load spaceranger/4.0.1
        
        # Cleanup to avoid 'pipestance' errors
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


rule complexity:
    input:
        bam="outs/possorted_genome_bam.bam",
        bai="outs/possorted_genome_bam.bam.bai"
    output:
        subsampled_bam=temp(os.path.join(SPACERANGER_OUTS, "subsampled_50M.bam")),
        sorted_bam=temp(os.path.join(SPACERANGER_OUTS, "subsampled_50M.sorted.bam")),
        complexity=os.path.join(SPACERANGER_OUTS, "50M_complexity_curve.txt")
    threads: 16
    resources:
        mem_mb=64000
    shell:
        """
        # 1. Get total count instantly from index
        TOTAL=$(samtools idxstats {input.bam} | awk '{{s+=$3+$4}} END {{print s}}')
        
        # 2. Calculate fraction for 50M (using bc for float math)
        # We use 50,000,000 / TOTAL
        FRACTION=$(echo "scale=6; 50000000 / $TOTAL" | bc)
        
        # 3. Subsample (Seed 42)
        samtools view -@ {threads} -b -s 42$FRACTION {input.bam} > {output.subsampled_bam}
        
        # 4. Sort (Required for Preseq)
        samtools sort -@ {threads} {output.subsampled_bam} -o {output.sorted_bam}
        
        # 5. Run Preseq
        preseq lc_extrap -B -o {output.complexity} {output.sorted_bam}
        """

rule deduplicate:
    input: rules.spaceranger.output.bam
    output: os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated.bam")
    threads: 10
    shell:
        """
        TEMP_DIR="{SPACERANGER_OUTS}/tmp_split_dedup"
        mkdir -p $TEMP_DIR
        CHROMS=$(samtools idxstats {input} | cut -f1 | grep -E "^chr|^[0-9XYM]" | grep -v '*')

        for CHROM in $CHROMS; do
            (
            samtools view -b -h -e " [CB] && [UB] " {input} $CHROM > $TEMP_DIR/${{CHROM}}_filtered.bam
            if [[ ! -s $TEMP_DIR/${{CHROM}}_filtered.bam ]]; then exit 0; fi
            samtools index $TEMP_DIR/${{CHROM}}_filtered.bam
            umi_tools dedup -I $TEMP_DIR/${{CHROM}}_filtered.bam --extract-umi-method=tag \
                --per-cell --cell-tag=CB --umi-tag=UB --ignore-umi --method=unique \
                --stdout=$TEMP_DIR/${{CHROM}}_dedup.bam
            ) &
            if [[ $(jobs -r | wc -l) -ge 10 ]]; then wait -n; fi
        done
        wait
        samtools merge -f -@ {threads} {output} $TEMP_DIR/*_dedup.bam
        samtools index -@ {threads} {output}
        rm -rf $TEMP_DIR
        """

rule fragments_2um:
    input: rules.deduplicate.output
    output: os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments2.bed.gz")
    threads: 8
    shell:
        """
        samtools view -@ {threads} -q 30 {input} | \
        awk 'BEGIN{{OFS="\\t"}} {{
            bc=""
            for(i=12; i<=NF; i++) {{
                if($i ~ /^CB:Z:/) {{ bc=substr($i,6); break }}
            }}
            if(bc != "") print $3, $4, $4+length($10), bc, "1"
        }}' | \
        sort -k1,1 -k2,2n -k3,3n -S 30G --parallel={threads} | \
        bgzip -c > {output}
        tabix -p bed {output}
        """

rule convert_bin:
    input:
        frag = rules.fragments_2um.output,
        mapping = rules.spaceranger.output.parquet
    output:
        frag_bin = os.path.join(SPACERANGER_OUTS, f"possorted_genome_bam_deduplicated_fragments_{BIN_STR}.bed.gz"),
        tbi = os.path.join(SPACERANGER_OUTS, f"possorted_genome_bam_deduplicated_fragments_{BIN_STR}.bed.gz.tbi")
    params:
        bin_col = f"square_{BIN_STR}"
    shell:
        """
python3 - <<EOF | sort -k1,1 -k2,2n -S 50G --parallel=5 | bgzip -c > {output.frag_bin}
import pandas as pd
import gzip
import sys
try:
    mappings = pd.read_parquet("{input.mapping}", columns=["square_002um", "{params.bin_col}"])
    translate = dict(zip(mappings["square_002um"], mappings["{params.bin_col}"]))
except Exception as e:
    sys.stderr.write(f"Error: {{e}}\\n")
    sys.exit(1)
with gzip.open("{input.frag}", 'rt') as f:
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) >= 4 and parts[3] in translate:
            parts[3] = translate[parts[3]]
            sys.stdout.write('\\t'.join(parts) + '\\n')
EOF
        tabix -p bed {output.frag_bin}
        """



rule convert_to_cells:
    input:
        frag = rules.fragments_2um.output,
        mapping = rules.spaceranger.output.parquet
    output:
        frag_cells = os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments_CELLS.bed.gz"),
        tbi = os.path.join(SPACERANGER_OUTS, "possorted_genome_bam_deduplicated_fragments_CELLS.bed.gz.tbi")
    shell:
        """
python3 - <<EOF | sort -k1,1 -k2,2n -S 50G --parallel=5 | bgzip -c > {output.frag_cells}
import pandas as pd
import gzip
import sys

try:
    # Load mapping but only keep bins that are assigned to a cell
    mappings = pd.read_parquet("{input.mapping}", columns=["square_002um", "cell_id", "in_cell"])
    mappings = mappings[mappings["in_cell"] == True]
    
    # Create the lookup dictionary: {{ 'square_002um_ID' : 'cellid_N' }}
    translate = dict(zip(mappings["square_002um"], mappings["cell_id"]))
except Exception as e:
    sys.stderr.write(f"Error loading parquet: {{e}}\\n")
    sys.exit(1)

with gzip.open("{input.frag}", 'rt') as f:
    for line in f:
        parts = line.strip().split('\\t')
        # parts[3] is the 2um barcode from the BAM
        if len(parts) >= 4 and parts[3] in translate:
            parts[3] = translate[parts[3]] # Replace bin ID with Cell ID
            sys.stdout.write('\\t'.join(parts) + '\\n')
EOF
        tabix -p bed {output.frag_cells}
        """