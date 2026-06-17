#!/bin/bash
#SBATCH --job-name=STagHD_pipeline
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=48:00:00

CONFIG="${1:?Usage: sbatch --chdir=<sample_dir> run_pipeline.sh <absolute/path/to/config.yaml>}"
WORKDIR=$(dirname "$(realpath "$CONFIG")")
SNAKEFILE=/gpfs/commons/home/cyan/data/7_space_tag_hd/pipeline/Snakefile

export MAMBA_EXE='/gpfs/commons/home/cyan/.local/bin/micromamba'
export MAMBA_ROOT_PREFIX='/gpfs/commons/home/cyan/miniconda3'
eval "$($MAMBA_EXE shell hook --shell bash)"
micromamba activate /gpfs/commons/home/cyan/miniconda3/envs/spacetag_hd_analysis
export PYTHONNOUSERSITE=1

echo "--- Environment Check ---"
echo "CONDA_PREFIX: $CONDA_PREFIX"
which snakemake
snakemake --version

cd "$WORKDIR"
mkdir -p logs

snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIG" \
    --cores 16 \
    --printshellcmds \
    --rerun-incomplete
