#!/bin/bash
#SBATCH --job-name=STagHD_pipeline
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
# 1. Path to micromamba binary
export MAMBA_EXE='/gpfs/commons/home/cyan/.local/bin/micromamba'

# 2. Point to the root directory where your envs are stored
export MAMBA_ROOT_PREFIX='/gpfs/commons/home/cyan/miniconda3'

# 3. Initialize the shell hook
eval "$($MAMBA_EXE shell hook --shell bash)"

# 4. Activate using the absolute path to the 'spacetag_hd_analysis' env
micromamba activate /gpfs/commons/home/cyan/miniconda3/envs/spacetag_hd_analysis

# 5. Debugging: confirm we are in the right place
echo "--- Environment Check ---"
echo "CONDA_PREFIX: $CONDA_PREFIX"
which snakemake
snakemake --version

# ==============================================================================
# EXECUTION
# ==============================================================================
# Ensure logs folder exists for snakemake's own logs
mkdir -p logs

# ==============================================================================
# EXECUTION
# ==============================================================================
# Remove --reason or --display-reasons to avoid the error you just saw
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cores 16 \
    --printshellcmds \
    --rerun-incomplete