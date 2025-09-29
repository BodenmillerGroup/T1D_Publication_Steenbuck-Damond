#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=3850
#SBATCH --time=72:00:00
#SBATCH --output=/home/nsteen/scripts/T1D_analysis/workflow/log/main_%j

module load mamba
module load singularityce/4.1.4
eval "$(conda shell.bash hook)"
conda activate snakemake_cluster

snakemake --profile /home/nsteen/scripts/T1D_analysis/workflow/env/slurm "$@"