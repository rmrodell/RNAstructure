#!/usr/bin/bash
#SBATCH --job-name=20250718_motif_%j
#SBATCH --output=20250718_motif_%j.out
#SBATCH --error=20250718_motif_%j.err
#SBATCH --time=12:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=64
#SBATCH --mem=32GB
#SBATCH --mail-user=rodell@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL

ml R

Rscript parametersweep_v2.R -o parametersweep_paralleled_output
