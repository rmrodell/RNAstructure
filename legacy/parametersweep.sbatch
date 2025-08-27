#!/usr/bin/bash
#SBATCH --job-name=sweep_20250731_%j
#SBATCH --output=sweep_%j.out
#SBATCH --error=sweep_%j.err
#SBATCH --time=24:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=64
#SBATCH --mem=64GB
#SBATCH --mail-user=rodell@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL

ml R
ml biology viennarna

Rscript /scratch/users/rodell/motifmatcher/RNAstructure/parametersweep_v5.R -i "/scratch/users/rodell/motifmatcher/input_pool1.csv" -f "/scratch/users/rodell/RNAfold_psipos" -o "/scratch/users/rodell/motifmatcher/sweep_20250731" -d "/scratch/users/rodell/motifmatcher/pool1_psipos_v2.csv"
