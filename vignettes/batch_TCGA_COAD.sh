#!/bin/sh
#SBATCH -A iicd
#SBATCH -J COAD
#SBATCH -t 5-00:00:00
#   SBATCH -t 12:00:00
#   SBATCH --cpus-per-task=32
#   SBATCH --mem=700gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=knd2127@columbia.edu

module load R

pwd

echo "Launching an R run"
date

R CMD BATCH --no-save --vanilla DECODE_for_TCGA_COAD.r routput
