#!/bin/bash
#SBATCH --job-name=output_pairwise
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16000M
#SBATCH --time=218:00:00
#SBATCH --output=output_pairwise.out
#SBATCH --mail-type=END

python3 ../scripts/pairwise_run.py woc-185-tsplib-b8dd6b77 6 21
