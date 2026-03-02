#!/bin/bash
#SBATCH --job-name=output_pairwise
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=218:00:00
#SBATCH --output=output_pairwise.out
#SBATCH --mail-type=END

python3 ../scripts/pairwise_linear_run.py random_instance_699_320_20 20 4
