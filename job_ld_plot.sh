#!/bin/bash
#SBATCH --partition=kamiak		# Partition (like a queue in PBS)
#SBATCH --job-name=gwas			# Job Name
#SBATCH --output=gwas.out
#SBATCH --error=gwas.err
#SBATCH --time=0-06:00:00		# Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=2			# Node count required for the job
#SBATCH --ntasks-per-node=1		# Number of tasks to be launched per Node
#SBATCH --cpus-per-task=20		# Number of threads per task (OMP threads)
#SBATCH --mem=220GB			# Amount of memory per node

module load r
Rscript --vanilla ld_plot.R
