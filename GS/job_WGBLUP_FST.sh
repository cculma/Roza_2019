#!/bin/bash
#SBATCH --partition=cahnrs		# Partition (like a queue in PBS)
#SBATCH --job-name=nph			# Job Name
#SBATCH --output=nph.out
#SBATCH --error=nph.err
#SBATCH --time=2-00:00:00		# Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=3			# Node count required for the job
#SBATCH --ntasks-per-node=1		# Number of tasks to be launched per Node
#SBATCH --cpus-per-task=20		# Number of threads per task (OMP threads)
#SBATCH --mem=220GB			# Amount of memory per node

module load r
Rscript --vanilla WGBLUP_FST.R