#!/bin/bash
#SBATCH --partition=msismall		# Partition (like a queue in PBS)
#SBATCH --job-name=svm			# Job Name
#SBATCH --output=svm.out
#SBATCH --error=svm.err
#SBATCH --time=2-00:00:00		# Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=3			# Node count required for the job
#SBATCH --ntasks-per-node=1		# Number of tasks to be launched per Node
#SBATCH --cpus-per-task=20		# Number of threads per task (OMP threads)
#SBATCH --mem=220GB			# Amount of memory per node

module load r
Rscript --vanilla SVM1.R
