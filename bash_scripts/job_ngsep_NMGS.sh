#!/bin/bash
#SBATCH --partition=kamiak	# Partition (like a queue in PBS)
#SBATCH --job-name=ngsep_%j	# Job Name
#SBATCH --output=ngsep_%j.out
#SBATCH --error=ngsep_%j.err
#SBATCH --time=5-0:00:00		# Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=2				# Node count required for the job
#SBATCH --ntasks-per-node=1		# Number of tasks to be launched per Node
#SBATCH --cpus-per-task=60		# Number of threads per task (OMP threads)

module load java
module load picard
#./ngsep_NMGS.sh
./runMapping_NGSEP.sh
