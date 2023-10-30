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

/home/cesar.medinaculma/NGSEP/NGSEPcore_4.3.2.jar

scp -r cesar.medinaculma@ceres.scinet.usda.gov:/home/cesar.medinaculma/NGSEP /home/cesar.medinaculma
scp -r cesar.medinaculma@ceres.scinet.usda.gov:/home/cesar.medinaculma/4_vcf_monoploid /home/cesar.medinaculma
scp -r medin297@mangi.msi.umn.edu:/home/samac/medin297/msi/6_DAl23-8031/py3_server .
 /project/xu_alfalfabreeding/system_from_home/MSI

scp -r py3_server cesar.medinaculma@ceres.scinet.usda.gov:/project/xu_alfalfabreeding/system_from_home/MSI

/home/cesar.medinaculma/NGSEP/
/project/xu_alfalfabreeding/system_from_home
/project/xu_alfalfabreeding/system_from_home/MSI
