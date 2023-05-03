#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=128g
#SBATCH -t 96:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simmo536@umn.edu
#SBATCH -p amdsmall
#SBATCH -o %A.out
#SBATCH -e %A.err

cd ~/2021Q2-PTOmodeling
module load matlab
matlab -nodisplay -r \
"parpool('local',$SLURM_JOB_CPUS_PER_NODE); \
iiSS = $SSlb:$SSub; \
study_coulombPTO_dampingStudy_multiSSmultiSeed"

# Useful commands:
# sbatch --export=SSlb=1,SSub=19 study_coulombPTO_dampingStudy_multiSSmultiSeed.sh
# sbatch --export=SSlb=20,SSub=38 study_coulombPTO_dampingStudy_multiSSmultiSeed.sh
# sbatch --export=SSlb=39,SSub=57 study_coulombPTO_dampingStudy_multiSSmultiSeed.sh
# sbatch --export=SSlb=58,SSub=76 study_coulombPTO_dampingStudy_multiSSmultiSeed.sh
# sbatch --export=SSlb=77,SSub=95 study_coulombPTO_dampingStudy_multiSSmultiSeed.sh
# sbatch --export=SSlb=96,SSub=114 study_coulombPTO_dampingStudy_multiSSmultiSeed.sh
# dos2unix  study_coulombPTO_dampingStudy_multiSSmultiSeed.sh

