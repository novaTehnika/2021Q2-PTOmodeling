#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=128g
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simmo536@umn.edu
#SBATCH -p amdsmall
#SBATCH -o %A.out
#SBATCH -e %A.err

cd ~/2021Q2-PTOmodeling
module load matlab
matlab -nodisplay -r \
"parpool('local',$SLURM_JOB_CPUS_PER_NODE); \
study_coulombPTO_dampingStudy_multiSS"

# sbatch  study_coulombPTO_dampingStudy_multiSS.sh
# dos2unix  study_coulombPTO_dampingStudy_multiSS.sh

