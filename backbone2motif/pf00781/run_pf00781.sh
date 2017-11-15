#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8000
#SBATCH --job-name="pfam"
#SBATCH --time 15:00:00
#SBATCH --array=0-1000
#SBATCH --output=/scratch/bonet/logs/pfam.%A_%a.out
#SBATCH --error=/scratch/bonet/logs/pfam.%A_%a.err

mkdir -p out
MODE=${1}
srun /scratch/lpdi/bin/Rosetta/devel/nubinitio/main/source/bin/rosetta_scripts.linuxiccrelease -parser:protocol pf00781.xml -s 3s40A.pdb -in:ignore_unrecognized_res -nstruct 10 -out:file:silent out/pf00781${MODE}_${SLURM_ARRAY_TASK_ID} -parser:script_vars mode=${MODE}



echo "CASTOR: RUN FINISHED"
