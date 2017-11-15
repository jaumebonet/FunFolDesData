#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8000
#SBATCH --job-name="pfam"
#SBATCH --time 10:00:00
#SBATCH --array=0-500
#SBATCH --output=/scratch/bonet/logs/pfam.%A_%a.out
#SBATCH --error=/scratch/bonet/logs/pfam.%A_%a.err

mkdir -p out
srun /scratch/lpdi/bin/Rosetta/devel/nubinitio/main/source/bin/rosetta_scripts.linuxiccrelease -parser:protocol pf00586.xml -s 3fd6B.pdb -in:ignore_unrecognized_res -nstruct 10 -out:file:silent out/pf00586_${SLURM_ARRAY_TASK_ID}

echo "CASTOR: RUN FINISHED"
