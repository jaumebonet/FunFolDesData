#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=serial
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4096
#SBATCH --time 03:00:00
#SBATCH --array=1-50
#SBATCH --output=log/abi.%A_%a.out
#SBATCH --error=log/abi.%A_%a.err


NSTRUCT="-nstruct 1"

# JOB
INPUTS="-in:file:s ../top7_full.pdb -in:file:native ../top7_full.pdb"

mkdir -p  relkoys/
srun /work/upcorreia/bin/rosetta_bin/stable/RosettaWeekly/rosetta_scripts.linuxiccrelease -parser:protocol relax.xml ${INPUTS} ${NSTRUCT} -out:file:scorefile relkoys.sc -out:path:pdb relkoys/ -out:suffix ${SLURM_ARRAY_TASK_ID}  -overwrite
