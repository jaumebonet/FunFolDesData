#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8192
#SBATCH --time 24:00:00
#SBATCH --array 1-50
#SBATCH --job-name="sd"
#SBATCH --output=./output.%A_%a.out
#SBATCH --error=./output.%A_%a.err


srun /work/upcorreia/bin/rosetta_future/stable/RosettaWeekly/rosetta_scripts.default.linuxiccrelease -restore_talaris_behavior -s 87_ffl_twostrand_1qys_Lys_0001.pdb -in:file:native 87_ffl_twostrand_1qys_Lys_0001.pdb -parser:protocol design_part.xml -ignore_unrecognized_res -ignore_zero_occupancy  -use_input_sc -flip_HNQ -no_optH false -out:prefix sd_top7_${SLURM_ARRAY_TASK_ID}  -out:file:silent sd_top7_nterm2_ -nstruct 2
echo "CASTOR: RUN FINISHED"
