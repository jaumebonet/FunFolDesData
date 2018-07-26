#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8192
#SBATCH --time 24:00:00
#SBATCH --array=1-200
#SBATCH --output=output.%A_%a.out
#SBATCH --error=output.%A_%a.err


srun /work/upcorreia/bin/rosetta_future/devel/RosettaNubInitio/rosetta_scripts.default.linuxiccrelease -parser:protocol fold_design.xml -in:file:s 188_1kx8_0033.pdb -in:file:native 188_1kx8_0033.pdb -out:nstruct 65 -in:ignore_waters -in:ignore_unrecognized_res -out:file:silent_struct_type binary -out:prefix ${SLURM_ARRAY_TASK_ID}_1kx8_ -out:file:silent 1kx8_silent2_${SLURM_ARRAY_TASK_ID}_ -parser:script_vars frags=frag_2
echo "CASTOR: RUN FINISHED"

