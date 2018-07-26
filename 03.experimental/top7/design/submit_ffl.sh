#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8192
#SBATCH --time 24:00:00
#SBATCH --array=1-200
#SBATCH --output=output.%A_%a.out
#SBATCH --error=output.%A_%a.err


srun /work/upcorreia/bin/rosetta_future/devel/RosettaNubInitio/rosetta_scripts.default.linuxiccrelease -parser:protocol ffl.xml -in:file:s 1QYS.pdb  -out:nstruct 4 -in:ignore_waters -in:ignore_unrecognized_res -out:file:silent_struct_type binary -out:prefix ${SLURM_ARRAY_TASK_ID}_ffl_ -out:file:silent ffl_silent_out_${SLURM_ARRAY_TASK_ID}_
echo "CASTOR: RUN FINISHED"

