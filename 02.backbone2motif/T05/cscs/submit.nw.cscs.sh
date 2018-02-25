#!/bin/bash -l
#SBATCH --ntasks=288
#SBATCH --ntasks-per-node=36
#SBATCH --nodes=8
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task 1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --time 02:00:00
#SBATCH --job-name="T05"
#SBATCH --mem=120GB
#SBATCH --output=/scratch/snx3000/jbonet/logs/T05.%A_%a.out
#SBATCH --error=/scratch/snx3000/jbonet/logs/T05.%A_%a.err

export OMP_NUM_THREADS=1

# Define Rosetta executable (not accessible form queues)
ROSETTACOMP=mpistatic.linuxiccrelease
ROSETTAAPP=rosetta_scripts

# The actual executable
ROSETTAEXE=/scratch/snx3000/${USER}/bin/${ROSETTAAPP}.${ROSETTACOMP}

# Tracer Outputs
TRACEDIR=/scratch/snx3000/${USER}/tracers/

# Fixed number of structures: 572
NSTRUCT="-nstruct 572"

# Common
COMMON="-overwrite -in:ignore_unrecognized_res -in:ignore_waters -out:file:silent_struct_type binary -out:mute protocols.abinitio protocols.abinitio.foldconstraints protocols.moves core.optimization"

# Redirect Database
ROSETTADB="-in:path:database /scratch/snx3000/"${USER}"/database"

# JOB
EXP="nubinitio"
FRAG="wauto"
PREFIX=${EXP}_${FRAG}_${SLURM_TASK_PID}
OUTDIR=${EXP}_${FRAG}
PDB="2pkoA"

# Trace log files
TRACER="-mpi_tracer_to_file "${TRACEDIR}"/"${PREFIX}"_tracer.log"

mkdir -p ${OUTDIR}

srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ${ROSETTAEXE} -parser:protocol run.${EXP}.xml -s ${PDB}.pdb -out:prefix ${PREFIX}_ -out:file:silent ${OUTDIR}/${PREFIX} ${TRACER} ${NSTRUCT} ${ROSETTADB} ${COMMON} -parser:script_vars frags=${FRAG} pdb=${PDB}
