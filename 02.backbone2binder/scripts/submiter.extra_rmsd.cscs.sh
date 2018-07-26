#!/bin/bash -l
#SBATCH --ntasks=288
#SBATCH --ntasks-per-node=36
#SBATCH --nodes=8
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task 1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --time 01:30:00
#SBATCH --job-name="bcl2_rmsd"
#SBATCH --output=/scratch/snx3000/jbonet/logs/bcl2_rmsd.%A_%a.out
#SBATCH --error=/scratch/snx3000/jbonet/logs/bcl2_rmsd.%A_%a.err

export OMP_NUM_THREADS=1

# Define Rosetta executable (not accessible form queues)
ROSETTACOMP=mpistatic.linuxiccrelease
ROSETTAAPP=rosetta_scripts

# The actual executable
ROSETTAEXE=/scratch/snx3000/${USER}/bin/${ROSETTAAPP}.${ROSETTACOMP}

# Tracer Outputs
TRACEDIR=/scratch/snx3000/${USER}/tracers/

# Fixed number of structures: 572
NSTRUCT="-nstruct 1"

# Common
COMMON="-overwrite -in:ignore_unrecognized_res -in:ignore_waters -keep_input_scores false"

# Redirect Database
ROSETTADB="-in:path:database /scratch/snx3000/"${USER}"/database"

# Arguments
INDIR=${1}
INFILE=${2}

# Naming conventions
JOBNAME="extra"
JOBTYPE="rmsd"
PREFIX=${JOBNAME}_${JOBTYPE}_${SLURM_TASK_PID}
SCRIPT="-parser:protocol "${JOBNAME}_${JOBTYPE}.xml
OUTDIR=${INDIR}

# Trace log files
TRACER="-mpi_tracer_to_file "${TRACEDIR}"/"${PREFIX}"_tracer_"${SLURM_TASK_PID}".log"

mkdir -p ${OUTDIR}

srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ${ROSETTAEXE} ${SCRIPT} -in:file:silent ${INFILE}  -out:prefix ${PREFIX}_ -out:file:score_only ${OUTDIR}/${PREFIX} ${TRACER} ${NSTRUCT} ${ROSETTADB} ${COMMON}
