#!/bin/bash -l
#SBATCH --ntasks=288
#SBATCH --ntasks-per-node=36
#SBATCH --nodes=8
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task 1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --time 05:00:00
#SBATCH --job-name="bcl2_classicFFL"
#SBATCH --output=/scratch/snx3000/jbonet/logs/bcl2_classicFFL.%A_%a.out
#SBATCH --error=/scratch/snx3000/jbonet/logs/bcl2_classicFFL.%A_%a.err

export OMP_NUM_THREADS=1

# Define Rosetta executable (not accessible form queues)
ROSETTABIN=/project/s731/bin/Rosetta
ROSETTACOMP=mpistatic.linuxiccrelease
ROSETTAVERSION=nubinitio/main/source/bin
ROSETTAAPP=rosetta_scripts
ROSETTAEXE=${ROSETTABIN}/${ROSETTAVERSION}/${ROSETTAAPP}.${ROSETTACOMP}

# The actual executable
ROSETTAEXE=/scratch/snx3000/jbonet/bin/${ROSETTAAPP}.${ROSETTACOMP}

# Tracer Outputs
TRACEDIR=/scratch/snx3000/jbonet/tracers/

# Fixed number of structures: 572
NSTRUCT="-nstruct 572"

# Redirect Database
ROSETTADB="-in:path:database /scratch/snx3000/"${USER}"/database"

# Arguments
FFL=${1}
DRC=${2}
ARGUMENTS="-parser:script_vars ffl="${FFL}" desrel="${DRC}

# Naming conventions
JOBNAME="4oydMimic"
JOBTYPE="classicFFL"
PREFIX=${JOBTYPE}_${FFL}_${DRC}_${SLURM_TASK_PID}
SCRIPT="-parser:protocol "${JOBNAME}.${JOBTYPE}.xml
OUTDIR="../output/"${FFL}"/"${DRC}

# Trace log files
TRACER="-mpi_tracer_to_file "${TRACEDIR}"/"${PREFIX}"_tracer_"${SLURM_TASK_PID}".log"

mkdir -p ${OUTDIR}

srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ${ROSETTAEXE} ${SCRIPT} @common_flags -out:prefix ${PREFIX}_ -out:file:silent ${OUTDIR}/${PREFIX} ${TRACER} ${NSTRUCT} ${ROSETTADB} ${ARGUMENTS}
