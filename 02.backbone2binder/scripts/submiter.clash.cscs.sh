#!/bin/bash -l
#SBATCH --ntasks=288
#SBATCH --ntasks-per-node=36
#SBATCH --nodes=8
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task 1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --time 02:00:00
#SBATCH --job-name="bcl2_clash"
#SBATCH --output=/scratch/snx3000/jbonet/logs/bcl2_clash.%A_%a.out
#SBATCH --error=/scratch/snx3000/jbonet/logs/bcl2_clash.%A_%a.err

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
NSTRUCT="-nstruct 1"

# Redirect Database
ROSETTADB="-in:path:database /scratch/snx3000/"${USER}"/database"

# Arguments
INDIR=${1}
INFILE=${2}

# Naming conventions
JOBNAME="binder"
JOBTYPE="evaluate"
PREFIX=${JOBNAME}_${JOBTYPE}_${SLURM_TASK_PID}
SCRIPT="-parser:protocol "${JOBNAME}_${JOBTYPE}.xml
OUTDIR=${INDIR}

# Trace log files
TRACER="-mpi_tracer_to_file "${TRACEDIR}"/"${PREFIX}"_tracer_"${SLURM_TASK_PID}".log"

mkdir -p ${OUTDIR}

srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ${ROSETTAEXE} ${SCRIPT} -in:file:silent ${INFILE} -keep_input_scores false -out:prefix ${PREFIX}_ -out:file:silent ${OUTDIR}/${PREFIX} ${TRACER} ${NSTRUCT} ${ROSETTADB}

# Naming conventions
JOBNAME="binder"
JOBTYPE="clash"
PREFIX=${JOBNAME}_${JOBTYPE}_${SLURM_TASK_PID}
SCRIPT="-parser:protocol "${JOBNAME}_${JOBTYPE}.xml
OUTDIR=${INDIR}

# Trace log files
TRACER="-mpi_tracer_to_file "${TRACEDIR}"/"${PREFIX}"_tracer_"${SLURM_TASK_PID}".log"

mkdir -p ${OUTDIR}

srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ${ROSETTAEXE} ${SCRIPT} -in:file:silent ${INFILE} -keep_input_scores false -out:prefix ${PREFIX}_ -out:file:score_only ${OUTDIR}/${PREFIX} ${TRACER} ${NSTRUCT} ${ROSETTADB}

# Naming conventions
JOBNAME="minimize"
JOBTYPE="evaluate"
PREFIX=${JOBNAME}_${JOBTYPE}_${SLURM_TASK_PID}
SCRIPT="-parser:protocol "${JOBNAME}_${JOBTYPE}.xml
OUTDIR=${INDIR}

# Trace log files
TRACER="-mpi_tracer_to_file "${TRACEDIR}"/"${PREFIX}"_tracer_"${SLURM_TASK_PID}".log"

mkdir -p ${OUTDIR}

srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ${ROSETTAEXE} ${SCRIPT} -in:file:silent ${INFILE} -keep_input_scores false -out:prefix ${PREFIX}_ -out:file:score_only ${OUTDIR}/${PREFIX} ${TRACER} ${NSTRUCT} ${ROSETTADB}
