#!/bin/bash -l
#SBATCH --ntasks=288
#SBATCH --ntasks-per-node=36
#SBATCH --nodes=8
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task 1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --time 05:00:00
#SBATCH --job-name="1kx8"
#SBATCH --output=/scratch/snx3000/jbonet/logs/1kx8_dgg.%A.out
#SBATCH --error=/scratch/snx3000/jbonet/logs/1kx8_ddg.%A.err

export OMP_NUM_THREADS=1

# Define Rosetta executable (not accessible form queues)
ROSETTACOMP=mpistatic.linuxiccrelease
ROSETTAAPP=rosetta_scriptsW

# The actual executable
ROSETTAEXE=/scratch/snx3000/${USER}/bin/${ROSETTAAPP}.${ROSETTACOMP}

# Tracer Outputs
TRACEDIR=/scratch/snx3000/${USER}/tracers/

# Common
COMMON="-overwrite -in:ignore_unrecognized_res -in:ignore_waters -out:file:silent_struct_type binary -ignore_zero_occupancy false -out:mute protocols.abinitio protocols.abinitio.foldconstraints protocols.moves core.optimization -restore_talaris_behavior"

# Redirect Database
ROSETTADB="-in:path:database /scratch/snx3000/"${USER}"/databaseW"

# JOB
OUTDIR="output/"${1}"_"${2}"/"
PREFIX=${3}
INSILENT=${OUTDIR}"/"${PREFIX}.silent.gz


# Trace log files
TRACER="-mpi_tracer_to_file "${TRACEDIR}"/"${PREFIX}"_tracer.log"

mkdir -p ${OUTDIR}

srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ${ROSETTAEXE} -parser:protocol ddg.xml -in:file:silent ${INSILENT} -out:file:silent ${OUTDIR}/${PREFIX}_ddg.silent ${TRACER} ${ROSETTADB} ${COMMON} -parser:script_vars protocol=${1}

