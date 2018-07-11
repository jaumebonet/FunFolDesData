#!/bin/bash -l
#SBATCH --ntasks=288
#SBATCH --ntasks-per-node=36
#SBATCH --nodes=8
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task 1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --time 03:00:00
#SBATCH --job-name="abinitio"
#SBATCH --output=/scratch/snx3000/jbonet/logs/abinitio.%A_%a.out
#SBATCH --error=/scratch/snx3000/jbonet/logs/abinitio.%A_%a.err

export OMP_NUM_THREADS=1

# Define Rosetta executable (not accessible form queues)
ROSETTACOMP=mpistatic.linuxiccrelease
ROSETTAAPP=AbinitioRelax

# The actual executable
ROSETTAEXE=/scratch/snx3000/${USER}/bin/${ROSETTAAPP}.${ROSETTACOMP}

# Tracer Outputs
TRACEDIR=/scratch/snx3000/${USER}/tracers/

# Fixed number of structures: 572
NSTRUCT="-nstruct 572"

# Common
COMMON="-abinitio:relax -relax:fast -abinitio:increase_cycles 10 -abinitio:rg_reweight 0.5 -abinitio:rsd_wt_helix 0.5 -abinitio:rsd_wt_loop 0.5"

# Redirect Database
ROSETTADB="-in:path:database /scratch/snx3000/"${USER}"/database"

# JOB
INPUTS="-in:file:fasta ../"${1}".fasta  -in:file:native ../"${1}".pdb"
FRAGS="-in:file:frag3 ../"${2}"/"${3}".200.3mers -in:file:frag9 ../"${2}"/"${3}".200.9mers"

# Trace log files
TRACER="-mpi_tracer_to_file "${TRACEDIR}"/"${1}"_tracer.log"

echo -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ${ROSETTAEXE} ${INPUTS} ${FRAGS} ${COMMON} ${NSTRUCT} ${ROSETTADB} ${TRACER} -out:prefix ${2}_${1} -out:file:silent abinitio_${2}.silent
srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ${ROSETTAEXE} ${INPUTS} ${FRAGS} ${COMMON} ${NSTRUCT} ${ROSETTADB} ${TRACER} -out:prefix ${2}_${1} -out:file:silent abinitio_${2}.silent 
