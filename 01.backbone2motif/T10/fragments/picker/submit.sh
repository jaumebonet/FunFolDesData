#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16384
#SBATCH --time 2:30:00
#SBATCH --output=output.%A_%a.out
#SBATCH --error=output.%A_%a.err

srun /work/upcorreia/bin/rosetta_future/stable/tools/fragment_tools/make_fragments.pl 2qdfA.fa
