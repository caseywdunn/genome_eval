#!/bin/sh
#SBATCH -J speciestree
#SBATCH -t 24:00:00
#SBATCH -c 16
#SBATCH --mem=60g
#SBATCH -C intel
#SBATCH --exclusive

set -e

module load examl/3.0.14

export AGALMA_DB="/gpfs/data/cdunn/analyses/genomes2016_DEBpreproposal.sqlite"
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"
export BIOLITE_TOOLS="raxml-pthreads=raxmlHPC-PTHREADS,examl=examl-OMP-AVX"
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

ID=ExpressionTree

agalma speciestree --id $ID --bootstrap 100
