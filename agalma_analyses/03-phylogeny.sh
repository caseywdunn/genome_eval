#!/bin/sh
#SBATCH -J phylogeny
#SBATCH -t 3-00:00:00
#SBATCH -N 4
#SBATCH -c 16
#SBATCH --mem=60g
#SBATCH -C intel
#SBATCH --exclusive

set -e

export AGALMA_DB="/gpfs/data/cdunn/analyses/genomes2016_DEBpreproposal.sqlite"
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"
export BIOLITE_HOSTLIST=$(hostlist -e -s, $SLURM_NODELIST)
export BIOLITE_TOOLS="raxml-pthreads=raxmlHPC-PTHREADS"

IMPORT_IDS=$(agalma diagnostics runid -n import)
echo $IMPORT_IDS

ID=AnimalTree

agalma homologize --id $ID $IMPORT_IDS --genome_type any --molecule_type any
agalma multalign --id $ID
agalma genetree --id $ID
agalma treeprune --id $ID
agalma multalign --id $ID
agalma supermatrix --id $ID
agalma supermatrix --id $ID --proportion 0.5
agalma supermatrix --id $ID --proportion 0.75
