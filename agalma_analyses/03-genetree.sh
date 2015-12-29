#!/bin/sh
#SBATCH -J genetree
#SBATCH -t 7-00:00:00
#SBATCH -N 8
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

ID=ExpressionTree

agalma homologize --id $ID $IMPORT_IDS --genome_type any --molecule_type any
agalma multalign --id $ID
agalma genetree --id $ID
