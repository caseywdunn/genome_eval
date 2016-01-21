#!/bin/sh
#SBATCH -J annotate
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH -C intel
#SBATCH --exclusive
#SBATCH --array=1-31

set -e

export AGALMA_DB="/gpfs/data/cdunn/analyses/genomes2016_DEBpreproposal.sqlite"
export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"

IDS=(
Monosiga
Salpinfoeca
Capsaspora
Sphaeroforma
Amphimedon
Mnemiopsis
Trichoplax
Nematostella
Acropora
Strongylocentrotus
Ptychodera
Saccoglossus
Ciona_int
Branchiostoma
Petromyzon
Takifugu
Xenopus
Anolis
Taeniopygia
Canis
Tribolium
Daphnia
Lottia
Octopus
Lingula
Schmidtea
Pristionchus
Hypsibius
Capitella
Tetranychus
Globodera
)

ID=${IDS[$SLURM_ARRAY_TASK_ID-1]}
echo $ID

agalma annotate --id $ID
