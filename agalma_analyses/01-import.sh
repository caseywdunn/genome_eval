#!/bin/sh
#SBATCH -t 60

set -e

export AGALMA_DB="/gpfs/data/cdunn/analyses/genomes2016_DEBpreproposal.sqlite"

agalma import --id "Monosiga" --seq_type aa
agalma import --id "Salpinfoeca" --seq_type aa
agalma import --id "Capsaspora" --seq_type aa
agalma import --id "Sphaeroforma" --seq_type aa
agalma import --id "Amphimedon" --seq_type aa
agalma import --id "Mnemiopsis" --seq_type aa
agalma import --id "Trichoplax" --seq_type aa
agalma import --id "Nematostella" --seq_type aa
agalma import --id "Acropora" --seq_type aa
agalma import --id "Strongylocentrotus" --seq_type aa
agalma import --id "Ptychodera" --seq_type aa
agalma import --id "Saccoglossus" --seq_type aa
agalma import --id "Ciona_int" --seq_type aa
agalma import --id "Branchiostoma" --seq_type aa
agalma import --id "Petromyzon" --seq_type aa
agalma import --id "Takifugu" --seq_type aa
agalma import --id "Xenopus" --seq_type aa
agalma import --id "Anolis" --seq_type aa
agalma import --id "Taeniopygia" --seq_type aa
agalma import --id "Canis" --seq_type aa
agalma import --id "Tribolium" --seq_type aa
agalma import --id "Daphnia" --seq_type aa
agalma import --id "Lottia" --seq_type aa
agalma import --id "Octopus" --seq_type aa
agalma import --id "Lingula" --seq_type aa
agalma import --id "Schmidtea" --seq_type aa
agalma import --id "Pristionchus" --seq_type aa
agalma import --id "Hypsibius" --seq_type aa
agalma import --id "Capitella" --seq_type aa
agalma import --id "Tetranychus" --seq_type aa
agalma import --id "Globodera" --seq_type aa
