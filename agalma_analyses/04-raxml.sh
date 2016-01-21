#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=120G
#SBATCH -C e5-2600
#SBATCH -t 72:00:00
set -e
module load raxml/8.2.0
SEED=2395182752288
raxmlHPC-PTHREADS -f a -T 20 -m PROTGAMMAWAG -p $SEED -x $SEED -s alignment.phy -# 100 -n AnimalTreeML
