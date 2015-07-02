#!/bin/bash -l
#SBATCH -D $HOME/AtCanary/scripts
#SBATCH -o $HOME/AtCanary/results/%j-out.txt
#SBATCH -e $HOME/AtCanary/results/%j-err.txt
#SBATCH -J dadi
#SBATCH -p idle_0064
set -e
set -u

# script usage:
# sbatch submit_model.sh <model.py>
model="${1}"

python "${model}"
