#!/bin/bash -l
#SBATCH -D /home/lv70590/Arun/AtCanary/scripts
#SBATCH -o /home/lv70590/Arun/AtCanary/results/%j-out.txt
#SBATCH -e /home/lv70590/Arun/AtCanary/results/%j-err.txt
#SBATCH -J dadi
#SBATCH -p mem_0064
set -e
set -u

# script usage:
# sbatch submit_model.sh <model.py>
model="${1}"

python "${model}"
