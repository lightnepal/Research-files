#!/bin/bash
#SBATCH --array=28-54
#SBATCH --mem=6G
source ~/.bashrc
srun ./octave1.sh


