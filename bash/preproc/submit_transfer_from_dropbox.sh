#!/bin/bash
#SBATCH --job-name=transfer_data
#SBATCH --time=48:00:00
#SBATCH -n 5
#SBATCH --mail-user=sagon151@stanford.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=normal,gse

module load system rclone 

rclone copy EdNeuroDropBox:2025_EAR_Groupitizing_EEG/Data/ /scratch/users/sagon151/eeg-dot-digit-data/raw_data/ --include "*.mat"
