#!/bin/bash
#SBATCH --account=PDNU0016
#SBATCH --job-name=S.eliza_transcriptome_TF_family_annotation
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --output=tf_family_eli_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=duong_m1@denison.edu

python3 S.elizabethae_transcriptome_TF_family_annotated.py
