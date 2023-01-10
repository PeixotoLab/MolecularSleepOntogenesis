#!/bin/sh

#SBATCH --partition=peixoto
#SBATCH --job-name=SalmonIndex
#SBATCH --error=SalmonIndex.err
#SBATCH --output=SalmonIndex.out
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1

module load salmon

salmon index -t /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/gentrome_transcripts_mouse.fa.gz \
             -d /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/decoys_mouse.txt \
             -i /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
             --gencode --threads 4 -k 31

