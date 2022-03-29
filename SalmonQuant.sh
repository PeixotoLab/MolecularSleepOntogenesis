#!/bin/sh

#SBATCH --partition=peixoto
#SBATCH --job-name=SalmonAlign_P90_2
#SBATCH --error=SalmonAlign_P90_2.err
#SBATCH --output=SalmonAlign_P90_2.out
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1


module load salmon 


salmon quant \
--threads 6 \
--index /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/WTSD5_PFC_1_R1.fastq.gz --mates2 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P90/WTSD5_PFC_1_R2.fastq.gz \
--output /data/peixoto/SleepDevelopmentStudy/SalmonResults/March2022/WTSDP90_1_quant \
--numBootstraps 30
--validateMappings &

salmon quant \
--threads 6 \
--index /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P90/WTSD5_PFC_2_R1.fastq.gz --mates2 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P90/WTSD5_PFC_2_R2.fastq.gz \
--output /data/peixoto/SleepDevelopmentStudy/SalmonResults/March2022/WTSDP90_2_quant \
--numBootstraps 30
--validateMappings &

salmon quant \
--threads 6 \
--index /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P90/WTSD5_PFC_3_R1.fastq.gz --mates2 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P90/WTSD5_PFC_3_R2.fastq.gz \
--output /data/peixoto/SleepDevelopmentStudy/SalmonResults/March2022/WTSDP90_3_quant \
--numBootstraps 30
--validateMappings &

salmon quant \
--threads 6 \
--index /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P90/WTSD5_PFC_4_R1.fastq.gz --mates2 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P90/WTSD5_PFC_4_R2.fastq.gz \
--output /data/peixoto/SleepDevelopmentStudy/SalmonResults/March2022/WTSDP90_4_quant \
--numBootstraps 30
--validateMappings &


salmon quant \
--threads 6 \
--index /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P90/WTSD5_PFC_5_R1.fastq.gz --mates2 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P90/WTSD5_PFC_5_R2.fastq.gz \
--output /data/peixoto/SleepDevelopmentStudy/SalmonResults/March2022/WTSDP90_5_quant \
--numBootstraps 30
--validateMappings &


wait

