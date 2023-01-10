#!/bin/sh

#SBATCH --partition=peixoto
#SBATCH --job-name=SalmonAlign_P16_1
#SBATCH --error=SalmonAlign_P16_1.err
#SBATCH --output=SalmonAlign_P16_1.out
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1


module load salmon 

salmon quant \
--threads 6 \
--index /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P16/P16_HC_15_R1.fastq.gz --mates2 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P16/P16_HC_15_R2.fastq.gz \
--output /data/peixoto/SleepDevelopmentStudy/SalmonResults/March2022/WTHCP16_15_quant \
--numBootstraps 30
--validateMappings &

salmon quant \
--threads 6 \
--index /data/peixoto/SleepDevelopmentStudy/SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P16/P16_HC_16_R1.fastq.gz --mates2 /data/peixoto/SleepDevelopmentStudy/FASTQ/WT/P30/P16/P16_HC_16_R2.fastq.gz \
--output /data/peixoto/SleepDevelopmentStudy/SalmonResults/March2022/WTHCP16_16_quant \
--numBootstraps 30
--validateMappings &



wait

