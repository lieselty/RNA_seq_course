#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000MB
#SBATCH --time=0-01:00:00

## change directory
cd ../reads

## fastqc
fastqc -o ../1_quality_checks/fastqc -f fastq ../reads/*fastq.gz
