#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000MB
#SBATCH --time=0-02:00:00
#SBATCH --job-name=fastqc_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/1_quality_checks/output_%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/1_quality_checks/error_%j.o

##load module 
module add UHTS/Quality_control/fastqc/0.11.9

## fastqc
fastqc -o ../1_quality_checks/fastqc -f fastq ../reads/*.fastq.gz
