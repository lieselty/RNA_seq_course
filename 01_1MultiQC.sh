#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000MB
#SBATCH --time=0-02:00:00
#SBATCH --job-name=multiqc_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/1_quality_checks/output_multiqc%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/1_quality_checks/error_multiqc%j.o

## load module MultiQC (version 1.8)
module load UHTS/Analysis/MultiQC/1.8

## define variables
WORKDIR="/data/users/okopp/rnaseq_course"
FASTQCDIR="$WORKDIR/1_quality_checks/fastqc"

## command line
multiqc $FASTQCDIR/*_fastqc.zip