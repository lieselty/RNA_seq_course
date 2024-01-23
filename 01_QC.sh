#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000MB
#SBATCH --time=0-02:00:00
#SBATCH --job-name=fastqc_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/1_quality_checks/output_%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/1_quality_checks/error_%j.o

## define variables
WORKDIR="/data/users/okopp/rnaseq_course"
READSDIR="$WORKDIR/reads/*.fastq.gz"
OUTDIR="$WORKDIR/1_quality_checks/fastqc"
mkdir -p $OUTDIR

## load module fastQC (version 0.11.9)
module add UHTS/Quality_control/fastqc/0.11.9

## command line
fastqc -o $OUTDIR -f $READSDIR
