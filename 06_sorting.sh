#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=25G
#SBATCH --time=0-01:00:00
#SBATCH --job-name=sambam_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/2_reference_genome/output/output_sorting%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/2_reference_genome/output/error_sorting%j.o

# define variables
WORKDIR="/data/users/okopp/rnaseq_course"
BAMDIR="$WORKDIR/2_reference_genome/mapping/bam"
OUTDIR="$WORKDIR/2_reference_genome/mapping"

#load module
module load UHTS/Analysis/samtools/1.10

samtools sort -m 24G -o $OUTDIR/sorted.bam -T temp -@ 4 $BAMDIR/*.bam


