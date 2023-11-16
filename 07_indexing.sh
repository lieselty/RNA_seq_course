#!/bin/bash
#SBATCH --array=1-16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=0-01:00:00
#SBATCH --job-name=indexing_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/2_reference_genome/output/output_indexing%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/2_reference_genome/output/error_indexing%j.o

# define variables
WORKDIR="/data/users/okopp/rnaseq_course"
SORTDIR="$WORKDIR/2_reference_genome/mapping/sorted"
OUTDIR="$WORKDIR/2_reference_genome/mapping/index"
SAMPLELIST="$WORKDIR/2_reference_genome/mapping/samplelist.tsv"

#load module
module load UHTS/Analysis/samtools/1.10

SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)

samtools index $SORTDIR/${SAMPLE}sorted.bam