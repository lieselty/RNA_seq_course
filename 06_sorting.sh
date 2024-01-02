#!/bin/bash
#SBATCH --array=1-16
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=25G
#SBATCH --time=0-01:00:00
#SBATCH --job-name=sorting_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/2_reference_genome/output/output_sorting%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/2_reference_genome/output/error_sorting%j.o

# define variables
WORKDIR="/data/users/okopp/rnaseq_course"
BAMDIR="$WORKDIR/2_reference_genome/mapping/bam"
OUTDIR="$WORKDIR/2_reference_genome/mapping/sorted"
SAMPLELIST="$WORKDIR/2_reference_genome/mapping/samplelist.tsv"
mkdir -p $OUTDIR
#load module
module load UHTS/Analysis/samtools/1.10

SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)

samtools sort -@ 4 -m 24G -o $OUTDIR/${SAMPLE}sorted.bam -T temp_bam $BAMDIR/${SAMPLE}_mapping.bam

