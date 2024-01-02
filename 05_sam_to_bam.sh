#!/bin/bash
#SBATCH --array=1-16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000MB
#SBATCH --time=0-01:00:00
#SBATCH --job-name=sambam_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/2_reference_genome/output/output_sambam%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/2_reference_genome/output/error_sambam%j.o

# define variables
WORKDIR="/data/users/okopp/rnaseq_course"
OUTDIR="$WORKDIR/2_reference_genome/mapping/bam"
SAMDIR="$WORKDIR/2_reference_genome/mapping/sam"
SAMPLELIST="$WORKDIR/2_reference_genome/mapping/samplelist.tsv"
mkdir -p $OUTDIR

#load module
module load UHTS/Analysis/samtools/1.10

SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)

samtools view -hbS $SAMDIR/${SAMPLE}_mapping.sam  -o $OUTDIR/${SAMPLE}_mapping.bam
