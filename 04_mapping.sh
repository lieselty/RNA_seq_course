#!/bin/bash
#SBATCH --array=1-16
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000MB
#SBATCH --time=0-02:00:00
#SBATCH --job-name=mapping_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/2_reference_genome/output_mapping%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/2_reference_genome/error_mapping%j.o

# define variables
WORKDIR="/data/users/okopp/rnaseq_course"
OUTDIR="$WORKDIR/2_reference_genome/mapping"
SAMPLELIST="$WORKDIR/2_reference_genome/mapping/samplelist.tsv"
INDEX="/data/users/okopp/rnaseq_course/2_reference_genome/reference/Mus_musculus"
READS="/data/users/okopp/rnaseq_course/reads"

#load module
module load UHTS/Aligner/hisat/2.2.1

SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)
READ1="$READS/${SAMPLE}_1.fastq.gz"
READ2="$READS/${SAMPLE}_2.fastq.gz"


hisat2 -x $INDEX -1 $READ1 -2 $READ2 -S "$OUTDIR/${SAMPLE}_mapping.sam" -p 4
