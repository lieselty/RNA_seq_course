#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000
#SBATCH --time=0-02:00:00
#SBATCH --job-name=readsgenes_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/3_count_reads_per_gene/output_count%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/3_count_reads_per_gene/error_count%j.e

module load UHTS/Analysis/subread/2.0.1

WORKDIR="/data/users/okopp/rnaseq_course"
BAMDIR="$WORKDIR/2_reference_genome/mapping/sorted"
OUTDIR="$WORKDIR/3_count_reads_per_gene"
ANNOTDIR="$WORKDIR/2_reference_genome/reference"
mkdir -p $OUTDIR

featureCounts -p -a $ANNOTDIR/Mus_musculus.GRCm39.110.gtf -o $OUTDIR/count_table.txt -s 2 -T 4 $BAMDIR/*.bam
#try without -s