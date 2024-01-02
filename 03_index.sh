#!/bin/bash
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8000MB
#SBATCH --time=0-03:00:00
#SBATCH --job-name=indexing_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/2_reference_genome/output_index%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/2_reference_genome/error_index%j.o

WORKDIR="/data/users/okopp/rnaseq_course"
REFDIR="$WORKDIR/2_reference_genome/reference"

module load UHTS/Aligner/hisat/2.2.1

hisat2-build -p 2 "$REFDIR/Mus_musculus_genome.fa" "$REFDIR/Mus_musculus"

