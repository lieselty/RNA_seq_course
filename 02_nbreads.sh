#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000MB
#SBATCH --time=0-00:30:00
#SBATCH --job-name=reads_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/1_quality_checks/output_reads_%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/1_quality_checks/error_reads_%j.o


## Define variables
READSDIR="/data/users/okopp/rnaseq_course/reads"
OUTPUTDIR="/data/users/okopp/rnaseq_course/1_quality_checks/nb_reads.txt"


## Loop through all .fastq files in the directory
for file in "$READSDIR"/*.fastq.gz; do
    # Count the number of 4-line groups in the file
    read_count=$(zcat "$file" | wc -l)
    read_count=$((read_count / 4))

    # Append the result to the output file
    echo "$file has $read_count reads" >> "$OUTPUTDIR"
done

echo "Results saved to $OUTPUTDIR"

#this code is not really neede because you get the info in the fastqc files