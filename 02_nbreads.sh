#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500MB
#SBATCH --time=0-00:30:00
#SBATCH --job-name=reads_okopp
#SBATCH --output=/data/users/okopp/rnaseq_course/1_quality_checks/output_reads_%j.o
#SBATCH --error=/data/users/okopp/rnaseq_course/1_quality_checks/error_reads_%j.o

directory="/data/users/okopp/rnaseq_course/reads"

# Define the output file
output_file="/data/users/okopp/rnaseq_course/1_quality_checks/nb_reads.txt"


# Loop through all .fastq files in the directory
for file in "$directory"/*.fastq.gz; do
    # Count the number of 4-line groups in the file
    read_count=$(zcat "$file" | wc -l)
    read_count=$((read_count / 4))

    # Append the result to the output file
    echo "$file has $read_count reads" >> "$output_file"
done

echo "Results saved to $output_file"

#this code is not really neede because you get the info in the fastqc files