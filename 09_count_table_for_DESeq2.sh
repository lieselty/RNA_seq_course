#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=0-00:05:00
#SBATCH --job-name=table_okopp

## define variables
WORKDIR="/data/users/okopp/rnaseq_course"
COUNTTABLE="$WORKDIR/3_count_reads_per_gene/count_table.txt"
OUTDIR="$WORKDIR/3_count_reads_per_gene"
FINAL_COUNT_TABLE="final_count_table.txt"
mkdir -p $OUTDIR

## define an array of sample names and corresponding file paths
declare -a SAMPLES=(
    "SRR7821918"
    "SRR7821919"
    "SRR7821920"
    "SRR7821921"
    "SRR7821922"
    "SRR7821937"
    "SRR7821938"
    "SRR7821939"
    "SRR7821949"
    "SRR7821950"
    "SRR7821951"
    "SRR7821952"
    "SRR7821953"
    "SRR7821968"
    "SRR7821969"
    "SRR7821970"
)

## use awk to extract relevant columns and perform substitutions
awk -F'\t' 'NR > 1 {
    printf "%s\t", $1; 
    for (i = 7; i <= NF; i++) {
        printf "%s", $i; 
        if (i < NF) printf "\t";
    }
    printf "\n";
}' "$COUNTTABLE" > "$OUTDIR/$FINAL_COUNT_TABLE"

## loop through the array and perform substitutions
for SAMPLE in "${SAMPLES[@]}"; do
    sed -i "s|/data/users/okopp/rnaseq_course/2_reference_genome/mapping/sorted/${SAMPLE}sorted.bam|${SAMPLE}|" "$OUTDIR/$FINAL_COUNT_TABLE"
done
