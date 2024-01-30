# DE Toxoplasma - rna-seq course
In this repository, you will find all the code to reproduce my analysis based on the reads from Singhania *et al.* (2019).
## Steps
### Qualtiy check
First use fastQC to assess the quality of the reads: [01_QC](01_QC.sh) to create the fastqc reports, and then use MultiQC to have all your fastqc in one report: [01_1MultiQC](01_1MultiQC.sh).
There is also a code to have the number of reads: [02_nbreads](02_nbreads.sh) it's not needed since we have this information in the fastqc files.
### Creating the index
First, you will need to download the reference genome from ensembl. Then run [03_index](03_index.sh).

### Mapping the reads against the reference genome
To do that, run [04_mapping](04_mapping.sh), [05_sam_to_bam](05_sam_to_bam.sh), [06_sorting](06_sorting.sh) and [07_indexing](07_indexing.sh)
Then you will have a .bam and a .bam.bai file for each of your reads.

### Counting the number of reads per genes
You will need to download the annotation corresponding to the reference genome you download before. Then run [08_reads_per_genes](08_reads_per_genes.sh) and [09_count_table_for_DESeq2](09_count_table_for_DESeq2.sh) to get the count table that you will use for the next step. Download it locally.

### Exploratory data analysis, Differential expression analysis & Overrepresentation analysis
To run the differential expression analysis, you will need another data table looking like that:
| Sample      | Group             |
|-------------|-------------------|
| SRR7821918  | Lung_WT_Case      |
| SRR7821919  | Lung_WT_Case      |
|     ...     |        ...        |


Then you can run the R scrip [10_DE](10_DE.R) and you will get everything.

## On this repository, you will also find:
### - All the [fastqc](fastqc) and [multiqc](fastqc/multiqc_report.html) reports 
### - A [recap from the quality of the mapping](mapping_hisat2_output.xlsx)
### - The [summary of featureCounts](count_table.txt.summary)
