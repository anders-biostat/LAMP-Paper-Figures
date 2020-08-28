## LAMP-sequencing to count workflow
####################################
## Exemplary workflow to process raw LAMP-sequencing reads to yield a count table of per-sample virus matching counts
## Author: Konrad Herbst (k.herbst@zmbh.uni-heidelberg.de)
## Version: 2020/08/26

## Step 01: assign wells to reads (based on offline i5 and i7 barcodes)
## Note: the used script also deals with offline UMIs which comes after the i5 barcode
julia assign_coordinates.jl --n_max_errors 1 --distance seqlev --outdir Coord_assignment \
      indices-P7-plates.csv indices-P5-wells.csv LAMP-sequencing_raw_sample100k.fastq.gz 2>&1 | tee 01-run_assign.log

## Step 02: run cutadapt
cutadapt --cores 7 --nextseq-trim 20 --minimum-length 43 --adapter 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC' \
         --output Coord_assignment/LAMP-sequencing_raw_sample100k_pass_trimmed.fastq.gz Coord_assignment/LAMP-sequencing_raw_sample100k_pass.fastq.gz \
         > 02-run_cutadapt.log

## Step 03: count reads
## Note: virus-specific sequences can be changed in the script
## Note: The script also extracts most prominent kmers from reads of individual wells which can be used for contamination analysis
julia count.jl Coord_assignment/LAMP-sequencing_raw_sample100k_pass_trimmed.fastq.gz 2>&1 | tee 03-run_counting.log

## Step 04: process count table
Rscript process_counts.R Counted/LAMP-sequencing_raw_sample100k_pass_trimmed_counted.tsv
