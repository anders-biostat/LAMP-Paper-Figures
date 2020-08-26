**LAMP-sequencing raw read processing**

This directory contains all scripts required to process LAMP-sequencing raw reads (example file 'LAMP-sequencing_raw_sample100k.fastq.gz') with the result of a R data object containing sample-wise virus-matching and -unmatching counts.

In order to run the workflow on the sample file above, run:
```bash
bash 00-run_workflow.sh
```

The system used to run this workflow had the following required software installed:
- julia (version 1.3.1) with the following packages installed
  - ArgParse (version 1.1.0)
  - CodecZlib (version 0.6.0)
  - BioSequences (version 2.0.1)
- cutadapt (version 2.8)
- R (version 4.0.2) with the following packages installed
  - readr (version 1.3.1)
  - dplyr (version 1.0.2)
  - stringr (version 1.4.0)
