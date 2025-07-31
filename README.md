# donertas-group/trimscreen


## Introduction

**donertas-group/trimscreen** is a bioinformatics pipeline that systematically evaluate the influence of different trimming strategies on 16S amplicon sequencing data. This pipeline screens all possible combinations of forward and reverse read trimming lengths, ranging from no trimming to aggressive trimming that still ensures a minimum overlap between paired-end reads. By processing each combination through a standard bioinformatics workflow, the pipeline aims to identify the trimming lengths that maximize observed taxonomic richness. This approach provides a data-driven method to optimize preprocessing parameters and improve the accuracy and resolution of microbial community profiling.

## Usage
First, prepare a samplesheet with the *full path* of your input data file that looks as follows. There should only be `A`-`Z`, `0`-`9` and `_` in sample names.

`samplesheet.csv` (the examples below assumes you have two samples `S10A`, `S10B` and a control `C01`):

```csv
sampleID,forwardReads,reverseReads
S10A,/<full_path_to>/E10A_R1.fastq.gz,/<full_path_to>/E10A_R2.fastq.gz
S10B,/<full_path_to>/E10B_R1.fastq.gz,/<full_path_to>/E10B_R2.fastq.gz
C01,/<full_path_to>/control_R1.fastq.gz,/<full_path_to>/control_R2.fastq.gz
```
Each row represents a pair of fastq files (paired end). Single-end is not enabled for this pipeline. Do not change the header line.

Then prepare a metadata sheet that looks as follows:

`metadata.csv`:

```csv
ID,condition
S10A,sample
S10B,sample
C01,control
```

Lastly prepare a parameter file:
`params.yaml`

```yaml
enable_filter: true
classification_bbduk: true
validation_blastn: true
fasta_blastn: '/scratch/shire/data/nj/reference/genome/nothobranchius_furzeri/NfurGRZ-RIMD1/GCF_043380555.1_NfurGRZ-RIMD1_genomic.fna.gz'
fasta_bbduk: '/scratch/shire/data/nj/reference/genome/nothobranchius_furzeri/NfurGRZ-RIMD1/GCF_043380555.1_NfurGRZ-RIMD1_genomic.fna.gz'

read_length: 300
FW_primer: 'CAATGGRSGVRASYCTGAHS'
RV_primer: 'AGGGTATCTAATCCT'
marker_size_min: 300
minimum_overlap: 50
step_size: 5

dada_ref_taxonomy: "silva=138"
min_len_asv: 300
max_len_asv: 500
```
`enable_filter`: whether host sequence removal should be performed. If set to `true`:
`classification_bbduk`: whether bbduk should be run.
`validation_blastn`: whether blastn should be run.
`fasta_bbduk` and `fasta_blastn` are the reference genome of the host.

Parameters required for setting trimming length screening:
`read_length`: number of cycles, i.e. length of raw sequencing read.
`marker_size_min`: minimum size of the expected marker amplicon. E.g. for 16s V4, expected marker size is 280 to 400 bp.
`step_size`: Difference between two adjacent trimming lengths to be screened, can be set to 1 to 10.
`minimum_overlap`: minimum overlapping length (bp) of paried reads, default 50, should not be smaller than 20.


Now, you can run the pipeline using:

```bash
mkdir <YOUR_DIR>; cd <YOUR_DIR>
git clone https://github.com/donertas-group/trimscreen.git
nextflow run trimscreen \
    -profile test,apptainer \
    --outdir <OUTDIR> \
    --input <INPUT_DIR>/samplesheet.csv \
    --run_host_removal true \
    -params-file <INPUT_DIR>/params.yaml \
    --metadata <INPUT_DIR>/metadata.tsv
```















