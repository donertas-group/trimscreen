# donertas-group/trimscreen


## Introduction

**donertas-group/trimscreen** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):
-->
First, prepare a samplesheet with the *full path* of your input data file that looks as follows:

`samplesheet.csv`:

```csv
sample,short_reads_fastq_1,short_reads_fastq_2,long_reads_fastq_1
E10A,/<full_path_to>/E10A_R1.fastq.gz,/<full_path_to>/E10A_R2.fastq.gz,
E10B,/<full_path_to>/E10B_R1.fastq.gz,/<full_path_to>/E10B_R2.fastq.gz,
CONTROL_REP1,/<full_path_to>/control_R1.fastq.gz,/<full_path_to>/control_R2.fastq.gz,
```
Each row represents a pair of fastq files (paired end). Single-end is not enabled for this pipeline. Don't forget the last comma at end of each line.

Then prepare a metadata sheet that looks as follows:

`metadata.tsv`:

```tsv
ID  condition
E10A    sample
E10B    sample
CONTROL_REP1    control
```
Note that metadata needs to be tab-separated file.

Lastly prepare a parameter file:
`params.yaml`

```yaml
fasta_bbduk: '/scratch/shire/data/nj/reference/genome/nothobranchius_furzeri/NfurGRZ-RIMD1/GCF_043380555.1_NfurGRZ-RIMD1_genomic.fna'
classification_bbduk: true
validation_blastn: true
fasta_blastn: '/scratch/shire/data/nj/reference/genome/nothobranchius_furzeri/NfurGRZ-RIMD1/GCF_043380555.1_NfurGRZ-RIMD1_genomic.fna'
enable_filter: true
FW_primer: 'CAATGGRSGVRASYCTGAHS'
RV_primer: 'AGGGTATCTAATCCT'
trunclenf_range: "235:5:240"
trunclenr_range: "200:5:205"
dada_ref_taxonomy: "silva=138"
min_len_asv: 200
max_len_asv: 500
```

the `fasta_bbduk` and `fasta_blastn` are the reference genome of the host.

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

















> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

donertas-group/trimscreen was originally written by Yi Wang.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use donertas-group/trimscreen for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
[![GitHub Actions CI Status](https://github.com/donertas-group/trimscreen/actions/workflows/ci.yml/badge.svg)](https://github.com/donertas-group/trimscreen/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/donertas-group/trimscreen/actions/workflows/linting.yml/badge.svg)](https://github.com/donertas-group/trimscreen/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/donertas-group/trimscreen)
