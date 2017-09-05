# Applied Human Computational Genomics Pipeline
**Chris Monaco \
  MS Bioinformatics 2018 | Georgia Institute of Technology**

Variant calling pipeline for genomic data analysis

## Class Server

Server: `gpuvannberg.biology.gatech.edu`
Data Folder: `/data2/AHCG2017FALL/`

## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```
# August 31 - Sept 5: Assessing a Standard Variant Call Pipeline

Based on data from this [paper](https://www.nature.com/nbt/journal/v31/n11/full/nbt.2696.html).

## Data Acquisition

## Bowtie Index

Pre-built index was supplied in the `../reference_genome` folder.

## Running Pipeline

  python ahcg_pipeline_v1.0.1.py -t /data2/AHCG2017FALL/bin/Trimmomatic-0.36 -b /data2/AHCG2017FALL/bin/bowtie2-2.2.9 -p /data2/AHCG2017FALL/bin/picard -g /data2/AHCG2017FALL/bin/GenomeAnalysisTK-3.8-0-ge9d806836 -i /data2/AHCG2017FALL/data/ -w /data2/AHCG2017FALL/reference_genome -r /data2/AHCG2017FALL/reference_genome -a /data2/AHCG2017FALL/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o /data2/AHCG2017FALL/output -d /data2/AHCG2017FALL/reference_genome/GATKResourceBundle


