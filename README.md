# Applied Human Computational Genomics Pipeline
**Chris Monaco \
  MS Bioinformatics 2018 | Georgia Institute of Technology**

Variant calling pipeline for genomic data analysis

## Mission Statement

To create an efficient and easy to use pipeline for the early detection of cancer through non-invasive liquid biopsy from next generation sequencing data.

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

- Reference genome: Human Reference GRCh38 from [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)
- SRR948994 acquired with SRA toolkit 

## Bowtie Index

Pre-built index was supplied in the `../reference_genome` folder.

## Running Pipeline

Python ahcg_pipeline_v1.0.1Cai.py \
-t /data2/AHCG2017FALL/bin/Trimmomatic-0.36/trimmomatic-0.36.jar \
-b /data2/AHCG2017FALL/bin/bowtie2-2.2.9/bowtie2 \
-p /data2/AHCG2017FALL/bin/picard/picard.jar \
-g /data2/AHCG2017FALL/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-i /data2/AHCG2017FALL/data/SRR948994_1.fastq /data2/AHCG2017FALL/data/SRR948994_2.fastq \
-w /data2/AHCG2017FALL/reference_genome/Bowtie2Index/genome \
-r /data2/AHCG2017FALL/reference_genome/genome.fa \
-a /data2/AHCG2017FALL/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa \
-o /data2/AHCG2017FALL/output \
-d /data2/AHCG2017FALL/reference_genome/GATKResourceBundle/dbsnp_146.hg38.vcf.gz

## Current Version 1.0.3 Notes

- Specified Java 1.8 for Picard tools

# Virtual Machine

Prebuilt VMs were stored in `/data2/VMbox_prebuilt`. 

  Monaco, Christopher                    cmonaco3         Ubuntu-64-DR-AHCG2017-p10025.ova    10025
  
## Adding and accessing VM

VM was imported with this command
 
    $ vboxmanage import Ubuntu-64-DR-AHCG2017-p10025.ova

To get list of VMs to find correct ID

    $ vboxmanage list vms

VM was started with

    $ vboxmanage startvm Ubuntu-64-DR-AHCG2017 --type headless

VM was logged into with

    $ ssh vannberglab@localhost -p 10025
    password: vanberglab
  
WM can be shutdown with

    $ vboxmanage controlvm Ubuntu-64-DR-AHCG2017 poweroff soft

## Copying files from GPUVannberg to VM

On GPUVannberg:

    $ scp -r -P 10025 /data2/AHCG2017FALL/bin/ vannebrglab@localhost:~/
  
On VM:
  
    $ scp -r cmonaco3@gpuvannberg.biology.gatech.edu:/data2/AHCG2017FALL/data/ .

## Cloning and Increasing Disk Size

Cloning disk

    $ vboxmanage clonehd Ubuntu-64-DR-AHCG2017-p10025-disk001.vmdk Ubuntu-64-DR-AHCG2017.vdi --format vdi
    
Resizing disk to 120 gig

    $  vboxmanage modifyhd Ubuntu-64-DR-AHCG2017.vdi --resize 120000


# Read Depth Coverage Calculations

## Extract Regions of Interest from BAM File

Extracting only certain regions significantly increases search time.

   $ samtools view input.bam "Chr10:18000-45500" > output.bam
   
## Compute Genome Coverage using Bed Tools

We will use the genomecov tool. [documentation](http://bedtools.readthedocs.io/en/stable/content/tools/genomecov.html)

   $ bedtools genomecov -d -ibam input.bam -g genome.bed > output.coverage
   
There's a script in bin/scripts

