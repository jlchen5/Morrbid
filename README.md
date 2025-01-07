# Morrbid
Bulk RNA-seq and scRNA analyisis for Morrbid (ncRNA gene)


Morrbid (myeloid RNA regulator of BCL2L11 induced cell death) in *Mus musculus* and its Orthologs is MIR4435-2HG in *Homo sapiens*


- Raw reads were filtered by `cutadpt v1.8.1,-m 20 --trim-n --max-n 0.7 -q 20`
- Contaminants were filtered by `fastq_screen (v0.5.1. 4)`
- Clean reads were mapped to UCSC `mm10` or `hg38` genome via `STAR v2.5.2b`
- Then assigned by `featureCounts (subread-1.4.6-p1 command line)`
- Genome_build: `mm9`
- Genome_build: `hg38`



## Step1 Bulk RNA-seq analysis

## Step2 scRNA-seq analysis
