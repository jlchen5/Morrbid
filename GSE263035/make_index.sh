#!/bin/bash

STAR --runThreadN 20  \
--runMode genomeGenerate  \
--genomeDir /public1/home/scb0639/reference/refdata-gex-GRCm39-2024-A/2.7.4a  \
--genomeFastaFiles /public1/home/scb0639/reference/refdata-gex-GRCm39-2024-A/fasta/genome.fa  \
--sjdbGTFfile /public1/home/scb0639/reference/refdata-gex-GRCm39-2024-A/genes/genes.gtf