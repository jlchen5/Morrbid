#!/bin/bash

cellranger count --id=heart_control_rep1 \
   --fastqs=/public1/home/scb0639/project/GSE263035 \
   --sample=heart_control_rep1 \
   --transcriptome=/public1/home/scb0639/reference/refdata-gex-GRCm39-2024-A \
   --create-bam false
