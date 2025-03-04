# Morrbid
Bulk RNA-seq and scRNA analyisis for Morrbid (ncRNA gene)


Morrbid (myeloid RNA regulator of BCL2L11 induced cell death) in *Mus musculus* and its Orthologs is MIR4435-2HG in *Homo sapiens*


- Raw reads were filtered by `cutadpt v1.8.1,-m 20 --trim-n --max-n 0.7 -q 20`
- Contaminants were filtered by `fastq_screen (v0.5.1. 4)`
- Clean reads were mapped to UCSC `mm10` or `hg38` genome via `STAR v2.5.2b`
- Then assigned by `featureCounts (subread-1.4.6-p1 command line)`
- Genome_build: `mm9`
- Genome_build: `hg38`



## 1 Bulk RNA-seq analysis

- step1 download the expression matrix
- step2 differential expressed genes analysis
- step3 plot volcano and heatmap
- step4 function analysis


## 2 scRNA-seq analysis
The single cell RNA-seq data was analyzed using [this script](https://github.com/jlchen5/Morrbid/blob/main/Untitled.ipynb).



## Methods

RNA-seq Data Processing
Raw RNA-seq data were obtained from the Gene Expression Omnibus (GEO) database under accession number GSEXXX. All genomic interval visualizations were generated using the R packages.

Single-Cell RNA-seq Analysis
Single-cell RNA-seq data (GEO:GSEXXX ) were processed as follows. Cells were filtered using thresholds informed by empirical distributions: Cells with <500 total detected genes or <1,000 total UMI counts (indicative of low RNA capture). Cells with >15% mitochondrial gene content (identifying stressed/apoptotic cells). Genes detected in <3 cells were excluded to reduce noise. The top 3,000 highly variable genes (HVGs) were identified using Seurat’s FindVariableFeatures. Dimensionality Reduction:Data were scaled (z-score normalization), and principal component analysis (PCA) was performed on HVGs. The top 20 principal components (PCs), selected via elbow plot and Jackstraw permutation testing (p < 0.001), were retained. Sample integration was achieved using Harmony with default parameters to mitigate batch effects. A shared nearest neighbor (SNN) graph was constructed and clustered via the Leiden algorithm. Uniform Manifold Approximation and Projection (UMAP) was applied to Harmony-adjusted PCs for visualization. All analyses were implemented in R (v4.1.3) .


## License

[MIT license](https://github.com/jlchen5/Morrbid/blob/main/LICENSE).
