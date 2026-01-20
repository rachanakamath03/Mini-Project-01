# Mini-Project-01

## Project Overview:
Type 2 Diabetes is a complex metabolic disorder characterisded by insulin resistence, chronic inflammation and systemic metabolic dysregulation. Transcriptomic profiling sing RNA sequencing (RNA-seq) provides a powerful approach to identify gene expression changes associated with disease pathology and to uncover poteintial biomarkers for early detection and risk stratification. This mini-project performs a comparative transcriptomic analysis of pancreatic islets RNA-seq data from individuals with Type Diabetes and healthy controls to identify differentially expressed genes (DEGs) and biologically relevant pathways associated with T2D.

# Objectives:
- To identify differentially expressed genes between T2D patients and healthy controls using bulk RNA-seq data
- To characterise transcriptomic patterns associated with inflammation, immune regulation, and metabolic dysfunction.
- To perform functional enrichment analysis to interpret disease-relevant biological pathways.

# Dataset Information:
Source: NCBI GENE EXPRESSION OMNIBUS (GEO) 
Accession ID: GSEGSE86468 
Organism: Homo sapiens 
Tissue: Pancreatic Islet 
Study Design: Case-Control (Type 2 Diabetes vs healthy controls) 
Technology: 
Data Used: Processed RSEM raw expected counts and sample metadata

Raw sequencuing files (FASTQ) were not used; publicaly available processed count data were analysed to ensure efficient and reproducible downstream analysis

# Analysis Workflow:
1. Data acqusition and organisation
2. Sample metadata curation and condition assignment
3. Preprocessing and low-count gene filtering
4. Differential expression analysis using DESeq2
5. Quality control and exploratory analysis (PCA, clustering)
6. Visulalization of DEGs ( volcano plot, heatmaps)
7. Functional enrichment analysis (GO and KEGG pathways)
8. Biological interpretation and disease relevance assessment

# Tools and Software
R (statistical computing)
