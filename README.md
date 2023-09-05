# NeuroGenomics Project: Analysis of Molecular Deficiencies in Neurological Disease

## Introduction
This project aims to study the molecular deficiencies associated with a specific neurological disease. Through the utilization of a mouse model with a knocked-out gene, the study delves deep into raw sequencing data to understand the implications of this gene and identify markers or indications that provide insights into the disease.

## Objectives:
- Estimate counts using Kallisto.
- Convert these estimated counts from transcript level to gene level.
- Conduct differential expression analysis with DESeq2.
- Analyze results to elucidate the molecular deficiencies.

## Methodology

### 1. Quantification with Kallisto:
Single-end sequencing reads were quantified using Kallisto. The mean was set to 300 with a standard deviation of 50.

### 2. Count Conversion:
Employed the `tximport` R package to transition from transcript-level counts to gene-level counts.

### 3. Differential Expression Analysis:
Using DESeq2, the data was subjected to differential analysis. This involved managing unrecognized genes and structuring the data frame for optimum analysis.

### 4. Results Exploration:
Results were dissected based on the p_adj value, subsequently being segmented into two categories based on the log2FoldChange value.

## Key Findings

- **Myelin's Role**: Significant molecular deficiencies related to the protein Myelin were identified.
  
- **Gene Analysis**: Genes like Neurod6, Myl4, Gtf2i, among others, played crucial roles related to myelin functionality, especially in membranes and myelin sheath integrity.
  
- **Potential Myelin Sheath Issues**: The analysis indicates potential problems with the myelin sheath in treated conditions.

- **Individual Gene Insights**: Genes such as Hapln4, Lingo3, LRG1, and TRF2 were identified as significant in the treated condition, hinting at their roles in stabilizing the extracellular matrix and potential issues with the myelin sheath.

## About Myelin

Myelin is a protective layer surrounding nerves, fundamental for the efficient transmission of electrical impulses in nerve cells. Damages or deficiencies in Myelin can lead to diseases like Multiple Sclerosis (MS), Guillain-Barre Syndrome (GBS), and more.

## Potential Treatments

There are no definitive cures for myelin damage, but treatments such as Interferon beta-1a, Fingolimod, and the use of OPCs can help control symptoms, reduce inflammation, and potentially promote remyelination.

---

For a detailed analysis, refer to the provided Jupyter Notebook, attached data, and the in-depth report.

