# Analysis of Sequencing Data

## Table of Contents
1. [Introduction](#introduction)
2. [PART 1 - Analysis of Bulk Sequencing](#part-1---analysis-of-bulk-sequencing)
     - [Methodology](#methodology)
     - [Quantification with Kallisto](#quantification-with-kallisto)
     - [Count Conversion](#count-conversion)
     - [Differential Expression Analysis](#differential-expression-analysis)
     - [Results Exploration](#results-exploration)
     - [Conclusion for Part 1](#conclusion-for-part-1)
3. [PART 2 - Analysis of Single-cell Sequencing In Situ](#part-2---analysis-of-single-cell-sequencing-in-situ)
     - [Steps](#steps)
     - [Code & Results](#code--results)
     - [Conclusions for Part 2](#conclusions-for-part-2)

---

## Introduction
# 🧬 Neuro-Genomics Course Project 🧠

Welcome to the Neuro-Genomics Course Project repository! Dive into the complex intersection of genomics and neurology, as we not only explore the mysteries of neurological diseases but also the therapeutic potential of targeted cancer drugs. 📊

## 🚀 Project Overview

- **Part 1**: Grasping Molecular Mechanisms of Neurological Diseases 🧩
  - Unravel the intricate connection between a specific mutated gene and its physiological implications in a neurological disorder.
  - Analyze raw sequencing data from mouse models to gain insights into the molecular deficiencies linked with the gene knockout.
  - Discover potential treatment pathways based on the molecular findings.

- **Part 2**: Dissecting Single Cell Sequencing in Breast Cancer Biopsy 🌌
  - Harness the potential of in situ sequencing to unveil the therapeutic relevance of immunotherapy in breast cancer patients.
  - Classify cells, determine their types, and assess their spatial distributions to gauge the effectiveness of PD-L1 inhibitor drugs.

## 💡 Highlights 

- **Bulk Sequencing Analysis** 🧫: Process and analyze sequencing data from mouse cortex tissues, differentiating between control (C) and knockout (KO) models.
  
- **Single Cell Sequencing** 🦠: Delve into the world of breast cancer at the cellular level, determining the potential effectiveness of immunotherapy drugs based on genetic markers and spatial cell arrangements.

- **Software & Scripting** 💻: Utilize powerful bioinformatics tools like Kallisto, Salmon, Bowtie2, and DeSeq2. We also ensure clarity by justifying software choices and detailing every step with appropriate comments in the scripts.

## 🛠 Tools and Languages

- **Languages**: R
  
- **Bioinformatics Tools**: Kallisto, Situ, and DeSeq2 🧬

## 📁 Dataset 

We work with raw sequencing files, including:
- Control mouse models: `C1.fastq`, `C2.fastq`, `C3.fastq`
- Knockout mouse models: `KO1.fastq`, `KO2.fastq`, `KO3.fastq`
- Breast cancer biopsy: `BreastCancerExpressionInCells.xlsx`, `BreastCancerLocationOfCells.xlsx`

---

## PART 1 - Analysis of Bulk Sequencing

### Methodology

#### Quantification with Kallisto:
[Provide the details here]

#### Count Conversion:
[Provide the details here]

#### Differential Expression Analysis:
[Provide the details here]

#### Results Exploration:
[Provide the details here]

### Conclusion for Part 1
[Provide the conclusion for Part 1 here]

---

## PART 2 - Analysis of Single-cell Sequencing In Situ

The goal of this part is to gain insights into the tissue's molecular characteristics and spatial organization.

### Steps

#### Preprocessing of Gene Expression Data: 
Loaded gene expression data into Seurat, followed by quality control, normalization, and scaling.

#### Principal Component Analysis (PCA):
Conducted PCA to reduce dimensionality and identify principal components capturing major variations.

#### k-NN Clustering:
Applied k-NN clustering to group cells based on gene expression profiles using principal components.

#### Differential Expression Analysis:
Conducted post-clustering differential expression analysis to identify significant patterns and marker genes.


### Code & Results

### Conclusions for Part 2


