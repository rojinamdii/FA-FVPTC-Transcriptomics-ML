# Microarray-Based Transcriptomic and Machine Learning Analysis of FA vs FVPTC

This repository contains the R scripts used for the transcriptomic, network-based, and machine-learning analyses reported in the following study:

> **Microarray-based Transcriptomic and Machine Learning Analysis Reveals Candidate Diagnostic Biomarkers Separating Follicular Adenoma and Follicular Variant Papillary Thyroid Carcinoma**

## Manuscript

**Authors:** Zhina Mohamadi, Rozhina Mohammadi, and Armin Zarinkhat  
**Equal contribution:** Zhina Mohamadi and Rozhina Mohammadi contributed equally and share first authorship.  
**Publication status:** Accepted. The journal citation and DOI can be added when available.

## Study overview

The study investigated molecular differences between follicular adenoma (FA) and source-annotated follicular variant papillary thyroid carcinoma (FVPTC) by integrating two public microarray cohorts. The workflow combined differential-expression analysis, functional enrichment, protein–protein interaction (PPI) network analysis, hub-gene prioritization, and internal machine-learning evaluation.

The analysis included **43 samples**: 24 FA and 19 FVPTC. FVPTC labels were retained from the original GEO annotations. The cases could not be retrospectively reclassified under contemporary WHO/NIFTP criteria because the required histological slides and detailed information on encapsulation, invasion, and tumor-interface pattern were unavailable. The findings therefore apply to the source-annotated, historically defined FVPTC group and should not be interpreted as representing one specific contemporary diagnostic subtype.

## Public datasets

The raw microarray data are publicly available from the NCBI Gene Expression Omnibus (GEO):

| GEO accession | Platform | FA samples | FVPTC samples | Total |
| --- | --- | ---: | ---: | ---: |
| [GSE29315](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29315) | [GPL8300](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8300) | 17 | 13 | 30 |
| [GSE54958](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54958) | [GPL6244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6244) | 7 | 6 | 13 |
| **Combined cohort** | — | **24** | **19** | **43** |

The raw CEL files are not included in this repository. They must be obtained from GEO using the accession numbers above.

## Analytical workflow

1. **Microarray preprocessing:** Raw CEL files were background-corrected and normalized within each dataset using Robust Multi-array Average (RMA). Probe identifiers were mapped to gene symbols using the corresponding platform annotation files. When multiple probe sets mapped to the same gene symbol, their expression values were averaged.
2. **Dataset integration and batch correction:** The two gene-level matrices were aligned to their shared genes and combined. Batch effects associated with the source dataset were adjusted using ComBat, while histological condition (FA vs FVPTC) was included in the model to preserve the biological contrast.
3. **Differential expression:** The limma framework was used to compare FVPTC with FA. P-values were adjusted using the Benjamini–Hochberg false-discovery-rate procedure. Differentially expressed genes were defined using |log2 fold change| > 1 and adjusted p-value < 0.05.
4. **Functional and network analyses:** Gene Ontology and KEGG enrichment analyses were performed. A PPI network was obtained from STRING using a medium-confidence score threshold of 0.4, followed by module detection and centrality-based hub-gene prioritization in Cytoscape.
5. **Machine-learning screening:** Seven classifiers were screened: logistic regression, linear support vector machine, decision tree, random forest, naive Bayes, k-nearest neighbors (KNN), and XGBoost.
6. **Final internal evaluation and candidate prioritization:** XGBoost, naive Bayes, and KNN were evaluated using nested repeated cross-validation (outer stratified five-fold cross-validation repeated 20 times, with inner four-fold tuning). Gene contributions were assessed using cross-fitted permutation importance on outer held-out samples. Model-normalized Borda ranking, rank-stability analysis, gene-correlation assessment, and panel-size sensitivity analysis were used to prioritize candidate genes.

## Main findings reported in the manuscript

- **74 differentially expressed genes** were identified: 61 upregulated and 13 downregulated in FVPTC compared with FA.
- Eight hub-derived candidate genes were carried forward to machine-learning evaluation: **FN1, KRT19, TIMP1, CHI3L1, CLU, SPP1, APOE, and EFEMP1**.
- The primary five-gene candidate panel was **CHI3L1, FN1, EFEMP1, CLU, and KRT19**.
- **APOE, TIMP1, and SPP1** were retained as secondary candidates.
- The final models had internal ROC-AUC values of 0.886 for XGBoost, 0.904 for naive Bayes, and 0.906 for KNN. No statistically significant pairwise differences in ROC-AUC were detected between these model families.
- In the panel-size sensitivity analysis, the five-gene panel had the highest mean internal AUC (0.916) among the evaluated panel sizes.

These results represent **internal, exploratory candidate prioritization**. The study did not include an independent external validation cohort or orthogonal experimental validation. The candidate genes should therefore not be presented as clinically validated diagnostic biomarkers.

## Repository contents

| File | Description |
| --- | --- |
| `Microarray preprocessing.R` | RMA preprocessing and probe-to-gene annotation for GSE29315. |
| `Dataset integration & DEG analysis.R` | Processing of GSE54958, shared-gene alignment, cohort integration, ComBat batch correction, and limma differential-expression analysis. |
| `Machine learning screening.R` | Initial seven-classifier screening and the preliminary top-three model workflow. |
| `Final ML evaluation & gene prioritization.R` | Corrected nested repeated cross-validation, cross-fitted permutation importance, consensus ranking, rank-stability, correlation, and panel-size sensitivity analyses. |
| `Publication figures.R` | Generation of publication figures from the analysis outputs. |
| `FA-FVPTC.Rproj` | RStudio project file. |
| `renv.lock` | Record of the R package environment. |
| `renv/activate.R` | Activates the project-specific `renv` environment. |

## Required input files

The following input files are required to reproduce the preprocessing steps but are not included in this repository:

- Raw CEL files for the selected samples in GSE29315 and GSE54958.
- `GPL8300-10787.csv`
- `GPL6244-17930-modified.csv`

The two annotation CSV files must be the same versions used in the analysis. The current R scripts read the annotation files by filename from the R working directory.

**Input-path note:** The preprocessing scripts call `ReadAffy()` without an explicit CEL-file path. Therefore, when running each preprocessing step, make sure the intended dataset's CEL files are available in the active R working directory, or update the CEL-file path in the script to match your local folder structure. Do not place both cohorts together in the working directory unless the script is configured to select the intended files.

## Suggested execution order

1. Open `FA-FVPTC.Rproj` in RStudio.
2. Run `Microarray preprocessing.R` with the GSE29315 CEL files and its annotation file available.
3. Run `Dataset integration & DEG analysis.R` with the GSE54958 CEL files and its annotation file available. This step uses the GSE29315 expression object created in the preceding step.
4. Run `Machine learning screening.R`.
5. Run `Final ML evaluation & gene prioritization.R`.
6. Run `Publication figures.R`.

The downstream scripts depend on intermediate objects and output files generated by earlier steps. In particular, the final evaluation and figure-generation scripts require outputs produced by the machine-learning screening workflow.

## Software environment

The analyses were performed in **R version 4.5.1**. The `renv.lock` file records the package environment. After opening the project in RStudio, restore the recorded packages with:

```r
renv::restore()
```

Restoring the R package environment does not download the GEO CEL files or the annotation CSV files.

## Citation

Please cite the associated manuscript:

Mohamadi Z, Mohammadi R, Zarinkhat A. *Microarray-based Transcriptomic and Machine Learning Analysis Reveals Candidate Diagnostic Biomarkers Separating Follicular Adenoma and Follicular Variant Papillary Thyroid Carcinoma.*

The journal, volume, pages, and DOI can be added here when the final publication details are available.
