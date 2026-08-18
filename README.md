# Sunflower Fatty Acid Stability GWAS

Analysis code for:

> **The stability of fatty acid composition in sunflower oil is dependent on environment and affected by structural variation**
> Markus Ingold, Qingming Gao, Jennifer R. Mandel, James P. McNellie, Kyle G. Keepers,
> Jessica G. Barb, John M. Burke, Loren H. Rieseberg, and Brent S. Hulke
> Submitted to *Theoretical and Applied Genetics*, Special Issue on Breeding for Improved
> Crop Quality. Manuscript in review (2026).

## Overview

In sunflower (*Helianthus annuus* L.), the seed fatty acid (FA) profile — primarily oleic,
linoleic, stearic and palmitic acid — determines oil quality, but the genetic basis of the trait
and its interaction with environment is poorly understood.

This repository contains the scripts used to relate FA composition and its **stability** to
genotype across an extensive geographic cline. Phenotypes come from the sunflower association
mapping (SAM) population — 287 varieties covering ~90% of sunflower diversity — grown in eight
field trials in North America: Vancouver BC (2010), Moorhead MN (2015; early and late plantings
2016), Ames IA (2010, 2013, 2014), and Athens GA (2010). Genotypes are 2.2M high-quality SNPs
(MAF > 5%) called against the HA412HOv2.0 reference.

Multivariate GWAS was run over four phenotype sets:

1. mean FA composition within each environment;
2. mean FA composition within each environment, omitting high oleic (HO) varieties;
3. **α stability** — trait stability *within* environments, quantified by the standard error among
   replicate samples;
4. **β stability** — trait stability *across* environments, via Eberhart and Russell's β.

Key results the code supports: HO varieties show high β stability, producing substantial marker
overlap between analyses (1) and (4); 147 candidate genes were identified around significant SNPs,
including 15 FA-metabolism genes, 6 heat shock proteins and 22 transcription factors; and a large
Chromosome 5 introgression with two flanking inverted sequences coincides with α stability in the
Georgia (hot, humid) trial.

## Repository layout

Directories follow the analysis stages. Scripts are listed roughly in the order they are run.

### `Phenotype_Analysis/`
| Script | Purpose |
| --- | --- |
| `phenotype_prep.R` | Reads raw FA measurements; drops phenotypes based on a single sample; computes per-environment means and standard errors (α stability); builds GGE biplots; writes the per-environment phenotype and standard-error CSVs used as GWAS input, plus the FAD2-1 mutant status list. |

### `Stability/`
| Script | Purpose |
| --- | --- |
| `Stability_regression.R` | Eberhart & Russell (1966) β stability. Regresses each genotype's FA value on the environmental index, extracting βᵢ per fatty acid; drops genotypes missing from more than three trials or from the Georgia trial; Yeo-Johnson power transformation of the output; renames samples to match the VCF. |
| `stability_boxplot.R` | Boxplots of β stability by FAD2-1 mutant status (FAD WT vs FAD mutant) with t-tests — the paper's Figure 3. |

### `Climate/`
| Script | Purpose |
| --- | --- |
| `chelsa_download.sh` | Downloads one month/year of CHELSA v2.1 monthly rasters (`clt cmi hurs pet pr rsds sfcWind tas tasmax tasmin vpd`). |
| `chelsa_download_files.sh` | Driver that calls `chelsa_download.sh` for every planting month and the two subsequent months across all trial years. |
| `climate_data_extraction.R` | Extracts the climate variables at each trial's coordinates from the CHELSA rasters (`raster`), producing `climate_data.csv`. |
| `climate_data_analysis.R` | Correlations (`corrplot`) and RDA/PCA of FA composition and α stability against the climate variables, latitude and longitude — the paper's Figure 1. |

### `Genotype_Analysis/`
| Script | Purpose |
| --- | --- |
| `ADMIXTURE.sh` | Runs ADMIXTURE for K = 1–10 with cross-validation on the LD-pruned PLINK bed. |
| `ADMIXTURE_results.R` | Plots CV error vs number of clusters (Figure S1) — no population structure was found, so the GEMMA kinship matrix was deemed sufficient. |
| `kinship_matrix.R` | Builds the VanRaden kinship matrix from the VCF dosage matrix (`vcfR`, `AGHmatrix`) for GEMMA. |
| `PCA_Dendrogram_etc.R` | Genotype PCA and hierarchical clustering / dendrograms of the SAM panel. |
| `fad-mutant_analysis.R` | Analyses the VCF subset around the FAD2-1 mutation used to classify high oleic varieties. |
| `gene_density.R` | Gene density histograms per chromosome from the HA412 GFF3, used to locate centromeres. |
| `post_GWAS_LD_blocks.R` | LD estimation (`ldsep`, `updog`) over the significant-variant VCF subset. |

### `GWAS_inputs/`
| Script | Purpose |
| --- | --- |
| `extract_high_oleic_samples.sh` | `bcftools view` subset of the 28 high oleic (FAD2-1 mutant) SAM accessions. |
| `gemma_script.sh` | The multivariate GWAS driver: builds the kinship matrix, then runs GEMMA `-lmm 1` with all four fatty acids (`-n 1 2 3 4`) for every environment, for the no-FAD subsets, for α stability (standard errors) and for β stability. |
| `vcf2gwas.sh` | Univariate GWAS via vcf2gwas for each environment. **Note:** vcf2gwas cannot do multivariate analysis — it silently uses only the first phenotype, which is why multivariate runs go through GEMMA directly. Output folders must be renamed between runs or they are overwritten. |

### `post-GWAS_files/`
| Script | Purpose |
| --- | --- |
| `post_GWAS_analysis.R` | The main results script: reads all GEMMA `.assoc.txt` outputs, applies the Bonferroni thresholds (over SNP count and over LD-block count), and draws the Manhattan/QQ/circular plots (`qqman`, `CMplot`) and the Venn diagram (`VennDiagram`) comparing analyses. |
| `create_results_table.R` | Assembles the significant-SNP results table across all environments and analyses (Table S4). |
| `create_genotype_matrices.R` | Extracts genotype matrices for the multivariate significant SNPs from the VCF. |
| `overlap_relaxed_beta_no_FAD.R` | Overlap of β stability signals with the no-FAD analyses under the relaxed (LD-block) threshold. |
| `DEG_analysis.R` | DESeq2 on the SRP092742 tissue expression matrix to find genes with log2 fold-change > 2 in developing seed vs 10 other tissues — the filter behind the candidate gene list. |
| `candidate_gene_function.R` | Pulls genes within a window of each significant SNP from the GFF3 and annotates their predicted function (MapMan, `biomaRt`). |
| `introgression_analysis.R` | Compares significant SNPs against the introgressions reported by Todesco et al. (2020). |
| `compare_to_Chernova.R` | **Not used in the current paper version** — Chernova et al. SNP positions did not lift over correctly. Kept for provenance. |

### `Introgression_Analysis/`
| Script | Purpose |
| --- | --- |
| `chr5_introgression_analysis.R` | Tests the Chromosome 5 introgression haplotypes against α stability (notably the Georgia trial), using two genotype PCs as covariates. |
| `Introgression_GWAS_analysis.R` | Visualises the GEMMA output from the introgression-marker GWAS. |

## Input data

The scripts expect these inputs, which are **not** included here (they are third-party or too
large for version control):

- **SNP genotypes** — downloaded from [HelianthOME](http://www.helianthome.org/download/#genotype)
  (Bercovitch et al., 2022), called against HA412HOv2.0 (Huang et al., 2023) with GATK best
  practices. Filtered for single-copy sites by depth, `minQ=100`, `MAF>=0.05`, `max-missing=0.9`,
  imputed with BEAGLE 5.3, converted to PLINK with PLINK v1.90b6.27. Referred to throughout as
  `IMPUTED_SAM_GATK_simple_biallelic_noSingleCopy_maxM0.9_maf0.05_minQ100_recode.vcf`.
- **Phenotypes** — gas chromatography FAME measurements (Hulke et al., 2010) of palmitic, stearic,
  oleic and linoleic acid, 20 seeds per head, for the SAM panel in the eight trials.
- **Climate** — CHELSA v2.1 monthly rasters (accessed 20 April 2026; Karger et al. 2017a, 2017b),
  downloaded by the scripts in `Climate/`.
- **Annotation** — `HAN412_Eugene_curated_v1_1.gff3`, and MapMan `X4.2_helianthus_annuus.txt` from
  [mapman.gabipd.org](https://mapman.gabipd.org).
- **Expression** — `SRP092742.salmon.count.tsv` tissue count matrix (Badouin et al., 2017).
- **Introgressions** — `SAM_introgression.pcadmix.txt` and the Chr5 haplotype table from
  Todesco et al. (2020).

## Software

External tools, at the versions used in the paper:

| Tool | Version |
| --- | --- |
| GEMMA | 0.98.5 |
| vcf2gwas | 0.8.7 |
| ADMIXTURE | 1.3.0 |
| PLINK | v1.90b6.27 |
| BEAGLE | 5.3 |
| BCFtools / SAMtools | — |
| BLAST (`blastn`) | — |
| LDBlockShow | 1.40 |

R packages used across the scripts:

*CRAN* — `AGHmatrix`, `ape`, `arsenal`, `cellWise`, `circlize`, `cluster`, `CMplot`, `corrplot`,
`data.table`, `dendextend`, `dplyr`, `factoextra`, `GGally`, `gge`, `ggfortify`, `ggplot2`,
`ggpubr`, `ggrepel`, `gridExtra`, `gtools`, `hrbrthemes`, `ldsep`, `metan`, `naniar`, `pegas`,
`PerformanceAnalytics`, `pheatmap`, `poppr`, `psych`, `purrr`, `qqman`, `qqplotr`, `raster`,
`readxl`, `reshape`, `reshape2`, `scales`, `sp`, `stringr`, `tidyr`, `tidyverse`, `updog`,
`vcfR`, `vegan`, `VennDiagram`, `viridis`, `xtable`, `zoo`

*Bioconductor* — `biomaRt`, `DESeq2`, `GenomicDistributions`, `GenomicRanges`, `GWASTools`,
`Repitools`, `SeqArray`, `snpStats`, `VariantAnnotation`

## Running the code

These are the analysis scripts as they were run for the paper, not a packaged pipeline. They are
kept in their original form for provenance, which means:

- **Paths are hardcoded.** Each script begins with a `setwd()` pointing at one of the original
  authors' machines (Linux and Windows both appear). Change the `setwd()` call and the relative
  input paths before running.
- Scripts are meant to be stepped through interactively in R rather than run start-to-finish.
- Some blocks are commented out — exploratory alternatives, or steps that were superseded.

Rough order of execution:

1. `Phenotype_Analysis/phenotype_prep.R` → per-environment means and standard errors
2. `Stability/Stability_regression.R` → β stability
3. `Genotype_Analysis/kinship_matrix.R` (and `ADMIXTURE.sh` to confirm no structure)
4. `GWAS_inputs/gemma_script.sh` (multivariate) and `GWAS_inputs/vcf2gwas.sh` (univariate)
5. `post-GWAS_files/post_GWAS_analysis.R` and the remaining `post-GWAS_files/` scripts
6. `Climate/`, `Introgression_Analysis/` and `Stability/stability_boxplot.R` for the supporting figures

## Citation

If you use this code, please cite the paper (see `CITATION.cff`). The citation entry will be
updated with the volume, pages and DOI once the paper is published.

## Funding

This work was supported by NSF grant 2308194; USDA-ARS CRIS projects 3060-21000-043-00D and
3060-21000-047-00D; USDA National Institute of Food and Agriculture grant 2008-35300-19263; and
Genome Canada and Genome BC's Applied Genomics Research in Bioproducts or Crops (ABC) Competition.

Mention of trade names or commercial products in this repository is solely for the purpose of
providing specific information and does not imply recommendation or endorsement by the U.S.
Department of Agriculture. USDA is an equal opportunity provider and employer.

## License

Released under the MIT License — see `LICENSE`.
