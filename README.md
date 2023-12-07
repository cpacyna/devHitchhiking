# devHitchhiking
This repo contains code for the paper "Multifocal, multiphenotypic tumours arising from an MTOR mutation acquired in early embryogenesis."

You can find Rmarkdowns that report our analysis methods:
- *pedRenalCaseStudy_1SNVfiltering.Rmd* uses several filtering tests to identify somatic SNVs.
- *pedRenalCaseStudy_2sharedVariants_fig1de.Rmd* takes corrected variant allele frequencies from Supplementary Table 2 and generates figures 1D and 1E.
- *pedRenalCaseStudy_3aneuploidyHeatmap_exDatFig1a.Rmd* creates extended data figure 1a based on ASCAT outputs for the four tumor samples.
- *pedRenalCaseStudy_4bulkRNAanalysis_exDatFig1d.Rmd* contains analysis for the bulk RNA-seq data, including differential expression and gene set enrichment analysis.
- *pedRenalCaseStudy_5scRNAanalysis.Rmd* contains analysis for the single cell RNA-seq data, including standard data processing, cell type annotation with logistic regression, and tumor cell identification with alleleIntegrator.

We also include some data:
- *ascat* contains the copy number variation results from ASCAT.
- *rna* contains both case and GTEx bulk RNA-seq counts as input for *pedRenalCaseStudy_4bulkRNAanalysis_exDatFig1d.Rmd*.

Lastly, we share two scripts with details on how we ran alleleCounter as an input for the shearwater-like filtering step in *pedRenalCaseStudy_1SNVfiltering.Rmd*.
- *alleleCounter_adultOncos.sh*
- *alleleCounter_pedsOncos.sh*

This project's raw data can be found on EGA and the Seurat object can be found on figshare (see paper for access details).

  


  
