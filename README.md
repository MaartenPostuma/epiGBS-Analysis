# epiGBS-Analysis
Perfrorms basic analysis for epiGBS, includes:
* reading data into methylKit
* filtering individuals on number of sites covered
* filtering sites based on mincoverage (10x)
* filtering sites based on individual missingness\
* Creating PCA's
* Creating overall methylation plots
* Loading filtered data into DSS

# prerequisits
This analysis requires
* Linux
* conda
Further R packages can be installed using conda
`conda env create -f src/env/methylkit.yaml`

# input

This pipeline runs on the output of epiGBS (https://github.com/nioo-knaw/epiGBS2). Next to this it requires a tab seperated Design file with a column with the sampleIDs (which should be the same as the samples input into epiGBS) is required.
This Design file can also include all treatments used in the experiment in order to colour the PCA and percentage methylation plots, and be used to run the DSS models.

```
sampleID treatment1 treatment2 etc.
sample_1 drought    heat
sample_2 rain       cold
```

# How to run
The scripts perform general analysis on the epiGBS pipeline, however input data, filtering contstraints and variables to colour plots / populate models are not general and need to be adjusted when running the scripts. 

The scripts in the pipeline need to be ran in the following order:
* src/readMethylKit.R
* src/filterMethylKit.R
* src/basicAnalysisMethylKit.R
* src/runDSS.R

# Output
The scripts output the following:
 * 4 methylKit directories: (methylKitDB (unfiltered), FiltCG, FiltCHG and FiltCHH)
 * results/nLociDistrib.png (distribution of number of loci per individual and the cutoff for filtering)
 * results/design_filtered.tsv (The design file filtered to only include samples passing filtering)
 * results/percMeth.png (a plot showing the different percentages of methylation in different contexts/treatments)
 * results/pca.png      (PCA plots for the different contexts which can be coloured using treatments from the Design.tsv)
  * 
