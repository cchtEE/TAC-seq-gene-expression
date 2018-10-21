# TAC-seq gene expression
This repository contains TAC-seq gene expression data analysis software.

## Input
* [TAC-seq data analysis](https://github.com/cchtEE/TAC-seq-data-analysis) output file with:
  * column `sample`
  * column `locus`
  * column `molecule_count`
  
* Target file with:
  * column `id` - values have to match with `locus` values in [TAC-seq data analysis](https://github.com/cchtEE/TAC-seq-data-analysis) output file
  * column `type` - possible values: `biomarker`, `housekeeper` or `spike_in`

## Output
Normalized count table with:
* column `sample`
* column for each target `id`
