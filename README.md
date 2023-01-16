# Stone

## Contents
### Code: [Preview](https://htmlpreview.github.io/?https://github.com/kkwock/Stone/blob/main/ncbiGeneSumm.html)
* `ncbi-ann-glyma.R`: Rscript with source functions
* `ncbiGeneSumm.Rmd`: R Notebook with raw code
* `ncbiGeneSumm.html`: HTML output of notebook
* `ncbiGeneSumm.RData`: Environmental Variables of pre-processed data.

### Data
* `data/Annotations/`: CSV outputs of Gene annotations and GO enrichments
* `data/Annotations/images`: Images of GO-enrichments

## How to use
1. Use `ncbiGeneSumm.Rmd` notebook in RStudios
2. Load in custom functions: `source('ncbi-ann-glyma.R')` 
3. Load variables used in rendering outputs: `load('ncbiGeneSumm.RData')`
