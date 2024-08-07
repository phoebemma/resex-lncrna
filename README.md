# Effects of resistance exercise training on the long non-coding RNA
landscape in young individuals


## Authors

Chidimma Phoebe Echebiri, Rafi Ahmad, Stian Ellefsen, Daniel Hammarstöm

## Repository organisation

`/Manuscript`

- My Library.bib
- lncRNA-manuscript.qmd
- science.csl

`/resources` - `rna-biology-csl.csl` Current Citations style file -

`/R`

- Contains the different scrpts used in analyses. From data extraction
  to final step

To reproduce the analyses, the scripts should be run in the following
order

*uncomment the “saveRDS..” lines to save the extracted files in the
locations needed for downstream analyses* 1. `Trainome_functions.R`.
This contains functions written from scripts that were repeatedly needed
for all Trainome analyses on R

2.  `contratrain_data_extraction.R`. This downloads the contratrain
    data, extracts and preprocesses the relevant data needed for
    downstream analyses

3.  `ct_data_exploration.R`. This runs basic exploratory anaylses on the
    extracted data

4.  `Model_using_seq_wrapper.R`. This runs and saves the outputs of the
    four models

5.  `Models_exploration.R`. This contains basic exploration of the the
    obtained models. The lower part of the script is commented out, and
    contains coexpression analyses of the lncRNAs and mRNA genes

6.  `Volume_models.R` : Basic exploration of the volume models

7.  `training_models.R` Contains basic exploration of the training
    models

8.  `co-expression_model.R` Contains the script for co-expression
    analyses of the lncRNAs and mRNA genes by modeling

9.  `/archived_scripts` Contains scripts not needed necessary for
    reproduction of the results, but might be relevant for future
    reference

`/README_files/libs`
