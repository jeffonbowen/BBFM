# BBFM

## Import

There are scripts for import and processing of the species data and environmental covariates.

Each of the two scripts write the processed data to "data\_processed" so that the import scripts don't have to run for each new session.

## Data Formats

Most or all of the biodiversity analysis require the data in vegan format. Will be using straight detection count (i.e. no QPAD corrections). However

For single-species analysis, can use the data in long format. However, need to add zeros. convert to wide, add zeros and then back to long. Could make a function!

## Random Notes

-   
