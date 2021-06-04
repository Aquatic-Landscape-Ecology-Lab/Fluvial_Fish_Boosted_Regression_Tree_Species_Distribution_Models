# BRTs
R code for generating fluvial fish species distribution models using Boosted Regression Trees (BRTs)

## Contacts:
Hao Yu (yuhao91@msu.edu); Arthur Cooper (coopera6@msu.edu)

## Purpose:
This R code was developed to generate species distribution models for fluvial fish species based on their native ranges using Boosted Regression Trees (BRTs) as the modeling approach (see Elith et al. 2008). This work was done in support of the USGS Aquatic GAP Project, with example applications of these BRT models found in Yu et al. 2020 and Cooper et al. 2019.

## Requirements:
R 4.0.2 (or newer versions) with ‘dismo’ and ‘labdsv’ packages installed.

## Inputs:
There are five input .csv tables required to run this code, including tables containing fish community samples linked to the NHDPlusV2.1 dataset, a list of fish species to be modeled, fish range data for USGS HUC8s, HUC8 designations for NHDPlusV2.1 COMIDs, and predictors used in model development. These files are described below.

•	Fish_data.csv: a community fish sample file containing species-level abundance in ‘stacked’ format with species data in rows (as opposed to columns). At a minimum this table should include the following fields:

     -	nfhp_id: unique sample identifier
     -	comid: unique flowline identifier from the NHDPlusV2.1 dataset
     -	common_name: fish species common name
     -	scientific_name: fish species scientific name
     -	itis_species: ITIS (Integrated Taxonomic Information System) species level taxonomic serial number (TSN)
     -	sp_count: the raw species abundance (count) for a sample
     -	restricted: value of 1 assigned to records that cannot be publicly shared (and not used in modeling) or null if not restricted
     
•	Fish_list.csv: a table containing a list of fish species to model with the following fields:

     -	Scientific_name: fish species scientific name
     -	ITIS: ITIS (Integrated Taxonomic Information System) species level taxonomic serial number (TSN)
     
•	HUC8_ranges.csv: a file containing species range data (for species in Fish_list.csv) from the USGS Non-indigenous Aquatic Species (NAS) program (https://nas.er.usgs.gov/).  This file should contain the following fields:

     -	Scientific_name: fish species scientific name
     -	HUC8: the level 8 USGS Hydrologic Unit Code (HUC) that the NHDPlusV2.1 flowline falls within
     -	OriginStatus: designates HUC8s as being ‘Native’ for native range and ‘Introduced’ for non-native range
     
•	HUC8_streams.csv: a table relating HUC8s to each NHDPlusV2.1 COMID for the conterminous U.S.:

     -	COMID: unique flowline identifier from the NHDPlusV2.1 dataset
     -	HUC8: the level 8 USGS Hydrologic Unit Code (HUC) that the NHDPlusV2.1 flowline falls within
     
•	Predictors.csv: a table of natural and anthropogenic fluvial predictors used in model development. For example, see: https://www.sciencebase.gov/catalog/item/58584d7ee4b0e40e53c24135

## Outputs:
There are seven output files generated when running this code, including X tables, a pdf file containing partial dependence plots, and an R data file that saves the R workspace. These files are described below.

•	BRT_Stats.csv: a table with basic BRT parameters, including species names, presences, absences, prevalence, tree complexity, learning rate, number of trees in the final model, deviance explained, training AUC, and cross validation AUC. This AUC is calculated by averaging 10 AUCs from 10 testing datasets in 10-fold cross-validation. 

•	BRT_VarContributions.csv: a table containing BRT predictor variable contributions for each species.

•	BRT_CV.csv: a file containing species model diagnostic and threshold information, including species name, deviance, presences, absences, AUC, correlation, threshold, sensitivity, specificity, and TSS. In this output file, AUC is calculated from one file by combining the 10 test datasets from 10-fold cross validation.This AUC value may be slightly different from that in the BRT_Stats.csv file above. 

•	BRT_CV_predict_all_[scientific name]_[common name].csv: a table containing  COMID, abundance, HUC8, common name, all predictor variables, PA(presence/absence), predicted probability, predicted PA within the sampled streams.

•	BRT_[scientific name]_[common name]_native_prediction.csv: logistic predictions for all stream reaches within the species’ native range.

•	BRT_[scientific name]_[common name]_plots.pdf: a file with partial dependence plots for the top 12 predictor variables for each species.

•	BRT_[scientific name]_[common name].RData: an R data file containing the BRT model outputs for each species.


## References:

Elith, J., Leathwick, J. R., & Hastie, T. (2008). A working guide to boosted regression trees. Journal of Animal Ecology, 77(4), 802-813.  https://doi.org/10.1111/j.1365-2656.2008.01390.x

Yu, H., Cooper, A. R., & Infante, D. M. (2020). Improving species distribution model predictive accuracy using species abundance: Application with boosted regression trees. Ecological Modelling, 432, 109202. https://doi.org/10.1016/j.ecolmodel.2020.109202.

Cooper, A. R., Tsang, Y. P., Infante, D. M., Daniel, W. M., McKerrow, A. J., & Wieferich, D. (2019). Protected areas lacking for many common fluvial fishes of the conterminous USA. Diversity and Distributions, 25(8), 1289-1303. https://doi.org/10.1111/ddi.12937.

## Copyright and License:

This is considered to be in the public domain, and is licensed under unlicense_

.. _unlicense: https://unlicense.org/

## Acknowledgements:

We thank the USGS Aquatic GAP program for funding this effort.
