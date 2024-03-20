# Calling ROHs

In this directory, we infer runs of homozygosity using two different approaches: using the `bcftools -roh` algorithm which uses a hidden Markov model to estimate autozygosity, and the `plink --homozyg` function which uses a sliding window approach. 

In script number 3 (`3_roh_calling_plink.R`), the function defined in `function_run_rohcalling.R` is called which not only calculates ROH/FROH using plink, but also plots the inferred ROHs in a circos plot. This function was defined to allows easy comparison between different parameter setting combinations. Additionally, we calculate genome-wide sMLH (standardised multi-locus heterozygosity) in script 4 `4_calculate_sMLH.R` (not reported in the manuscript), and in script number 5 we compare all the different approaches.

In the manuscript, we only work with bcftools inferred ROHs/FROH as all measures were highly positively correlated but plink appeared highly sensitive to alternative (realistic) parameter settings.
