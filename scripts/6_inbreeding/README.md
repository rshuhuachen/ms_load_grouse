# Calling ROHs

In this repository, we infer runs of homozygosity using two different approaches: using the `bcftools -roh` algorithm which uses MCMC, and the `plink --homozyg` function which uses a sliding window approach. In script number 3, the function defined in `function_run_rohcalling.R` is called which not only calculates ROH/FROH using plink, but also plots the inferred ROHs in a circos plot. This function was defined to allows easy comparison between different parameter setting combinations. Additionally, we calculate genome-wide sMLH (standardised multi-locus heterozygosity) in script 4, and in script number 5 we compare all the different approaches.

In the manuscript, we only work with bcftools inferred ROHs/FROH as all measures were highly positively correlated but plink appeared highly sensitive to other parameter settings.
