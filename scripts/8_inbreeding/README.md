In this folder, we calculate the genomic inbreeding coefficient FROH. We used the software BCFtools to identify runs of homozygosity (ROHs).

In the first script, we run BCFtools within R, extract the ROHs and calculate FROH (1_run_bcftools_rohs)
In the second script, we test whether inbreeding differs among lek sites (2_test_lek_effect)
In the third script, we test for inbreeding depression using Bayesian modelling. The model set up is the same as the LMS models in the next series of scripts (in 9_models/lms). We also run model diagnostics and test for prior sensitivity.