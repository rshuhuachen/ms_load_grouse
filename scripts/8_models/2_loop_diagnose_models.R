pacman::p_load(brms, bayesplot, tidyverse)
extrafont::loadfonts(device="all")

#adjust pattern here to what models you want to diagnose
#note that the plotting takes a while so this is better to do in smaller batches

output_ams <- list.files(path = "output/models/annual/ams",
                     pattern = "model*", full.names=T)
output_trait <- list.files(path = "output/models/annual/traits/",
                         pattern = "model*", full.names=T)
output_total <- list.files(path = "output/models/total_hom_het/",
                         pattern = "lms*", full.names=T)
output_region <- list.files(path = "output/models/per_gene_region/",
                         pattern = "lms*", full.names=T)

output <- c(output_ams, output_trait, output_total, output_region)

diagnose_summary <- list()
for (i in 1:length(output)){
  #load fit
  load(file = output[[i]])
  #get posteriors
  posterior <- as.array(fit)
  log_ps <- log_posterior(fit)
  nuts <- nuts_params(fit) #divergence
  #get only beta and sd
  betas <- variables(fit)[grep("b_", variables(fit))]
  sd <- variables(fit)[grep("sd_", variables(fit))]
  
  #global patterns in divergence
  diverge_beta <- mcmc_parcoord(posterior, np = nuts, pars= betas)
  diverge_sd <- mcmc_parcoord(posterior, np = nuts, pars= sd)
  
  #identify collinearity between parameters
  collin_beta <- mcmc_pairs(posterior, np = nuts, pars= betas)
  collin_sd <- mcmc_pairs(posterior, np = nuts, pars= sd)
  
  #traceplot
  trace_beta <- mcmc_trace(posterior, pars = betas, np = nuts)
  trace_sd <- mcmc_trace(posterior, pars = sd, np = nuts)
  
  #rhat
  rhat <- mcmc_rhat(brms::rhat(fit))
  
  #effective sample size
  neff <- mcmc_neff(neff_ratio(fit))
  
  #autocorrelation
  autocor_beta <- mcmc_acf(posterior, pars = betas)
  autocor_sd <- mcmc_acf(posterior, pars=sd)
  
  #quick glance results
  areas <- mcmc_areas(fit, pars=betas)
  
  #combine in list
  diagnosis <- list(diverge_beta = diverge_beta, 
                    diverge_sd = diverge_sd, 
                    collin_beta = collin_beta, 
                    collin_sd = collin_sd, 
                    trace_beta = trace_beta, 
                    trace_sd = trace_sd, 
                    rhat = rhat, 
                    neff = neff, 
                    autocor_beta = autocor_beta, 
                    autocor_sd = autocor_sd,
                    areas = areas)
  
  
  modelname <- sub(".*/", "", output[i]) 
  modelname <- sub(".RData", "", modelname)
  
  # save output
  pdf(file=paste0("output/models/diagnosis/", modelname, ".pdf"))
  diverge_beta
  diverge_sd
  collin_beta
  collin_sd
  trace_beta
  trace_sd
  rhat
  neff
  autocor_beta
  autocor_sd
  areas
  dev.off()
  
  # add to summary
  diagnose_summary[[modelname]] <- diagnosis
}

#Look at individual plots like this:
diagnose_summary$attend_froh$rhat

