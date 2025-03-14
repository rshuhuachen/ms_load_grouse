### load packages ####

pacman::p_load(brms, bayesplot, dplyr, data.table)

### load models ###

#### gerp ####
load(file = "output/models/annual/traits/model_attend_gerp45.RData")
fit_gerp_attend <- fit
load(file = "output/models/annual/traits/model_fight_gerp45.RData")
fit_gerp_fight <- fit
load(file = "output/models/annual/traits/model_dist_gerp45.RData")
fit_gerp_dist <- fit
load(file = "output/models/annual/traits/model_eyec_gerp45.RData")
fit_gerp_eyec <- fit
load(file = "output/models/annual/traits/model_blue_gerp45.RData")
fit_gerp_blue <- fit
load(file = "output/models/annual/traits/model_lyre_gerp45.RData")
fit_gerp_lyre <- fit

load(file = "output/models/annual/ams/model_trait_ams_gerp45.RData")
fit_gerp_ams <- fit

rm(fit)

### snpeff 
load(file = "output/models/annual/traits/model_attend_high.RData")
fit_high_attend <- fit
load(file = "output/models/annual/traits/model_fight_high.RData")
fit_high_fight <- fit
load(file = "output/models/annual/traits/model_dist_high.RData")
fit_high_dist <- fit
load(file = "output/models/annual/traits/model_eyec_high.RData")
fit_high_eyec <- fit
load(file = "output/models/annual/traits/model_blue_high.RData")
fit_high_blue <- fit
load(file = "output/models/annual/traits/model_lyre_high.RData")
fit_high_lyre <- fit

load(file = "output/models/annual/ams/model_trait_ams_high.RData")
fit_high_ams <- fit

rm(fit)

### indirect effects loop #####
get_indirect <- function(mediator, method, trait_model, ams_model){
  treatment = "b_scaletotal_load"
  path1 <- as_draws_df(trait_model, variable =treatment)
  path1 <- path1$b_scaletotal_load
  
  path2 <- as_draws_df(ams_model, variable = mediator)
  path2 <- unlist(c(path2[,1]))
  
  indirect <- path1*path2
  
  direct <- as_draws_df(ams_model, variable =treatment)
  direct <- direct$b_scaletotal_load
  
  total <- indirect + direct
  
  effect_attend <- data.frame(treatment = treatment,
                              mediator = mediator,
                              method = method,
                              indirect_median = round(median(indirect), 2),
                              indirect_lower = round(quantile(indirect, probs = c(.025)), 2),
                              indirect_upper = round(quantile(indirect, probs = c(.975)), 2),
                              direct_median = round(median(direct), 2),
                              direct_lower = round(quantile(direct, probs = c(.025)), 2),
                              direct_upper = round(quantile(direct, probs = c(.975)), 2),
                              total_median = round(median(total), 2),
                              total_lower = round(quantile(total, probs = c(.025)), 2),
                              total_upper = round(quantile(total, probs = c(.975)), 2),
                              path1_median = round(median(path1), 2),
                              path1_lower = round(quantile(path1, probs = c(.025)), 2),
                              path1_upper = round(quantile(path1, probs = c(.975)), 2),
                              path2_median = round(median(path2), 2),
                              path2_lower = round(quantile(path2, probs = c(.025)), 2),
                              path2_upper = round(quantile(path2, probs = c(.975)), 2),
                              indirect_lower_80 = round(quantile(indirect, probs = c(.1)), 2),
                              indirect_upper_80 = round(quantile(indirect, probs = c(.9)), 2),
                              direct_lower_80 = round(quantile(direct, probs = c(.1)), 2),
                              direct_upper_80 = round(quantile(direct, probs = c(.9)), 2))
  
  return(effect_attend)
}

effects <- data.frame(rbind(get_indirect(mediator="b_scalelyre", method = "gerp45", 
                                         trait_model=fit_gerp_lyre, ams_model = fit_gerp_ams),
                            get_indirect(mediator="b_scaleeyec",  method = "gerp45", 
                                         trait_model=fit_gerp_eyec, ams_model = fit_gerp_ams),
                            get_indirect(mediator="b_scaleblue",  method = "gerp45", 
                                         trait_model=fit_gerp_blue, ams_model = fit_gerp_ams),
                            get_indirect(mediator="b_scaleattend",  method = "gerp45", 
                                         trait_model=fit_gerp_attend, ams_model = fit_gerp_ams),
                            get_indirect(mediator="b_scalefight",  method = "gerp45", 
                                         trait_model=fit_gerp_fight, ams_model = fit_gerp_ams),
                            get_indirect(mediator="b_scaledist",  method = "gerp45", 
                                         trait_model=fit_gerp_dist, ams_model = fit_gerp_ams),
                            get_indirect(mediator="b_scalelyre", method = "high", 
                                         trait_model=fit_high_lyre, ams_model = fit_high_ams),
                            get_indirect(mediator="b_scaleeyec",  method = "high", 
                                         trait_model=fit_high_eyec, ams_model = fit_high_ams),
                            get_indirect(mediator="b_scaleblue",  method = "high", 
                                         trait_model=fit_high_blue, ams_model = fit_high_ams),
                            get_indirect(mediator="b_scaleattend",  method = "high", 
                                         trait_model=fit_high_attend, ams_model = fit_high_ams),
                            get_indirect(mediator="b_scalefight",  method = "high", 
                                         trait_model=fit_high_fight, ams_model = fit_high_ams),
                            get_indirect(mediator="b_scaledist",  method = "high", 
                                         trait_model=fit_high_dist, ams_model = fit_high_ams)))


write.csv(effects, file = "output/models/annual/direct_indirect_summary.csv", quote=F, row.names = F)
