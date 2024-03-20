
## function to calculate load
calculate_load_snpef <- function(file){
  # only get GT info, PL and DP are filtered by already anyway 
  gt <- c(10:(ncol(file)))
  select_n3 <- function(x){x = substr(x,1,3)}
  file[gt] <- lapply(file[gt], select_n3)
  
  # replace genotype with RS value but separate per zygosity, do per ID
  load <- list()
  for( id in 10:(ncol(file))){
    subset_id <- file[,c(1:9, id)]
    potential_data <- subset(subset_id, subset_id[[10]] == "1/0" | subset_id[[10]] == "0/1")
    realized_data <- subset(subset_id, subset_id[[10]] == "1/1")
    total_data <- subset(subset_id, subset_id[[10]] == "1/1" | subset_id[[10]] == "1/0" | subset_id[[10]] == "0/1")
    potential_load_sum <- nrow(potential_data)
    realized_load_sum <- nrow(realized_data)
    total_load_sum <- nrow(total_data)
    n_genotyped <- nrow(subset_id) - nrow(subset(subset_id, subset_id[[10]] == "./."))
    n_total <- nrow(subset_id)
    df <- data.frame(id = colnames(file[id]),
                     n_total = n_total,
                     n_genotyped = n_genotyped,
                     L_masked = potential_load_sum / n_genotyped,
                     L_expressed = realized_load_sum / n_genotyped,
                     Lt_additive = (potential_load_sum*0.5 + realized_load_sum) / n_genotyped)
    load[[id]] <- df
  }
  load <- do.call(rbind.data.frame, load)
  return(load)
  
}

