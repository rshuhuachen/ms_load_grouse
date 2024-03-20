#### Function to run ROH, create a plot to summarise called ROHs and to visualise their distribution with a circosplot ####

run_roh <- function(vcf,
                    choose_homozyg_snp,
                    choose_homozyg_kb,
                    choose_homozyg_density,
                    choose_homozyg_gap,
                    choose_homozyg_het, #set to same as homozyg_snp for the default being unlimited
                    choose_homozyg_window_snp,
                    choose_homozyg_window_het, 
                    choose_homozyg_window_missing,
                    choose_homozyg_window_threshold, 
                    runname){
  pacman::p_load(detectRUNS, tidyverse, gridExtra,egg, BioCircos)
  
  ### Step preprocessing - Convert vcf to plink file
  system(paste0("/prj/blackgrouse/bin/mambaforge/envs/wgr_snake/bin/plink --vcf ", vcf, " --recode --out ", getwd(), " --double-id --allow-extra-chr"))
  
  ### Step 1 - Call ROH ###
  
  #vcffile_loc <- "P:\\Black Grouse PhD\\Projects\\ROHcalling\\rawSNPcalls"
  system(paste0("/prj/blackgrouse/bin/plink --file ", vcf,
                " --homozyg",
                " --homozyg-snp ", choose_homozyg_snp,
                " --homozyg-kb ", choose_homozyg_kb,
                " --homozyg-density ", choose_homozyg_density,
                " --homozyg-gap ", choose_homozyg_gap,
                " --homozyg-het ", choose_homozyg_het,
                " --homozyg-window-snp ", choose_homozyg_window_snp,
                " --homozyg-window-het ", choose_homozyg_window_het,
                " --homozyg-window-missing ", choose_homozyg_window_missing,
                " --homozyg-window-threshold ", choose_homozyg_window_threshold,
                " --allow-extra-chr",
                " --out ", getwd(), "/output/inbreeding/", runname))
  
  
  ### Step 2 - Load in ROH call and put in individual IDs ###
  run <- read.table(paste0(getwd(), "/output/inbreeding/", runname, ".hom"), header = TRUE)
  files_ids <- read.csv(paste0(getwd(), "/data/metadata/file_list_all_bgi_clean.csv"))
  run <- left_join(run, files_ids[,c(2,9)], by = c("FID" = "loc"))
  write.csv(run, paste0(getwd(), "/output/inbreeding/", runname, "_output.csv"), quote=F, row.names = F)
  
  ### Step 3 - calculate per ID: nroh, min, max, mean
  if(nrow(run)>0){
    theme_set(theme_classic() + theme(panel.background = element_rect(fill = "white")))
    
    roh_summary_perid_nroh <- run %>% 
      group_by(id) %>%
      summarise(n_roh = n()) 
    
    roh_summary_perid_mean_roh <- aggregate(run[, 9], #choose KB
                                            list(run$id), mean)
    names(roh_summary_perid_mean_roh)[2] <- "mean_roh_kb"
    
    roh_summary_perid_min_roh <- aggregate(run[, 9], 
                                           list(run$id), min)
    names(roh_summary_perid_min_roh)[2] <- "min_roh_kb"
    
    roh_summary_perid_max_roh <- aggregate(run[, 9], 
                                           list(run$id), max)
    names(roh_summary_perid_max_roh)[2] <- "max_roh_kb"
    
    roh_summary_perid <- as.data.frame(left_join(roh_summary_perid_nroh, roh_summary_perid_mean_roh, by = c("id" = "Group.1")) %>%
                                         left_join(roh_summary_perid_min_roh, by = c("id" = "Group.1")) %>%
                                         left_join(roh_summary_perid_max_roh, by = c("id" = "Group.1")) )
    
    ### Step 4 - calculate FROH
    #total genome length = 1004266063
    froh <- aggregate(run[,9], list(run$id), sum)
    names(froh) <- c("id", "sum")
    froh$sum <- froh$sum * 1000 #in bases not kb
    froh$Froh_genome <- froh$sum / 1004266063
    write.csv(froh, paste0(getwd(), "/output/inbreeding/", runname, "_froh.csv"), quote=F, row.names = F)
    # add to roh_summary_per_id
    roh_summary_perid <- left_join(roh_summary_perid, froh, by = "id")
    
    ### Step 5 - Visualize FROH distribution
    
    froh_hist <- ggplot(froh, aes(x = Froh_genome)) + 
      geom_density(colour = "black",fill = "#5E7085") + 
      labs(title = "FROH distribution", y = "Frequency", x = (bquote(F[ROH])), subtitle = 
             paste0("Minimum FROH = ", round(min(froh$Froh_genome), 4),
                    ", maximum FROH = ", round(max(froh$Froh_genome), 4),
                    ", mean FROH = ", round(mean(froh$Froh_genome), 4),
                    ", 
N IDs with ROH = ", length(unique(froh$id))))
    
    ### Step 6 - Visualize ROH length distribution
    # add a category 'ROH length' to run df
    run <- run %>% mutate(ROH_CAT = case_when(
      KB > 0 & KB < 1000 ~ "0-1 Mb",
      KB > 1000 & KB < 2000 ~ "1-2 Mb",
      KB > 2000 & KB < 3000 ~ "2-3 Mb",
      KB > 3000 & KB < 4000 ~ "3-4 Mb",
      KB > 4000 & KB < 5000 ~ "4-5 Mb",
      KB > 5000 & KB < 6000 ~ "5-6 Mb",
      KB > 6000 & KB < 7000 ~ "6-7 Mb",
      KB > 7000 & KB < 8000 ~ "7-8 Mb",
      KB > 8000 & KB < 9000 ~ "8-9 Mb",
      KB > 9000 ~ "> 9Mb"
    ))
    run$ROH_CAT <- factor(run$ROH_CAT, levels = c("0-1 Mb", "1-2 Mb", "2-3 Mb",
                                                  "3-4 Mb","4-5 Mb","5-6 Mb",
                                                  "6-7 Mb","7-8 Mb","8-9 Mb","> 9Mb"))
    
    roh_length_hist <- ggplot(run, aes(x=ROH_CAT)) + geom_histogram(fill = "#5E7085", stat = "count")+
      labs(title = "ROH frequency by length category", y = "Frequency") + theme(axis.title.x = element_blank())
    
    write.csv(roh_summary_perid, paste0(getwd(), "/output/inbreeding/", runname, "_summary_rohs.csv"), quote=F, row.names = F)
    
    ### Step 7 - Visualise ROH distribution with Circos
    library(data.table); library(BioCircos); library(grDevices)
    
    refgenome <- fread(paste0(getwd(), "/data/genomes/refgenome/tom-uni3242-mb-hirise-s8mkv__08-23-2022__hic_output.fasta.fai"))
    
    names(refgenome) <- c("contig", "bases", "start", "bases_per_line", "bytes_per_line")
    
    #select only scaffolds with at least one roh
    bg_select <- refgenome[refgenome$contig %in% run$CHR]
    bg_select$contig_short <- gsub("scaffold_", "", bg_select$contig)
    bg_select_circo <- setNames(as.list(bg_select$bases), c(bg_select$contig_short))
    run$chr <- gsub("scaffold_" , "", run$CHR)
    
    #set to real id names
    run <- run %>% group_by(id) %>%
      mutate(id_no = row_number())
    run$id_no = as.integer(factor(run$id))
    rohs_chromosomes = c(run$chr)
    rohs_begin = c(run$POS1)
    rohs_end = c(run$POS2)
    rohs_colours = data.frame(id = c(1:length(unique(run$id))),
                              col = c(rainbow(length(unique(run$id)))))
    rohs_colours$col <- as.character(rohs_colours$col)
    run_cols <- left_join(run, rohs_colours, by = c("id_no" = "id"))
    
    tracklist_perid = NULL
    
    for (i in 1:length(unique(run$id_no))){
      rohs_chromosomes_id = c(run$chr[which(run$id_no == i)])
      rohs_begin_id = c(run$POS1[which(run$id_no == i)])
      rohs_end_id = c(run$POS2[which(run$id_no == i)])
      tracklist_perid = tracklist_perid + 
        BioCircosArcTrack(paste0('track_', i),
                          rohs_chromosomes_id,
                          rohs_begin_id, 
                          rohs_end_id,
                          colors = c(rohs_colours[i,2]),
                          minRadius = 1.16 + i*0.02, 
                          maxRadius = 1.18 + i*0.02)
      
    }
    plot.new()
    circos <- BioCircos(tracklist = tracklist_perid, genome = bg_select_circo, 
                        genomeFillColor = "Oranges",
                        chrPad = 0.02, displayGenomeBorder = FALSE, yChr =  FALSE,
                        genomeTicksDisplay = FALSE,  genomeLabelTextSize = 0)
    
    
    ### Return the other plots
    rohplots <- egg::ggarrange(froh_hist, roh_length_hist, nrow = 2, ncol = 1,
                               top = paste0("Visualising ROH calls: distribution of FROH and ROH length"),
                               bottom = paste0("Parameter settings: 
                           homozygosity_window_snp = ",  choose_homozyg_window_snp,
                                               ", homozyg_window_het = ", choose_homozyg_window_het,
                                               ", 
homozyg_window_missing = ", choose_homozyg_window_missing,
                                               ", homozyg_window_threshold = ", choose_homozyg_window_threshold,
                                               ", homozyg_gap = ", choose_homozyg_gap,
                                               ", 
homozyg_kb = ", choose_homozyg_kb,
                                               ", homozyg_density = ", choose_homozyg_density,
                                               ", homozyg_snp = ", choose_homozyg_snp,
                                               ", homozyg_het = ", choose_homozyg_het))
    
    ggsave(rohplots, file = paste0(getwd(), "/output/inbreeding/", runname, "_rohplots.png"), width=1000, height=700, units=c("px"))
    library(htmltools)
    save_html(circos, file = paste0(getwd(), "/output/inbreeding/", runname, "_circos.html"))# not working
    ### Return Circos
    return(circos)
    
    
  }
  else(print("No ROHs detected, removing run files"))
  if(nrow(run)==0){
    system(paste0("rm ", getwd(), "/output/inbreeding/", runname, "*"))
  }
  system(paste0("rm ", getwd(), "/output/inbreeding/", runname, ".hom.summary")) #too big a size
  }

