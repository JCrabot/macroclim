##    BUILDING TABLE WITH SELECTION PARAMETER IF DOES NOT EXIST   ##

if (!dir.exists(paste0(wd, "output/"))){
  dir.create(paste0(wd, "output/"), recursive = T)
}

if(!file.exists(paste0(wd,"output/taxa_selection_parameters_", file_name, ".csv"))) {
  
  print("-- Building distribution change table --")
  
  #### Calculate the change in environmental conditions for each taxon ####
  
  var_final = c("HYRIV_ID", "wt_min_mean", "wt_wm_mean", "wt_range_mean", "wt_ztw_mean",
                "q_dq_mean", "q_zfw_mean", "q_si_mean",
                "bio04_mean","bio15_mean","bio18_mean",
                "lc3_mean", "lc4_mean", "lc5_mean", "lc7_mean", "lc8_mean", "lc9_mean",
                "lc14_mean", "lc15_mean",
                "gla_pc_use", "pop_ct_csu", "hft_ix_c09", "ria_ha_csu",
                "riv_tc_usu", "gwt_cm_cav", "slt_pc_cav", "snd_pc_uav", "ele_mt_cmn", "aet_mm_cyr",
                "ero_kh_cav",
                "yr_irclass"
  )
  
  niche_evol <- data.frame(matrix(nrow=0, ncol = length(var_final)+2))
  colnames(niche_evol) <- c("taxon","Period", var_final)
  
  
  if(!file.exists(paste0(wd, "output/niche/", niche_folder, "niche_evol.csv"))){
    list_files <- list.files(paste0(wd, "output/niche/", niche_folder), pattern= "^niche_evol_")
    list_files <- list_files[!grepl("frag", list_files)]
    for (i in list_files) {
       # Loads files previously created
       niche_species <- read.csv(file=paste0(wd,"output/niche/",niche_folder, i), header = TRUE, sep = ";")
       niche_species$taxon <- str_remove(str_remove(i, "niche_evol_"),".csv")
    
       print(i)
    
       niche_species <- niche_species %>%
         select(- HYRIV_ID) %>% 
         group_by(taxon, Period)%>%
         summarise_if(is.numeric, list(mean = mean, min = min, max = max, sd = sd), na.rm = TRUE)
    
       # Combine in one big file
       niche_evol <- rbind(niche_evol,niche_species)
       rm(niche_species)
    }
    
    write.table(niche_evol, file=paste0(wd, "output/niche/",niche_folder,"niche_evol.csv"), row.names= FALSE, sep=";")
  } else {
    niche_evol <- read.csv(file=paste0(wd, "output/niche/", niche_folder, "niche_evol.csv"), header = TRUE, sep = ";")
    }
  
  niche_evol <- gather(niche_evol, variable, value, 3:ncol(niche_evol),
                                   factor_key=TRUE) %>%
    arrange(value) %>%
    mutate(Period = as.factor(Period))

  ## difference calculation
  niche_evol <- niche_evol %>%
    spread(Period, value) %>%
    mutate(diff = future - present)
  
  niche_evol <- pivot_longer(niche_evol, cols = future:diff, names_to = "period") %>%
    mutate(variable = str_replace(variable, "_mean_","_")) %>%
    filter(!((period == "diff" & grepl("_sd",variable))|
    (period == "future"))) %>%
    mutate(variable = paste(variable, period, sep = "_")) %>%
    select(-period) %>%
    filter(!(str_detect(variable,"^lc|snw|gla|pop|hft|ria|riv|gwt|slt|snd|ero"))) %>%
    spread(variable, value) %>%
    select(where(~ n_distinct(.) > 1))
  
  #### Load variable importance ####
  variable_importance <- read.csv(paste0(wd,"output/niche/", niche_folder, "variable_importance.csv"),
                                  header = TRUE, sep = ";", dec= ".")
  
  #### Load summary on spatial distribution ####
  distrib_summary <- read.csv(paste0(wd,"output/niche/", niche_folder, "distrib_summary.csv"),
                                  header = TRUE, sep = ";", dec= ".")

  #### Combine interesting parameters in one summary table ####
  taxa_selection_parameters <-  niche_evol %>%
    select(taxon, wt_min_mean_present, wt_range_mean_present, wt_min_mean_diff, wt_range_mean_diff) %>%
    left_join(variable_importance %>%
                select(taxon, OOB, wt_min_mean, wt_wm_mean, wt_range_mean, wt_ztw_mean, ele_mt_cmn),
              by = "taxon") %>%
    left_join(distrib_summary %>%
                select(taxon, pc_reaches_added, pc_reaches_lost, size_area_change,  change_continuity, contains("change")),
              by = "taxon") %>%
    left_join(summary_occurences %>%
                select(taxon, nb_sites),
              by = "taxon") %>%
    distinct() %>%
    mutate(latitude_max_change = as.numeric(latitude_max_change),
           latitude_range_change = as.numeric(latitude_range_change),
           altitude_max_change = as.numeric(altitude_max_change),
           altitude_range_change = as.numeric(altitude_range_change))
  
    write.table(taxa_selection_parameters,
                file=paste0(wd,"output/taxa_selection_parameters_", file_name,".csv"),
                row.names= FALSE, sep=";", dec=".")

}

rm(list=setdiff(ls(), c("wd","wd_script","env", "bio", "sites", "rivers", "river_segments", "bbox_terra", # raw data to keep
                        "geom_rivers","list_less_polluted_reaches", "list_taxa", "taxo", 
                        "log_time", "rescale", #functions
                        "list_datasets", "summary_occurences", # to check conditions to run script at genus level
                        "stats_table_present","stats_table_future", #environmental variables for SDMs
                        "rivers_to_keep", "hyriv_hybas", "basins12", "rivers_df", "rivers", "rivers_selec",# rivers data
                        "length", "dams_intersected", "dams_id",# used for fragmentation
                        "type", "proba", "export", "random", "discard", "fragmentation", "downsample", "selec", "rf_parameters", # parameters for functions
                        "prediction_folder", "niche_folder", "file_name", "file_name_frag" #output names
))) 
gc()


