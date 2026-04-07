##          Model output exploration and summary table   FOR NON FLYING-TAXA      ##
if (!file.exists(paste0(wd,"output/taxa_selection_parameters_", file_name_frag, ".csv"))){
  
  niche_evol <- read.csv(file=paste0(wd,"output/niche/", niche_folder, "niche_evol.csv"), header = TRUE, sep = ";")
  
  if (!exists("stats_table_present")) { # same table for all species
    
    # Read the previously created file
    stats_table_zon <- fread(paste0(wd, "environment/environment.csv"))
    
    # We will keep only the mean of each variable of the stats_table,
    # to use as predictors in the species distribution model.
    stats_table_zon <- stats_table_zon %>%
      dplyr::select(-ends_with("_sd") & -ends_with("_min") & -ends_with("_max"))
    
    stats_table <- stats_table_zon
    
    # Some climatic variables need re-scaling before the modelling.
    # We define the following function:
    offset <- function(x, na.rm = F) (x - 273.15)
    
    # … and apply them to rescale the values
    
    stats_table <- stats_table  %>%
      mutate(across(starts_with("WT-c"), offset)) %>%
      mutate(across(starts_with("WT-d"), offset)) %>%
      mutate(across(starts_with("WT-h"), offset)) %>%
      mutate(across(starts_with("WT-m"), offset)) %>%
      mutate(across(starts_with("WT-wm_"), offset)) %>%
      mutate(across(starts_with("WT-wq"), offset))
    
    # We split the dataset into two datasets, according to present and future climatic variables.
    # The first will be used in the training of the model and the second in the projection of future suitable habitats
    
    stats_table_present <- stats_table %>%
      dplyr::select(!contains("204")) %>% # exclude all future variables
      dplyr::select(!contains("future")) %>%
      rename_with(~str_remove(., '_present')) %>%
      rename_with(~str_remove(., '_1981.2010')) %>%
      rename_with(~str_remove(., 'LUCAS_LUC_v1.1_historical_Europe_0.1deg_2010_2015_')) %>%
      rename_with(~str_remove(., '_ipsl_hist_1976-01-31_to_2005-12-31'))
    
    stats_table_future <- stats_table %>%
      rename_with(~str_remove(., '_future')) %>%
      dplyr::select(!contains("_19") & !contains("_2010")) %>% # exclude all present variables
      rename_with(~str_remove(., '_ipsl-cm6a-lr_ssp126_2041-2070')) %>%
      rename_with(~str_remove(., 'LUCAS_LUC_v1.1_ssp126_Europe_0.1deg_2046_2055_')) %>%
      rename_with(~str_remove(., '_ipsl_rcp2p6_2041-01-31_to_2060-12-31'))
    
    stats_table_present <- janitor::clean_names(stats_table_present) 
    stats_table_future <- janitor::clean_names(stats_table_future) 
    
    colnames(stats_table_present)[1] <- "HYRIV_ID"
    colnames(stats_table_future)[1] <- "HYRIV_ID"
    
    var_final = c("HYRIV_ID", "wt_min_mean", "wt_wm_mean", "wt_range_mean", "wt_ztw_mean",
                  "q_dq_mean", "q_zfw_mean", "q_si_mean",
                  "bio04_mean","bio15_mean","bio18_mean",
                  "lc3_mean", "lc4_mean", "lc5_mean", "lc7_mean", "lc8_mean", "lc9_mean",
                  "lc14_mean", "lc15_mean",
                  "gla_pc_use", "pop_ct_csu", "hft_ix_c09", "ria_ha_csu",
                  "riv_tc_usu", "gwt_cm_cav", "slt_pc_cav", "snd_pc_uav", "ele_mt_cmn",
                  "aet_mm_cyr", "ero_kh_cav",
                  "yr_irclass"
    )
    
    stats_table_present <- stats_table_present[,..var_final]
    stats_table_future <- stats_table_future[,..var_final]
    
    # Clear up memory
    rm(list=c("stats_table_zon","stats_table")) ; gc()
    
  }
  
  #### Selecting non-flying taxa ####
  list_taxa_nf <- unique(bio$taxon[bio$taxon %in% (taxo %>% filter(class !="Insecta"))$taxon])
  list_taxa_nf <- list_taxa[order(list_taxa_nf)]
  
  #### Compute new SELECTION PARAMETER IF changes has been brought to predictions ####
  if(!file.exists(paste0(wd,"output/niche/", niche_folder, "distrib_summary_frag.csv"))){
    log_time("Building distribution change table for fragmentation scenario")
    for (taxa in list_taxa_nf){
      if(file.exists(paste0(wd,"output/predictions/", prediction_folder, taxa,"_", discard, "discard_diff_frag.csv"))){
  
        pred <- read.csv(paste0(wd,"output/predictions/", prediction_folder, taxa,"_", discard, "discard_diff_frag.csv"), header=T, sep=",", stringsAsFactors = FALSE) %>%
          mutate(future = present + diff)
        
        if(nrow(pred %>% filter(!is.na(change)))>0){
          
          # print(taxa)
      
          #### Step 1 - Calculate new niche evol####
          
          species.niche.present <- stats_table_present %>%
            as.data.frame() %>%
            filter (HYRIV_ID %in% (pred$HYRIV_ID[pred$present==1])) %>%
            select(where(~ n_distinct(.) > 1)) %>%
            mutate(across(where(is.numeric), function(x) ifelse(x <= - 999999, NA, x))) %>%
            drop_na %>%
            mutate(Period = "present")
          
          species.niche.future <- stats_table_present %>%
            as.data.frame() %>%
            filter (HYRIV_ID %in% (pred$HYRIV_ID[pred$future==1])) %>%
            select(where(~ n_distinct(.) > 1)) %>%
            mutate(across(where(is.numeric), function(x) ifelse(x <= - 999999, NA, x))) %>%
            drop_na %>%
            mutate(Period = "future")
          
          if ((nrow(species.niche.present) > 1) & (nrow(species.niche.future) > 1)){
            species.niche <- rbind(species.niche.present %>% select(any_of(colnames(species.niche.future))),
                                   species.niche.future %>% select(any_of(colnames(species.niche.present)))) %>%
              select(Period, everything()) %>% #reorder 
              mutate(Period = factor(Period, levels=c("present","future")))
            
    
          write.table(species.niche,
                      file=paste0(wd, "output/niche/", niche_folder, "niche_evol_", taxa, "_frag.csv"),
                      row.names= FALSE, sep=";")
          
          #### Step 2: exploring spatial distribution of each taxon ####
          
          if (!file.exists(paste0(wd,"output/niche/", niche_folder, "distrib_summary_frag.csv"))){
            distrib_summary <- data.frame(matrix(nrow = 1, ncol = 31)) 
            colnames(distrib_summary) <- c("taxon",
                                           "nb_reaches_pres", "nb_reaches_fut",
                                           "pc_reaches_added","pc_reaches_lost", "size_area_change", "change_continuity",
                                           
                                           "latitude_mean_present", "latitude_min_present","latitude_max_present",  "latitude_range_present",
                                           "altitude_mean_present", "altitude_min_present","altitude_max_present",  "altitude_range_present",
                                           
                                           "latitude_mean_future", "latitude_min_future","latitude_max_future",  "latitude_range_future",
                                           "altitude_mean_future", "altitude_min_future","altitude_max_future",  "altitude_range_future",
                                           
                                           "latitude_mean_change", "latitude_min_change","latitude_max_change",  "latitude_range_change",
                                           "altitude_mean_change", "altitude_min_change","altitude_max_change",  "altitude_range_change"
                                           
            )
            
            write.table(distrib_summary,
                        file=paste0(wd,"output/niche/", niche_folder, "distrib_summary_frag.csv"),
                        row.names= FALSE, sep=";", dec=".")
            
          }
          
          distrib_summary <- read.csv(paste0(wd,"output/niche/", niche_folder, "distrib_summary_frag.csv"), header = TRUE, sep = ";", dec=".")
          
          if (!exists("geom_rivers", envir = globalenv())) {
            river_segments <-terra::vect(paste0(wd,"/environment/hydroatlas/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp"),
                                         extent=bbox_terra) #already cropped to the right extent
            geom_rivers <- terra::geom(river_segments)%>%
              as.data.frame()%>%
              select(geom, x, y)%>%
              group_by(geom) %>%
              summarise(x = mean(x, na.rm = TRUE),
                        y = mean(y, na.rm = TRUE)) %>%
              bind_cols(river_segments$HYRIV_ID) %>% rename (HYRIV_ID = 4)
            
            river_segments <- river_segments %>%
              as.data.frame()
          }
        
          nb_reaches_pres <- nrow(species.niche.present)
          nb_reaches_fut <- nrow(species.niche.future)
          pc_reaches_added <- nrow(species.niche.future %>% filter(!(HYRIV_ID %in% species.niche.present$HYRIV_ID)))/nb_reaches_pres*100
          pc_reaches_lost <- nrow(species.niche.present %>% filter(!(HYRIV_ID %in% species.niche.future$HYRIV_ID)))/nb_reaches_pres*100
          size_area_change <- (nb_reaches_fut - nb_reaches_pres)/nb_reaches_pres*100
          
          latitude_mean_present <- mean((geom_rivers %>% filter(HYRIV_ID %in% species.niche.present$HYRIV_ID))$y, na.rm=T)
          latitude_min_present <- min((geom_rivers %>% filter(HYRIV_ID %in% species.niche.present$HYRIV_ID))$y, na.rm=T)
          latitude_max_present <- max((geom_rivers %>% filter(HYRIV_ID %in% species.niche.present$HYRIV_ID))$y, na.rm=T)
          latitude_range_present <- latitude_max_present - latitude_min_present
          altitude_mean_present <- mean((river_segments %>% filter(HYRIV_ID %in% species.niche.present$HYRIV_ID))$ele_mt_cav, na.rm=T)
          altitude_min_present <- min((river_segments %>% filter(HYRIV_ID %in% species.niche.present$HYRIV_ID))$ele_mt_cav, na.rm=T)
          altitude_max_present <- max((river_segments %>% filter(HYRIV_ID %in% species.niche.present$HYRIV_ID))$ele_mt_cav, na.rm=T)
          altitude_range_present <- altitude_max_present - altitude_min_present
          
          latitude_mean_future <- mean((geom_rivers %>% filter(HYRIV_ID %in% species.niche.future$HYRIV_ID))$y, na.rm=T)
          latitude_min_future <- min((geom_rivers %>% filter(HYRIV_ID %in% species.niche.future$HYRIV_ID))$y, na.rm=T)
          latitude_max_future <- max((geom_rivers %>% filter(HYRIV_ID %in% species.niche.future$HYRIV_ID))$y, na.rm=T)
          latitude_range_future <- latitude_max_future - latitude_min_future
          altitude_mean_future <- mean((river_segments %>% filter(HYRIV_ID %in% species.niche.future$HYRIV_ID))$ele_mt_cav, na.rm=T)
          altitude_min_future <- min((river_segments %>% filter(HYRIV_ID %in% species.niche.future$HYRIV_ID))$ele_mt_cav, na.rm=T)
          altitude_max_future <- max((river_segments %>% filter(HYRIV_ID %in% species.niche.future$HYRIV_ID))$ele_mt_cav, na.rm=T)
          altitude_range_future <- altitude_max_future - altitude_min_future
          
          latitude_mean_change <- latitude_mean_future - latitude_mean_present
          latitude_min_change <- latitude_min_future - latitude_min_present
          latitude_max_change <- latitude_max_future - latitude_max_present
          latitude_range_change <- latitude_range_future - latitude_range_present
          altitude_mean_change <- altitude_mean_future - altitude_mean_present
          altitude_min_change <- altitude_min_future - altitude_min_present
          altitude_max_change <- altitude_max_future - altitude_max_present
          altitude_range_change <- altitude_range_future - altitude_range_present
          
          # looking for continuity of the distribution area
          species.hybas.present <- species.niche.present %>%
            left_join(hyriv_hybas, by = "HYRIV_ID") %>%
            select(HYBAS_L12) %>%
            distinct()
          basin_present <- basins12 %>%
            filter(HYBAS_ID %in% species.hybas.present$HYBAS_L12)
          basin_present <- aggregate(basin_present, by = "MAIN_BAS", dissolve = TRUE) |> disagg()
          basin_present$area <- expanse(basin_present, unit="km", transform=TRUE)
          
          species.hybas.future <- species.niche.future %>%
            filter(Period == "future") %>%
            left_join(hyriv_hybas, by = "HYRIV_ID") %>%
            select(HYBAS_L12) %>%
            distinct()
          basin_future <- basins12 %>%
            filter(HYBAS_ID %in% species.hybas.future$HYBAS_L12)
          basin_future <- aggregate(basin_future, by = "MAIN_BAS", dissolve = TRUE) |> disagg()
          basin_future$area <- expanse(basin_future, unit="km", transform=TRUE)
          
          continuous_area_present <- mean(basin_present$area)
          continuous_area_future <- mean(basin_future$area)
          change_continuity <- (continuous_area_future - continuous_area_present)/ continuous_area_present * 100
          
          distrib_summary <- rbind(distrib_summary,
                                   c(taxa,
                                     nb_reaches_pres, nb_reaches_fut,
                                     pc_reaches_added,pc_reaches_lost, size_area_change, change_continuity,
                                     
                                     latitude_mean_present, latitude_min_present,latitude_max_present,  latitude_range_present,
                                     altitude_mean_present, altitude_min_present,altitude_max_present,  altitude_range_present,
                                     
                                     latitude_mean_future, latitude_min_future,latitude_max_future,  latitude_range_future,
                                     altitude_mean_future, altitude_min_future,altitude_max_future,  altitude_range_future,
                                     
                                     latitude_mean_change, latitude_min_change,latitude_max_change,  latitude_range_change,
                                     altitude_mean_change, altitude_min_change,altitude_max_change,  altitude_range_change
                                   )
          )
          
          write.table(distrib_summary, file=paste0(wd,"output/niche/", niche_folder, "distrib_summary_frag.csv"),
                      row.names= FALSE, sep=";", dec=".")
        
          
          }  else if ((nrow(species.niche.present) > 1) & (nrow(species.niche.future) ==0 )){
            
            if (!file.exists(paste0(wd,"output/niche/", niche_folder , "disappearing_taxa.csv"))){
              disappearing <- data.frame(matrix(nrow = 1, ncol = 2))
              colnames(disappearing) <- c("taxon", "nb_reaches_present")
              
              write.table(disappearing,  file=paste0(wd,"output/niche/", niche_folder , "disappearing_taxa.csv"),
                          row.names= FALSE, sep=";", dec=".")
              
            }
            
            disappearing <- read.csv(paste0(wd,"output/niche/", niche_folder , "disappearing_taxa.csv"), header = TRUE, sep = ";", dec=".")
            
            disappearing <- rbind(disappearing, c(taxa, nrow(species.niche.present)))
            
            write.table(disappearing,file=paste0(wd,"output/niche/", niche_folder , "disappearing_taxa.csv"),row.names= FALSE, sep=";", dec=".")
            
          }
        }
      }
    }
  }
  #### Creating selection parameter table only for non-insect taxa ####

  if(!file.exists(paste0(wd,"output/taxa_selection_parameters_", file_name_frag ,".csv"))){
    #### Niche evol ####
  
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
  
  niche_evol_frag <- data.frame(matrix(nrow=0, ncol = length(var_final)+2))
  colnames(niche_evol_frag) <- c("taxon","Period", var_final)
  
  if(!file.exists(paste0(wd, "/output/niche/", niche_folder, "niche_evol_frag.csv"))){
    
    list_files_nicheevol <- list.files(paste0(wd,"/output/niche/", niche_folder), pattern = "_frag.csv$", full.names = TRUE)
    list_files_nicheevol_selec <- list_files_nicheevol[str_remove_all(basename(list_files_nicheevol), "niche_evol_|_frag.csv") %in% bio$taxon]
    
    for (i in list_files_nicheevol_selec) {
      # Loads files previously created
      niche_species <- read.csv(i, header = TRUE, sep = ";")
      niche_species$taxon <- str_remove(str_remove(basename(i), "niche_evol_"),"_frag.csv")
      
      # print(i)
      
      niche_species <- niche_species %>%
        select(- HYRIV_ID) %>% 
        group_by(taxon, Period)%>%
        summarise_if(is.numeric, list(mean = mean, min = min, max = max, sd = sd), na.rm = TRUE)
      
      # Combine in one big file
      niche_evol_frag <- rbind(niche_evol_frag,niche_species)
      rm(niche_species)
    }
    
    write.table(niche_evol_frag, file=paste0(wd, "output/niche/",niche_folder,"niche_evol_frag.csv"), row.names= FALSE, sep=";")
  } else {
    niche_evol_frag <- read.csv(file=paste0(wd, "output/niche/", niche_folder, "niche_evol_frag.csv"), header = TRUE, sep = ";")
    }
  
  
  niche_evol_frag <- gather(niche_evol_frag, variable, value, 3:ncol(niche_evol_frag),
                       factor_key=TRUE) %>%
      arrange(value) %>%
      mutate(Period = as.factor(Period)) %>%
      spread(Period, value) %>%
      mutate(diff = future - present)%>%
      pivot_longer(cols = future:diff, names_to = "period") %>%
      mutate(variable = str_replace(variable, "_mean_","_")) %>%
      filter(!((period == "diff" & grepl("_sd",variable))|
                 (period == "future"))) %>%
      mutate(variable = paste(variable, period, sep = "_")) %>%
      select(-period) %>%
      filter(!(str_detect(variable,"^lc|snw|gla|pop|hft|ria|riv|gwt|slt|snd|ero"))) %>%
      spread(variable, value) #%>%
      # select(where(~ n_distinct(.) > 1))
    
    #### Load variable importance ####
    variable_importance <- read.csv(paste0(wd,"output/niche/", niche_folder, "variable_importance.csv"),
                                    header = TRUE, sep = ";", dec= ".")
    #### Load summary on spatial distribution ####
    distrib_summary <- read.csv(paste0(wd,"output/niche/", niche_folder, "distrib_summary_frag.csv"),
                                header = TRUE, sep = ";", dec= ".")
  
    #### Combine interesting parameters in one summary table ####
    taxa_selection_parameters_frag <-  niche_evol_frag %>%
      select(taxon, wt_min_mean_present, wt_range_mean_present, wt_min_mean_diff, wt_range_mean_diff) %>%
      left_join(variable_importance %>%
                  select(taxon, OOB, wt_min_mean, wt_wm_mean, wt_range_mean, wt_ztw_mean, ele_mt_cmn),
                by = "taxon") %>%
      left_join(distrib_summary %>%
                  select(taxon, pc_reaches_added, pc_reaches_lost, size_area_change, change_continuity, contains("change")),
                by = "taxon") %>%
      left_join(summary_occurences %>%
                  select(taxon, nb_sites),
                by = "taxon") %>%
      distinct() %>%
      mutate(latitude_max_change = as.numeric(latitude_max_change),
             latitude_range_change = as.numeric(latitude_range_change),
             altitude_max_change = as.numeric(altitude_max_change),
             altitude_range_change = as.numeric(altitude_range_change))
  
    write.table(taxa_selection_parameters_frag,
                file=paste0(wd,"output/taxa_selection_parameters_", file_name_frag ,"_aquatic_only.csv"),
                row.names= FALSE, sep=";", dec=".")
  } else {
    log_time("Loading distribution change table for fragmentation")
    taxa_selection_parameters_frag <- read.csv(paste0(wd,"output/taxa_selection_parameters_", file_name_frag ,"_aquatic_only.csv"),
                                               header = TRUE, sep = ";", dec= ".")
  }
  #### Modifying previous table ####
  taxa_selection_parameters <- read.csv(paste0(wd,"output/taxa_selection_parameters_", file_name, ".csv"),
                                          header = TRUE, sep = ";", dec= ".")
  
  taxa_selection_parameters <- taxa_selection_parameters %>%
    filter(! taxon %in% taxa_selection_parameters_frag$taxon) %>%
    bind_rows(taxa_selection_parameters_frag) %>%
    arrange(taxon)
  
  write.table(taxa_selection_parameters,
              file=paste0(wd,"output/taxa_selection_parameters_", file_name_frag ,".csv"),
              row.names= FALSE, sep=";", dec=".")
  
  rm(list=setdiff(ls(), c("wd","wd_script","env", "bio", "sites", "rivers", "river_segments", "bbox_terra", # raw data to keep
                          "geom_rivers","list_less_polluted_reaches", "list_taxa", "taxo", 
                          "log_time", "rescale", #functions
                          "list_datasets", "summary_occurences", # to check conditions to run script at genus level
                          "stats_table_present","stats_table_future", #environmental variables for SDMs
                          "rivers_to_keep", "hyriv_hybas", "basins12", "rivers_df", "rivers", "rivers_selec", # rivers data
                          "length", "dams_intersected", "dams_id",# used for fragmentation
                          "type", "proba", "export", "random", "discard", "fragmentation", "downsample", "selec", "rf_parameters", # parameters for functions
                          "prediction_folder", "niche_folder", "file_name", "file_name_frag" #output names
  ))) 
  gc()
  
}