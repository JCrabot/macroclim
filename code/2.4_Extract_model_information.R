##               Extract information out of distribution models ##
if (exists("varimp_RA_random", envir = globalenv())) {
  ########  Writing to disk variable importance by taxon   ----
  if (!dir.exists(paste0(wd, "output/niche/", niche_folder))){
    dir.create(paste0(wd, "output/niche/", niche_folder), recursive = T)
  }
  
  if (!file.exists(paste0(wd,"output/niche/", niche_folder, "variable_importance.csv"))){
    variable_importance <- data.frame(matrix(nrow = 1,
                                             ncol = ncol(data_train_RA %>% select(wt_min_mean:yr_irclass, random))+2))
    colnames(variable_importance) <- c("taxon", "OOB", colnames(data_train_RA %>% select(wt_min_mean:yr_irclass, random)))
    write.table(variable_importance,
                file=paste0(wd,"output/niche/", niche_folder, "variable_importance.csv"),
                row.names= FALSE, sep=";")
  
  }
  
  variable_importance <- read.csv(paste0(wd,"output/niche/", niche_folder, "variable_importance.csv"), header = TRUE, sep = ";")
  
  vector_var<-vector(length=ncol(variable_importance))
  vector_var[1] <- taxa_selected
  for (i in 3:ncol(variable_importance)) {
    vector_var[i] <- match(colnames(variable_importance)[i], names(varimp_RA_random))
  }
  
  vector_var[2] <- model_RA_random$prediction.error
  variable_importance <- rbind(variable_importance, vector_var)
  
  write.table(variable_importance,
              file=paste0(wd,"output/niche/", niche_folder, "variable_importance.csv"),
              row.names= FALSE, sep=";")
  
  rm(vector_var)
  
  #### Step 1: creating table with all variables for present and future ####
    species.niche.present <- stats_table_present %>%
      as.data.frame() %>%
      filter (HYRIV_ID %in% (current_hab_RA$HYRIV_ID[current_hab_RA$pred_occur == 1])) %>%
      select(where(~ n_distinct(.) > 1)) %>%
      mutate(across(where(is.numeric), function(x) ifelse(x <= - 999999, NA, x))) %>%
      drop_na %>%
      mutate(Period = "present")
    
    species.niche.future <- stats_table_present %>%
      as.data.frame() %>%
      filter (HYRIV_ID %in% (prediction_RA$HYRIV_ID[prediction_RA$pred_occur == 1])) %>%
      select(where(~ n_distinct(.) > 1)) %>%
      mutate(across(where(is.numeric), function(x) ifelse(x <= - 999999, NA, x))) %>%
      drop_na %>%
      mutate(Period = "future")
    
    if ((nrow(species.niche.present) > 1) & (nrow(species.niche.future) > 1)){
        species.niche <- rbind(species.niche.present %>% select(any_of(colnames(species.niche.future))),
                               species.niche.future %>% select(any_of(colnames(species.niche.present)))) %>%
          select(Period, everything()) %>% #reorder 
          mutate(Period = factor(Period, levels = c("present","future")))
        
      write.table(species.niche,
                  file=paste0(wd, "output/niche/", niche_folder, "niche_evol_",
                              taxa_selected,".csv"), row.names= FALSE, sep=";")
    
      #### Step 2: exploring spatial distribution of each taxon ####
      
      if (!file.exists(paste0(wd,"output/niche/", niche_folder, "distrib_summary.csv"))){
        distrib_summary <- data.frame(matrix(nrow = 1, ncol = 31))
        colnames(distrib_summary) <- c("taxon",
                                     "nb_reaches_pres", "nb_reaches_fut",
                                     "pc_reaches_added","pc_reaches_lost", "size_area_change","change_continuity",
                                     
                                     "latitude_mean_present", "latitude_min_present","latitude_max_present",  "latitude_range_present",
                                     "altitude_mean_present", "altitude_min_present","altitude_max_present",  "altitude_range_present",
                                     
                                     "latitude_mean_future", "latitude_min_future","latitude_max_future",  "latitude_range_future",
                                     "altitude_mean_future", "altitude_min_future","altitude_max_future",  "altitude_range_future",
                                     
                                     "latitude_mean_change", "latitude_min_change","latitude_max_change",  "latitude_range_change",
                                     "altitude_mean_change", "altitude_min_change","altitude_max_change",  "altitude_range_change"
                                     
                                     )
        
        write.table(distrib_summary,
                    file=paste0(wd,"output/niche/", niche_folder, "distrib_summary.csv"),
                    row.names= FALSE, sep=";", dec=".")
        
      }
      
      distrib_summary <- read.csv(paste0(wd,"output/niche/", niche_folder, "distrib_summary.csv"), header = TRUE, sep = ";", dec=".")
      
      # build object  with geographic coordinates of river segments
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
                             c(taxa_selected,
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
      
      write.table(distrib_summary, file=paste0(wd,"output/niche/", niche_folder, "distrib_summary.csv"),
                  row.names= FALSE, sep=";", dec=".")
      
    } else if ((nrow(species.niche.present) > 1) & (nrow(species.niche.future) ==0 )){
      
      if (!file.exists(paste0(wd,"output/niche/", niche_folder, "disappearing_taxa.csv"))){
        disappearing <- data.frame(matrix(nrow = 1, ncol = 2))
        colnames(disappearing) <- c("taxon", "nb_reaches_present")
        
        write.table(disappearing,  file=paste0(wd,"output/niche/", niche_folder, "disappearing_taxa.csv"),
                    row.names= FALSE, sep=";", dec=".")
        
      }
      
      disappearing <- read.csv(paste0(wd,"output/niche/", niche_folder, "disappearing_taxa.csv"), header = TRUE, sep = ";", dec=".")
      
      disappearing <- rbind(disappearing, c(taxa_selected, nrow(species.niche.present)))
      
      write.table(disappearing,file=paste0(wd,"output/niche/", niche_folder, "disappearing_taxa.csv"),row.names= FALSE, sep=";", dec=".")
      
    }
    
}