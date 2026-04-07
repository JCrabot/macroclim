##                  Species distribution models                ##
#################################################################
#### Generating the presence table ####

## We join the species occurrences with the present environmental variables table

presence <- left_join(spdata_ids, stats_table_present, by = "HYRIV_ID")

if (discard != "no"){ # discard polluted sites
  presence <- presence %>%
    filter((HYRIV_ID %in% list_less_polluted_reaches$HYRIV_ID))
}

#### Generating the absence table ####

# to get the id of the river reaches where the taxon was not observed
absence_ids <- sites %>%
  filter(site_id %in% (sites %>% filter(! (site_id %in% spdata$site_id)))$site_id)

# to remove "false absences" from datasets where the taxon was not identified at all
datasets_to_use <-  names(list_datasets %>% filter(genus == taxa_selected) %>% select(-c(genus, sum))
      )[((list_datasets %>% filter(genus == taxa_selected)%>% select(-c(genus, sum)))> 0)]

absence_ids <- absence_ids %>%
  filter(site_id %in% (sites %>% filter(dataset %in% datasets_to_use))$site_id)

absence <- left_join(absence_ids, stats_table_present, by = "HYRIV_ID")

if (discard != "no"){ # discard polluted sites
  absence <- absence %>%
    filter((HYRIV_ID %in% list_less_polluted_reaches$HYRIV_ID)) 
}

#### Combining the training data and compute SDMs ####

# Preparing data to train the model on
data_model_RA <- data.table::rbindlist(list(presence, absence), fill = TRUE)

presence$occurrence <- 1
absence$occurrence <- 0
data_model_RA <- data.table::rbindlist(list(presence, absence), fill = TRUE)
data_model_RA$occurrence <- as.factor(data_model_RA$occurrence)

data_train_RA <- data_model_RA

if (nrow(presence)> 10){
  #### SDM species distribution model =  random forest ####
  
  set.seed(17)
  
  num.trees = 1000
  mtry = 7
  
  rf_param <- rf_parameters[which(rf_parameters$taxon == taxa_selected),]
  min.node.size = rf_param[1, "min.node.size"]
  min.bucket = rf_param[1, "min.bucket"]
  max.depth = rf_param[1,"max.depth"]
  
  if (downsample) {   # with down-sampling
    pres_num <- as.numeric(table(data_train_RA$occurrence)["1"])
    sample_size <- c("0" = pres_num / nrow(data_train_RA),
                     "1" = pres_num / nrow(data_train_RA))
  } else {  # without down-sampling
    sample_size <- 1
  }

  for (random in c(TRUE,FALSE)){
    if (random == TRUE) {
      
      # Generate a random variable
      data_train_RA$random <- runif(dim(data_train_RA)[1], min=0, max=100)
      
      # Run the random forest model
      model_RA_random <- ranger(data_train_RA$occurrence ~ .,
                         data = data_train_RA %>% select(wt_min_mean:yr_irclass, random), # with random factor
                         num.trees = num.trees, 
                         mtry= mtry,
                         min.node.size = min.node.size,
                         min.bucket = min.bucket,
                         max.depth = max.depth,
                         sample.fraction = sample_size,
                         replace = T,
                         oob.error = T,
                         keep.inbag = T,
                         num.threads = 4,
                         importance ='impurity',
                         probability = proba
      )
      
      varimp_RA_random <- importance(model_RA_random)[order(importance(model_RA_random),decreasing=TRUE)] # variables in decreasing order of importance in the model
      random_index <- which(names(varimp_RA_random) == "random")
      if (random_index > mtry & random_index){
        n_lim = min((random_index+3), length(varimp_RA_random))
        selec_var <- names(varimp_RA_random)[1:n_lim] #tolerance of 3 ranks to include the variability between models 
        selec_var <- selec_var[selec_var != "random"]
      } else {
      selec_var <- names(varimp_RA_random)[1:(mtry+1)]
      selec_var <- selec_var[selec_var != "random"]
      }
    
    }  else {
      
      if (selec){ # considers only environmental variables more important than the random factor
        model_RA <- ranger(data_train_RA$occurrence ~ .,
                           data = data_train_RA %>% select(any_of(selec_var)), #without random factor, keeping variables > random
                           num.trees = num.trees, 
                           mtry= mtry,
                           min.node.size = min.node.size,
                           replace = T,
                           sample.fraction = sample_size,
                           oob.error = T,
                           keep.inbag = T,
                           num.threads = 4,
                           importance ='impurity',
                           probability = proba
        )
      } else { # considers all environmental variables
        model_RA <- ranger(data_train_RA$occurrence ~ .,
                           data = data_train_RA %>% select(wt_min_mean:yr_irclass), # without random factor
                           num.trees = num.trees, 
                           mtry= mtry,
                           min.node.size = min.node.size,
                           replace = T,
                           sample.fraction = sample_size,
                           oob.error = T,
                           keep.inbag = T,
                           num.threads = 4,
                           importance ='impurity',
                           probability = proba
        )
      }
      
    }
  }
  
  ## Inspect the model
  # model_RA
  
  #### Predictions ----
if (random == FALSE) {  # because cannot be computed when the random variable has been generated
  #### Model prediction for the present RA ####
  
  current_pred_RA <- predict(model_RA, data = stats_table_present[,!1])
  
    name_val <- "pred_occur" # necessary if the raster is exported afterwards
  
    if (proba) {
      # probability of occurrence
  
      current_hab_RA <- data.table(HYRIV_ID = stats_table_present$HYRIV_ID,
                                pred_occur = as.integer((round(current_pred_RA$predictions[,2],2) * 100)))
      pred_name <-paste0(wd, "gis/predictions/", taxa_selected, "_", discard, "discard_current_occur_proba.shp") # if raster exported afterwards
  
    } else {
      # binary response present 1 or absent 0
  
      current_hab_RA <- data.table(HYRIV_ID = stats_table_present$HYRIV_ID,
                                   pred_occur = as.integer(as.character(current_pred_RA$predictions)))
      pred_name <- paste0(wd, "gis/predictions/", taxa_selected, "_", discard, "discard_current_occur_01.shp")
    }
  
  # write.csv(current_hab_RA, pred_name, row.names = FALSE)
  
  if (export) {
    
    if (! ("rivers" %in% ls(envir = .GlobalEnv))) {
      rivers <-terra::vect(paste0(wd,"hydroatlas/RiverATLAS_v10_shp/RiverATLAS_v10_eu.shp"),
                                   extent=bbox_terra) #already cropped to the right extent
    }
    # create shapefile with predictions
    river_pred <- merge(rivers,current_hab_RA, by = "HYRIV_ID")
    # plot(river_pred, y  ="pred_occur")
    writeVector(river_pred,
                filename = pred_name,
                overwrite=TRUE)
    rm(river_pred)
  }
  
  #### Model prediction_RA for the future RA ####
  
  #  Now we use the model with future variables
  pred_RA <- predict(model_RA, data = stats_table_future[,!1])
  
    name_val = "pred_occur"
  
    if (proba) {
      # probability of occurrence
  
    prediction_RA <- data.table(HYRIV_ID = stats_table_future$HYRIV_ID,
                                pred_occur = as.integer((round(pred_RA$predictions[,2],2) * 100)))
  
    pred_name = paste0(wd, "gis/predictions/", taxa_selected, "_", discard, "discard_future_occur_proba.shp")
  
      }  else {
        prediction_RA <- data.table(HYRIV_ID = stats_table_future$HYRIV_ID,
                                    pred_occur = as.integer(as.character(pred_RA$predictions)))
  
        pred_name = paste0(wd, "gis/predictions/", taxa_selected, "_", discard, "discard_future_occur_01.shp")
  
    }

  if (export) {
    # create shapefile with predictions
    river_pred <- merge(rivers,prediction_RA, by = "HYRIV_ID")
    writeVector(river_pred,
                filename = pred_name,
                overwrite=TRUE)
    rm(river_pred)
  }
  
  
  #### Map difference in distribution RA ####
  
  
  #adapt file name to the model parameters

    prediction_RA$diff=prediction_RA$pred_occur-current_hab_RA$pred_occur
    if (proba) {
      diff_name <- paste0(wd, "gis/predictions/", taxa_selected, "_", discard, "discard_diff_occur_proba.shp")
    } else {
      diff_name <- paste0(wd, "gis/predictions/", taxa_selected, "_", discard, "discard_diff_occur_01.shp")
    }
  
  # print(paste("min: ",min(prediction_RA$diff), " max: ", max(prediction_RA$diff)))
    if (!dir.exists(paste0(wd, "output/predictions/", prediction_folder))){
      dir.create(paste0(wd, "output/predictions/", prediction_folder), recursive=T)
    }
    write.csv(prediction_RA, paste0(wd, "output/predictions/", prediction_folder, taxa_selected, "_", discard, "discard_diff.csv"),
            row.names = FALSE)
  
  if (export) {
    # create shapefile with predictions
    river_pred <- merge(rivers,prediction_RA[,c(1,3)], by = "HYRIV_ID")
    writeVector(river_pred,
                filename = diff_name,
                overwrite=TRUE)
    rm(river_pred)
  }
  

}

}