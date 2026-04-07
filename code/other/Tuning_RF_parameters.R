########### Parameter tuning for random forests ######

#### Function to generate the training datasets, compute the RF model, and 
#### compute indicators of the model performance
calculate_tss <- function(taxa_selected){
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
  
  #### Combining the training data ####
  
  # data to train the model on
  data_model_RA <- data.table::rbindlist(list(presence, absence), fill = TRUE)
  
  presence$occurrence <- 1
  absence$occurrence <- 0
  data_model_RA <- data.table::rbindlist(list(presence, absence), fill = TRUE)
  data_model_RA$occurrence <- as.factor(data_model_RA$occurrence)
  
  data_train_RA <- data_model_RA
  
  #### Running RF ####
    recap <- vector(mode = "list", length = 13)
    
    set.seed(17)
    num.trees = 1000
    mtry = 6
    proba = FALSE
    
    for( i in 1:3){
    for (min.node.size in c(1, 5, 50, 100)){
      for (min.bucket in c(1, 5, 10, 20)){
          for (max.depth in c(3, 5, 20, Inf)){
          #### SDM species distribution model =  random forest ####
  
          # Set sample fraction size
          pres_num <- as.numeric(table(data_train_RA$occurrence)["1"])
          sample_size <- c("0" = pres_num / nrow(data_train_RA),
                           "1" = pres_num / nrow(data_train_RA))
          
          # Run the random forest model
          model_RA <- ranger(data_train_RA$occurrence ~ .,
                             data = data_train_RA %>% select(wt_min_mean:yr_irclass),
                            num.trees = num.trees, 
                              mtry= mtry,
                              min.node.size = min.node.size,
                              sample.fraction = sample_size,
                              min.bucket = min.bucket,
                              max.depth = Inf,
                              replace = T,
                              oob.error = T,
                              keep.inbag = T,
                              num.threads = 4,
                              importance ='impurity',
                              probability = proba
          )
          
          #### Predictions ----
          current_pred_RA <- predict(model_RA, data = stats_table_present[,!1])
          current_hab_RA <- data.table(HYRIV_ID = stats_table_present$HYRIV_ID,
                                       pred_occur = as.integer(as.character(current_pred_RA$predictions)))
          
          #### Calculating True Skill Statistic and other indicators of performance -----
          # contingency table
          conting <- current_hab_RA %>%
            right_join(data_test %>% select(HYRIV_ID, occurrence), by = "HYRIV_ID")
          
          presence_presence <- sum(conting$pred_occur == 1 & conting$occurrence == 1)
          presence_absence <- sum(conting$pred_occur == 1 & conting$occurrence == 0)
          absence_presence <- sum(conting$pred_occur == 0 & conting$occurrence == 1)
          absence_absence <- sum(conting$pred_occur == 0 & conting$occurrence == 0)
          
          sensitivity <- presence_presence / (presence_presence + absence_presence)
          specificity <- absence_absence / (absence_absence + presence_presence)
          false_pos_rate <- presence_absence / (presence_absence + absence_absence)
          accuracy <- (presence_presence + absence_absence) / (presence_presence + presence_absence + absence_presence + absence_absence)
          TSS <- sensitivity + specificity - 1
          ORSS <- (presence_presence * absence_absence - presence_absence * absence_presence) / (presence_presence * absence_absence + presence_absence * absence_presence)
          
          recap <- rbindlist(list(recap, list(taxa_selected,
               num.trees, mtry, min.node.size, max.depth, min.bucket,
               sensitivity, specificity, false_pos_rate, accuracy, TSS, ORSS,
               model_RA$prediction.error)),  use.names = F)
          }
        }
      }
    }
    return(recap)
}


#### Setting parameters for parallelization
cl<-makeCluster(30)
void<-  capture.output(clusterEvalQ(cl, c(library(dplyr),
                                          library(data.table), library("ranger"))))
clusterExport(cl, c("bio", "sites", "intersect", "stats_table_present", "stats_table_future",
                    "calculate_tss", "taxo", "list_datasets", "resolution"))

#### Launching function in parallel for all taxa and compile results in a table
recap_all_taxa <- parLapply(cl, taxa_list, function(taxa) {
  calculate_tss(taxa_selected = taxa)}
) %>% 
  rbindlist(use.names=F)

colnames(recap_all_taxa)=c("taxon",
                           "num.trees","mtry","min.node.size", "max.depth",
                           "min.bucket", "sensitivity", "specificity",
                           "false_pos_rate", "accuracy", "tss", "orss", "oob")

