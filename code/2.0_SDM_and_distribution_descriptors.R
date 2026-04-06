#### STARTING GIANT LOOP ####

list_taxa <- unique((bio %>% filter(resolution == resolution))$taxon)
list_taxa <- list_taxa[order(list_taxa)]

for (taxa_selected in list_taxa) { 

  print(taxa_selected)

  source(paste0(wd_script,"2.1_Biological_data_wrangling.R"))
  source(paste0(wd_script,"2.2_Environmental_data_wrangling.R"))
  if(nrow(spdata_ids)>20){
    source(paste0(wd_script,"2.3_SpeciesDistributionModels.R"))
    source(paste0(wd_script,"2.4_Extract_model_information.R"))
  }

  rm(list=setdiff(ls(), c("wd","wd_script","env", "bio", "sites", "rivers", "river_segments", "bbox_terra", # raw data to keep
                          "geom_rivers","list_less_polluted_reaches", "list_taxa", "taxo", 
                          "log_time", "rescale", #functions
                          "list_datasets", "summary_occurences", # to check conditions to run script at genus level
                          "stats_table_present","stats_table_future", #environmental variables for SDMs
                          "rivers_to_keep", "hyriv_hybas", "basins12", "rivers_df", "rivers", "rivers_selec", "rivers_selec_df", # rivers data
                          "length", "dams_intersected", "dams_id",# used for fragmentation
                          "type", "proba", "export", "random", "discard", "fragmentation", "downsample", "selec", "rf_parameters", # parameters for functions
                          "prediction_folder", "niche_folder", "file_name", "file_name_frag" #output names
                          ))) 
  gc()
  
}

