##        Environmental data wrangling before modeling         ##
#################################################################
#### EXTRACT HYDROSHED ID OF RIVER SEGMENTS WHERE THE TAXA HAS BEEN OBSERVED ####

spdata_ids <- as.data.frame(sites[sites$site_id %in% species$site_id,])

#### Extract HYDROSHED ID of polluted river segments if necessary ####
if (discard != "no"){ 
  if (!file.exists(paste0(wd,"environment/list_polluted_segments_", as.integer(discard),".csv"))){
    rivers <-terra::vect(paste0(wd,"hydroatlas/RiverATLAS_v10_shp/RiverATLAS_v10_eu.shp"),
                         extent=bbox_terra) #already cropped to the right extent
    pc = (100 - as.integer(discard))/100
    list_less_polluted_reaches <- (rivers %>%
                                     as.data.frame() %>%
                                     filter(HYRIV_ID %in% (intersect$HYRIV_ID)) %>%
                                     select(HYRIV_ID, pop_ct_csu, pop_ct_usu, urb_pc_cse, urb_pc_use,
                                            hft_ix_c09, hft_ix_u09) %>%
                                     filter(hft_ix_c09 < quantile(hft_ix_c09, pc))
    )$HYRIV_ID
  
    write.csv(list_less_polluted_reaches, paste0(wd,"environment/list_polluted_segments_", discard,".csv"), row.names = FALSE)
  
  } else {
    list_less_polluted_reaches <- read.csv(paste0(wd,"environment/list_polluted_segments_", discard,".csv"), header=T, sep=",")
    names(list_less_polluted_reaches) <- "HYRIV_ID"
  }
}

#### Preparing data for modeling ####

if (!exists("stats_table_present")) { # same table for all species
  
  # Read the file of environmental variables
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
                "riv_tc_usu", "gwt_cm_cav", "slt_pc_cav", "snd_pc_uav", "ele_mt_cmn", "aet_mm_cyr",
                "ero_kh_cav",
                "yr_irclass"
  )

  stats_table_present <- stats_table_present[,..var_final]
  stats_table_future <- stats_table_future[,..var_final]

  # Clear up memory
  rm(list=c("stats_table_zon","stats_table")) ; gc()
  
}
