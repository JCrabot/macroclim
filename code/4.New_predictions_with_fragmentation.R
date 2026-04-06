
# intersection of AMBER dams and Hydrosehd, obtained in QGis
if (!exists("dams_intersected", envir = globalenv())) {
  dams_intersected  <- fread(paste0(wd,
                                    "environment/amber-project-data/intersect_amber_dams_hydroshed_0.005.csv")) %>%
    filter(!is.na(HYRIV_ID))
  dams_id <- dams_intersected$HYRIV_ID
}

if (!exists("river_segments", envir = globalenv())) {
  river_segments <-terra::vect(paste0(wd,"environment/hydroatlas/RiverATLAS_v10_shp/RiverATLAS_v10_eu.shp"),
                             extent=bbox_terra)#already cropped to the right extent
}

if (!exists("rivers_selec_df", envir = globalenv())) {
  rivers_selec <-terra::vect(paste0(wd,"environment/hydroatlas/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp"),                          
                           extent=bbox_terra) #already cropped to the right extent
  rivers_selec_df <- rivers_selec %>%
    as.data.frame() %>%
    select(HYRIV_ID, NEXT_DOWN, MAIN_RIV, LENGTH_KM, ORD_STRA) %>%
    as.data.table()
  
  rivers_selec_df[, nsegs_network := .N, by=MAIN_RIV]
  rivers_inland <- merge(rivers_selec_df[nsegs_network > 1,], 
                         rivers_selec_df[nsegs_network > 1, .(HYRIV_ID, LENGTH_KM, ORD_STRA)],
                         by.x='NEXT_DOWN', by.y='HYRIV_ID', 
                         all.x=T,
                         suffixes=c("", "_down")) %>%
    .[, .(HYRIV_ID, NEXT_DOWN, MAIN_RIV, ORD_STRA, ORD_STRA_down,
          LENGTH_KM, LENGTH_KM_down, nsegs_network)] 
  
  rivers_coastal <- rivers_selec_df[nsegs_network == 1,]
}

#### Selecting non-flying taxa ####
list_taxa <- unique(bio$taxon[bio$taxon %in% (taxo %>% filter(class !="Insecta"))$taxon])
list_taxa <- list_taxa[order(list_taxa)]

#### Set parameters or parallelization ####
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#### Fragmentation loop on all non-flying taxa - using distance length ####

# Bound is a file that was used to look at the "boundary" as in Da Silva et al. 2024
# From this "bound" variable, we get the river network distance in km between the outlets
# of two river reaches along the HydroRIVERS network 
length <- fread(paste0(wd,"/environment/bound.csv"), sep  = ",", dec=".", header= T) %>%
  left_join(rivers_selec_df, by = c("id1" = "HYRIV_ID")) %>% 
  mutate(dist_km = 1/(bound * bound)) %>%
  rename(HYRIV_ID_from = id1,
         HYRIV_ID_to = id2) %>%
  as.data.table()

log_time("Starting fragmentation loop")

for (taxa in list_taxa){ 
  if(!file.exists(paste0(wd,"output/predictions/",prediction_folder, taxa, "_" , discard, "discard_diff_frag.csv"))){
   if(file.exists(paste0(wd,"output/predictions/", prediction_folder, taxa,  "_", discard, "discard_diff.csv"))){
    
    pred <- fread(paste0(wd,"output/predictions/", prediction_folder, taxa, "_", discard,  "discard_diff.csv")) %>%
      setnames(c('HYRIV_ID', 'future', 'diff')) %>%
      .[, present := future - diff] 
    
    reaches_present_today <- pred[present ==1, HYRIV_ID]
    
    ### Identify newly colonized reaches
    id_new_reaches <-  pred[diff ==1, HYRIV_ID]
    
    unique_main <- rivers_inland[HYRIV_ID %in% id_new_reaches, unique(MAIN_RIV)]
    
    if (length(id_new_reaches) > 0 & length(unique_main) >0){
      
      id_no_colonization <- foreach (i = 1:length(unique_main),
                                     .packages=c("dplyr", "data.table", "igraph"),
                                     .combine = c) %dopar% { # for each new colonized reach
                                       
                                       in_main_riv = unique_main[i]
                                       length_sub = length[MAIN_RIV == in_main_riv,]
                                       
                                       # selectioning in the distance matrices all pairs of sites where (HYRIV_ID_to) there was a colonization downstream (HYRIV_ID_from) a presence today
                                       dist_mat_to <- length_sub[as.character(HYRIV_ID_to) %in% id_new_reaches,]
                                       dist_mat_to <- dist_mat_to[as.character(HYRIV_ID_from) %in% reaches_present_today,]
                                       colonization_effective <- dist_mat_to$HYRIV_ID_to # keeping colonizations if there was a presence upstream
                                       
                                       # removing colonization if there was no presence upstream nor downstream
                                       dist_mat_from <- length_sub[as.character(HYRIV_ID_from) %in% id_new_reaches,]
                                       dist_mat_from <- dist_mat_from[!as.character(HYRIV_ID_from) %in% colonization_effective,]
                                       no_col_downstream<- dist_mat_from[!as.character(HYRIV_ID_to) %in% reaches_present_today,HYRIV_ID_from]
                                       dist_mat_from <- dist_mat_from[as.character(HYRIV_ID_to) %in% reaches_present_today,] %>%
                                         mutate(dam = case_when(HYRIV_ID_to %in% dams_id ~ "dam",
                                                                .default = "no_dam")) %>%
                                         arrange(dist_km, dam) %>%
                                         group_by(HYRIV_ID_from)%>%
                                         slice(1) %>%
                                         ungroup
                                       no_col_dam <- dist_mat_from[which(dist_mat_from$dam == "dam"),]$HYRIV_ID_from
                                       
                                       
                                       if (length(id_new_reaches[id_new_reaches %in% rivers_coastal$HYRIV_ID]) > 0) {
                                         id_no_colonization <- c(no_col_downstream, no_col_dam, id_new_reaches[id_new_reaches %in% rivers_coastal$HYRIV_ID])
                                       } else {
                                         id_no_colonization <- c(no_col_downstream, no_col_dam)
                                       }
                                       
                                       return(id_no_colonization)
                                       
                                     }
    
      pred <- pred %>%
        mutate(future = case_when(HYRIV_ID %in% id_no_colonization ~ 0,
                                  .default = future),
               diff = case_when(HYRIV_ID %in% id_no_colonization ~ 0,
                                .default = diff),
               change = case_when(HYRIV_ID %in% id_no_colonization ~ "CHANGED",
                                  .default = NA))
      write.csv(pred, paste0(wd,"output/predictions/", prediction_folder, taxa, "_" , discard, "discard_diff_frag.csv"), row.names = FALSE)
    
    }
   }
  }
}

stopCluster(my.cluster)  # closing the nodes used in parallelization

log_time("Fragmentation files saved") # giving time