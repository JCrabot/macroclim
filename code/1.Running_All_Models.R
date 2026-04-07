#### Working directory and packages ####
pacman::p_load(here, #working directory
               data.table, dplyr, janitor, tidyverse, stringr, tibble, #data wrangling
               terra, tidyterra, # gis
               foreach, doParallel, # parallelize fragmentation calculation
               ranger #random forests
)

wd <- paste0(here(), "/data/")
wd_script <- paste0(here(), "/code/")

setwd(wd)

#### Specify model parameters used in the random forest ####
proba <- FALSE # TRUE to get a probability of occurrence in the random forests (else = binary result presence/absence)
export <- FALSE # if TRUE writes predictions to disk with to display maps later in GIS
downsample <- TRUE # if TRUE, down-sampling carried out for the RF, if FALSE, all sites considered for the RF

#### Loading data ####

### Biological information
# Community matrix
bio  <- read.csv(paste0(wd,"biological/bio_genus.csv"),
                 header=T, sep=",", stringsAsFactors = FALSE)
# Information on genus occurrences across sites
summary_occurences <- read.csv(paste0(wd,"biological/summary_occurences.csv"),
                               header=T, sep=",", stringsAsFactors = FALSE)
# List of datasets to be considered for the SDM for each genus
list_datasets <- read.csv(paste0(wd,"biological/datasets_by_genus.csv"),
                          header=T, sep=";",  stringsAsFactors = FALSE)
# Parameters to use for each genus in the random forests, obtained after parameter tuning
rf_parameters <- read.csv(paste0(wd,"biological/RF_parameters.csv"),
                          header=T, sep=",", stringsAsFactors = FALSE) %>%
  filter(!is.na(tss))
# Classification table
taxo  <- read.csv(paste0(wd,"/biological/classification_macroclim.csv"),
                  header=T, sep=",", stringsAsFactors = TRUE)

### Site information
sites  <- read.csv(paste0(wd,"biological/sites.csv"),
                   header=T, sep=";", stringsAsFactors = FALSE)

### HydroATLAS data
# Equivalence of ID of river reaches from HydroRivers and ID of catchments of HydroBasins
hyriv_hybas <- read.csv(paste0(wd,"environment/hydroatlas/hyriv_hybas_equivalence.csv"),header=T,sep=",")
bbox_terra <- c(-10.5,32.702, 34.856, 71.31) # geographic extent of the study area
# Catchment polygons of HydroBasins at the &é_th hierarchical level
basins12 <- vect(paste0(wd, "environment/hydroatlas/hybas_eu_lev01-12_v1c/hybas_eu_lev12_v1c.shp"),
                 extent=bbox_terra)

rivers <-terra::vect(paste0(wd,"/environment/hydroatlas/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp"),
                     extent=bbox_terra) #already cropped to the right extent
rivers_df <- rivers %>%
  as.data.frame() %>%
  select(HYRIV_ID, NEXT_DOWN, MAIN_RIV, LENGTH_KM)

### Short function to monitor the running time during the loop below
log_time <- function(msg = "") {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n")
}

# Rescaling function
rescale <- function(x){(x-min(x))/(max(x)-min(x))}

#### Running all models ####

for (discard in c("no", "10", "50")){  # #whether to discard "polluted" sites, can be "no", 10 for 10%, 50 for 50%

  ## Specifying name of output files and folders ##
  file_name <- paste0(discard, "discard_nofrag")
  file_name_frag <- paste0(discard, "discard_frag")
  niche_folder <- prediction_folder <- paste0(discard, "discard/")

  cat("----------------------- Discard parameter:", discard, "---------------------")

  log_time("--------- Running the RF model for all taxa
           and extracting descriptors of distribution present and future ---------")
  source(paste0(wd_script,"2.0_SDM_and_distribution_descriptors.R"))
  
  log_time("--------- Calculating descriptors of distribution change ---------")
  source(paste0(wd_script,"3.Calculate_distribution_change.R"))
    
  log_time("--------- Calculating new predictions considering fragmentation ---------")
  source(paste0(wd_script,"4.New_predictions_with_fragmentation.R"))
  
  log_time("--------- Calculating new descriptors of distribution change with fragmentation ---------")
  source(paste0(wd_script,"5.Calculate_distribution_change_after_fragmentation.R"))
  
  log_time("--------- Calculating taxa vulnerability scores ---------")
  source(paste0(wd_script,"6.Calculate_vulnerability_scores.R"))
  
  log_time("--------- End of model analysis ---------")
      
}
