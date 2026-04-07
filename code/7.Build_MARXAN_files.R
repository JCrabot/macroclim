### Working directory and packages ####
wd <- paste0(getwd(), "/data/")
setwd(wd)
pacman::p_load(dplyr, stringi, stringr, tidyr, # data wrangling
               terra, tidyterra #gis
               )

#### Set parameters of the model we are considering ####

#whether to discard "polluted" sites, can be "no", 10 for 10%, 50 for 50%
discard <- "no" # no or 10 or 50
# whether the model is considering fragmentation by dams or not
frag = "nofrag" # "" o "_frag"

# Set the names of the associated files and folders
file_name <- paste0(discard, "discard_nofrag")
file_name_frag <- paste0(discard, "discard_frag")
niche_folder <- prediction_folder <- paste0(discard, "discard/")

print(paste("___________________", file_name, "___________________" ))

#### Load data  ####

### Hydroatlas
if(!exists("hyriv_hybas", envir = globalenv())){
  
  # Getting the equivalence between ID of river reaches and ID of catchments
  hyriv_hybas <- read.csv(paste0(wd,"/environment/hydroatlas/hyriv_hybas_equivalence.csv"),
                          header=T, sep=",", stringsAsFactors = FALSE) %>%
    select(HYRIV_ID, HYBAS_L12)
 
  # Loading the catchments with additional column indicating the country
   hybas_country <-terra::vect(paste0(wd, "environment/hydroatlas/hybas_eu_lev01-12_v1c/hybas12_europe.shp")) %>%
    as.data.frame() %>%
    select(HYBAS_ID, NAME_EN) %>%
    rename(country = NAME_EN)
  
  # River network
  bbox_terra <- c(-10.5,32.702, 34.856, 71.31) # geographic extent of the study area
  rivers <-terra::vect(paste0(wd,"/environment/hydroatlas/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp"),
                       extent=bbox_terra) #already cropped to the right extent
  
  # Dataframe of the rivers with only required columns
  rivers_df <- rivers %>%
    as.data.frame() %>%
    select(HYRIV_ID, NEXT_DOWN, MAIN_RIV, LENGTH_KM)
}

### Load file of community temporal beta diversity for all rivers
temporal_turn <- read.csv(paste0(wd, "/biological/temporal_turnover.csv"),
                          header=T, sep=",", dec =".")

### Function for the fair allocation of catchment between countries
# to distribute evenly across countries the % to conserve across europe
allocate_fair <- function(nb_catch, target_tot) {
  
  # Where "nb_catch" is a vector of the total number of catchments where a genus
  # is predicted to occur in different countries and "target_tot" is the total
  # number of catchments that should ideally by monitored across Europe for the
  # given genus.
  # This function aims at distributing this number (target_tot) evenly
  # between the different countries where the genus occurs.
  
  n <- length(nb_catch) # number of countries
  remaining_target <- target_tot[1] # remaining catchments to allocate for monitoring
  remaining_cap <- nb_catch # remaining catchments available (where the genus "occurs")
  allocation <- rep(0, n) # initializing vector giving the number of catchments to monitor by country
  
  
  if (remaining_target[1] < n ) {
    
    # When there are very few catchments to allocate (less catchments
    # than there are countries of occurrence), it attributes randomly these few
    # catchments across countries.
    
    p <- round(remaining_target[1], 0)
    rand_country <- sample.int(n, p)
    allocation[rand_country] <- 1
    
  } else { 
    
    # In all other situations, as long as there are catchments to allocate to
    # countries for monitoring
    
    while (remaining_target >= 1) {
      
      # Identifying countries where there are available catchments
      active <- which(remaining_cap > 0)
      
      # Stop if there is no country with available catchments
      if (length(active) == 0) break
      
      # Splitting evenly the number of remaining catchments to monitor across
      # countries that have available catchments
      L <- remaining_target / length(active)
      
      for (i in active) { # for countries with available catchments
        
        # If it is possible, take L catchments in the country considered
        # If there are less, take all remaining catchments available of this country
        give <- min(L, remaining_cap[i])
        
        # Input this valor to the vector of catchments to monitor by country
        allocation[i] <- allocation[i] + give
        
        # Substract the number of catchments newly allocated to the 
        # number of catchments that remain to allocate
        remaining_target <- remaining_target - give
        
        # Substract the number of catchments newly allocated to the 
        # number of catchments that remain available in the country considered
        remaining_cap[i] <- remaining_cap[i] - give
        
      }
    }
    
  }
  
  return(allocation)
}

#### Building or loading table of all taxa occurrences by river for present and future ####

if (!file.exists(paste0(wd,"output/predictions/", prediction_folder,
                        "hydroriver_",  file_name, "_rich_present.csv"))){
  
  if(!exists("bio", envir = globalenv())){
    bio  <- read.csv(paste0(wd,"biological/bio_genus.csv"),
                     header=T, sep=",", stringsAsFactors = FALSE)
    rf_parameters <- read.csv(paste0(wd,"/biological/RF_parameters.csv"),
                              header=T, sep=",", stringsAsFactors = FALSE) %>%
      filter(!is.na(tss))
    
  }

  # Listing all files (one file per taxon) describing difference in predictions
  # between present and future
  list_file_predictions <- list.files(paste0(wd,"/output/predictions/",
                                             prediction_folder, ""),
                                      pattern = paste0("_diff"), full.names = TRUE)
  # Selecting only the files for the scenario considered
  list_file_predictions_selec <- list_file_predictions[str_remove(basename(list_file_predictions),
                                                                  paste0("_", discard,"discard_diff.csv")) %in% bio$taxon]
  
  # Prepare table for the present occurrences, with all river IDs
  hydroriver_taxa_present <- as.data.frame(matrix(nrow = nrow(rivers_df), ncol = 1)) %>%
    setNames("HYRIV_ID") %>%
    mutate(HYRIV_ID = rivers_df$HYRIV_ID)

  # Prepare table for the future occurrences, with all river IDs
  hydroriver_taxa_future <- as.data.frame(matrix(nrow = nrow(rivers_df), ncol = 1)) %>%
    setNames("HYRIV_ID") %>%
    mutate(HYRIV_ID = rivers_df$HYRIV_ID)

  for (i in 1:length(list_file_predictions_selec)){ # for all taxa
    
    # Reading the prediction file 
    prediction <- read.csv(list_file_predictions_selec[i], header=T, sep=",",
                           stringsAsFactors = FALSE) %>%
      # because in the "diff" file, pred_occur is the future occurrence (1 or 0)
      # and "diff" is the future-present occurrence
      rename(future = pred_occur) %>% 
      mutate(present = future - diff)

    # Merging table with the relevant column for present
    hydroriver_taxa_present <- hydroriver_taxa_present %>%
      left_join(prediction %>% select(HYRIV_ID, present) %>% # present occurrence
                  # getting the name of the genus
                  rename_with(~str_remove(basename(list_file_predictions_selec[i]),
                                          paste0("_", discard, "discard_diff.csv")), 2),
                by = "HYRIV_ID")

    # Merging table with the relevant column for future
    hydroriver_taxa_future <- hydroriver_taxa_future %>%
      left_join(prediction %>% select(HYRIV_ID, future) %>% # future occurence
                  # getting the name of the genus
                  rename_with(~str_remove(basename(list_file_predictions_selec[i]),
                                          paste0("_", discard,"discard_diff.csv")), 2),
                by = "HYRIV_ID")

  }

  # if fragmented scenario, take base file and change only future columns for aquatic taxa
  if (frag == "frag"){

    list_file_predictions_frag=list_file_predictions[grepl("frag", list_file_predictions)]
    for (i in 1:length(list_file_predictions_frag)){
      prediction <- read.csv(list_file_predictions_frag[i], header=T, sep=",",
                             stringsAsFactors = FALSE)
      taxo_frag <- str_remove(basename(list_file_predictions_frag[i]),
                              paste0("_", discard,"discard_diff_frag.csv"))
      hydroriver_taxa_future <- hydroriver_taxa_future %>%
        replace(grep(taxo_frag, colnames(hydroriver_taxa_future)),
                as.vector(prediction$future))
    }
  }

  # Filter taxa with TSS > 0.4 for present
  hydroriver_taxa_present <- hydroriver_taxa_present %>%
    mutate(HYRIV_ID = as.factor(HYRIV_ID)) %>%
    mutate_if(is.numeric, ~1 * (. != 0)) %>%
    select(any_of(c("HYRIV_ID", (rf_parameters %>% filter(tss >=0.4))$taxon))) %>%
    mutate(sum = rowSums(across(-HYRIV_ID)),
           HYRIV_ID = as.integer(as.character(HYRIV_ID)))

  # Write table for present
  write.csv2(hydroriver_taxa_present,
             paste0(wd,"output/predictions/", prediction_folder, "hydroriver_", 
                    file_name, "_rich_present.csv"), row.names = FALSE)

  # Filter taxa with TSS > 0.4 for future
  hydroriver_taxa_future <- hydroriver_taxa_future %>%
    mutate(HYRIV_ID = as.factor(HYRIV_ID)) %>%
    mutate_if(is.numeric, ~1 * (. != 0)) %>%
    select(any_of(c("HYRIV_ID", (rf_parameters %>% filter(tss >=0.4))$taxon))) %>%
    mutate(sum = rowSums(across(-HYRIV_ID)))

  # Write table for future
  write.csv2(hydroriver_taxa_future,
             paste0(wd,"output/predictions/", prediction_folder, "hydroriver_", 
                    file_name, "_rich_future.csv"), row.names = FALSE)


  # Get the overall table with difference in occurrences for all taxa
  hydroriver_taxa_diff <- cbind(hydroriver_taxa_future$HYRIV_ID,
                                (hydroriver_taxa_future[,-1] - hydroriver_taxa_present[,-1])) %>%
    rename("HYRIV_ID" = 1) %>%
    mutate(count_ext =  rowSums(across(-c(HYRIV_ID, sum)) == -1),
           count_col = rowSums(across(-c(HYRIV_ID, sum)) == 1))

  # Calculate percentage of change
  pc_change <- round((hydroriver_taxa_future$sum - hydroriver_taxa_present$sum)/ hydroriver_taxa_present$sum * 100, 2)
  hydroriver_taxa_diff <- hydroriver_taxa_diff %>%
    mutate(pc_diff = pc_change)

  # Write table for difference in occurrence
  write.csv2(hydroriver_taxa_diff,
             paste0(wd,"output/predictions/", prediction_folder, "hydroriver_",
                    file_name, "_rich_diff.csv"), row.names = FALSE)


} else { # if tables already exist

  hydroriver_taxa_present <- read.csv(paste0(wd,"output/predictions/",
                                             prediction_folder, "hydroriver_",  file_name, "_rich_present.csv"), header=T, sep=";", stringsAsFactors = FALSE)
  hydroriver_taxa_future <- read.csv(paste0(wd,"output/predictions/",
                                            prediction_folder, "hydroriver_",  file_name, "_rich_future.csv"), header=T, sep=";", stringsAsFactors = FALSE)
  hydroriver_taxa_diff <- read.csv(paste0(wd,"output/predictions/",
                                          prediction_folder, "hydroriver_",  file_name, "_rich_diff.csv"), header=T, sep=";", dec=",", stringsAsFactors = FALSE)

}

#### Start puvsp file - long community matrix ####

# MARXAN always requires at least three files: one of these three files is the
# community matrix in a long format (occurrence of each taxon at each site)
# This file is often called "puvsp".

# Starting to build the community matrix based on the present occurrences
puvsp_raw <- hydroriver_taxa_present %>%
  left_join(hyriv_hybas, by = "HYRIV_ID") %>% # Getting catchment ID
  select(- c(sum, HYRIV_ID)) %>% # Removing the species richness column
  # Grouping at the catchment scale
  select(HYBAS_L12, everything()) %>%
  group_by(HYBAS_L12) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE) %>%
  ungroup() %>%
  mutate(across(- HYBAS_L12,  ~1 * (. != 0))) %>%
  # Changing from wide to long format
  pivot_longer(!HYBAS_L12, names_to = "species", values_to = "amount") %>%
  rename(pu = HYBAS_L12) %>%
  mutate(species = as.factor(species)) %>%
  filter(amount>0)

# Adding information on the country, which will be used to compute a 
# fair allocation of the monitoring effort between countries.
puvsp <- puvsp_raw %>%
  right_join(hybas_country, by = c("pu" = "HYBAS_ID")) %>%
  # Creating a different ID for each genus in each country
  mutate(sp_country = paste0(species, country),
         sp_country = as.integer(as.factor(sp_country)))

# Table giving the equivalent between taxa name and taxa ID such as defined in
# the community matrix "puvsp"
equiv_spec_id <- puvsp %>%
  rename(species_id = sp_country) %>%
  select(-c(pu, amount)) %>%
  distinct()

# Small formating puvsp
puvsp <- puvsp %>%
  select(sp_country, pu, amount) %>%
  rename(species = sp_country)

#### Build pu file - site characteristics  ####

# MARXAN always requires at least three files: one of these three files is the 
# list of the planning units (pu), here the catchments, for which is indicated
# the cost it represents when being selected in a solution, and the "status", 
# which indicates whether to always ignore or always select (or neither) a 
# given planning unit.

# Getting the average temporal beta diversity across catchments
hybas_temp_turn <- temporal_turn %>%
  left_join(hyriv_hybas, by = "HYRIV_ID") %>%
  select(HYBAS_L12, everything()) %>%
  group_by(HYBAS_L12) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup()

# Setting the cost of the planning unit as the opposite of temporal beta
# diversity, because the higher is the temporal beta diversity, the lower we 
# want the cost of this catchment (i.e. it is interesting to monitor)
max_turn <- max(hybas_temp_turn$temp_turn)
cost_pu <- hybas_temp_turn %>%
  select(HYBAS_L12, temp_turn) %>%
  mutate(cost = max_turn - temp_turn)

# Formating the "pu" file exactly as required by MARXAN
# starting from the community matrix puvsp to have the same catchment IDs
pu <- puvsp %>%
  select(pu) %>%
  distinct() %>%
  mutate(pu = as.integer(as.character(pu))) %>%
  left_join(cost_pu %>%
              select(HYBAS_L12, temp_turn, cost), by = c("pu" = "HYBAS_L12")) %>%
  mutate(
    # the following line "locks in" (tells MARXAN to always select)
    # the units that have a temporal beta diversity greater than 0.50
    # by setting the status to "2"
    status = case_when(temp_turn >= 0.50 ~ 2, .default = 0) 
    ) %>%
  filter(!is.na(cost)) %>%
  select(pu, cost, status) %>%
  rename(id = pu)

#### Build spec file - species characteristics - and complete puvsp file ####

# MARXAN always requires at least three files: one of these three files is the
# "spec" file which specifies for each taxon what proportion (prop) of its 
# distribution area should be ideally selected in a solution, and how much
# a solution that does not select it should be penalized (spf = species
# penalty factor).

# Loading the macroclim score of each genus
taxa_weight <- read.csv(paste0(wd,"output/taxa_vuln_", file_name,".csv"),
                        header = TRUE, sep = ",", dec= ".")  %>%
  select(taxa, taxa_vuln) %>%
  rename(weight = taxa_vuln) %>%
  # calculating the exponential of the macroclim score
  mutate(exp_weight = exp(weight))

# Scaling the exponential of the macroclim score
min_score_exp = min(taxa_weight$exp_weight)
max_score_exp = max(taxa_weight$exp_weight)

taxa_weight <- taxa_weight %>%
  mutate(weight_calibrated =
           round((exp_weight - min_score_exp) / (max_score_exp - min_score_exp), 3))

# Merge the vulnerability scores with the ID of the genus of each country
spec <- equiv_spec_id %>%
  left_join(taxa_weight, by = c("species" = "taxa")) %>%
  mutate(spf = 100) %>% # input a SPF of 100 = maximum penalty for all taxa
  filter(!is.na(weight_calibrated)) %>%
  rename(name = species,
         id = species_id)

## Now starting to count how many catchments to monitor by country for each
# genus, with an effort of allocating evenly the number of catchments between
# countries. 

# First, counting the total number of catchments where the genus is predicted
# to occur
count_catch <- puvsp_raw %>%
  select(-c(amount, pu)) %>%
  group_by(species) %>%
  summarize(nb_catch_tot = length(species)) %>%
  left_join(taxa_weight, by = c("species" = "taxa")) %>%
  filter(!is.na(weight)) %>%
  rename(id = species)

# Now counting the number of catchments per country where the genus is predicted
# to occur
count_catch_by_country <- puvsp %>%
  select(-c(amount, pu)) %>%
  group_by(species) %>%
  summarize(nb_catch = length(species)) %>%
  # joining the table with the vulnerability score
  left_join(spec %>% select(name, id, spf), by = c("species" = "id")) %>%
  # joining the table with the name and ID of all genera
  left_join(equiv_spec_id %>% select(species_id, country),
            by = c("species" = "species_id")) %>%
  # joining the table with the total number of catchments
  left_join(count_catch, by = c("name" = "id")) %>%
  filter(!is.na(weight)) %>%
  # Defining a first target number of catchments to monitor as 
  # (scaled exponential vulnerability) * (1 percent of all European catchments)
  mutate(target_tot = weight_calibrated*459) %>%
  # Now distributing this number of catchments fairly between country using
  # the "allocate_fair" function defined in the "Load data" chunk
  group_by(name) %>%
  mutate(allocation = allocate_fair(nb_catch, first(target_tot)),
         # transforming this number of catchments into a proportion
         prop_alloc = allocation/nb_catch) %>% 
  ungroup()

# Final formatting of the "spec" table for MARXAN
spec <- count_catch_by_country %>%
  select(species, prop_alloc, spf, name) %>%
  rename(prop = prop_alloc,
         id = species)

### Making sure the puvsp file only contains planning units  that are in the pu
# table and taxa that are in the spec table
pu_id <- pu$id
puvsp <- puvsp %>%
  filter(species %in% spec$id) %>%
  filter(pu %in% pu_id)

#### Export files ####
# Create a subfolder for the MARXAN analysis of the considered model
wd_marxan <- paste0(wd,"/output/marxan/", file_name, "/input/")
dir.create(wd_marxan, recursive = T)

# Write tables
# Planning units
write.csv(pu, paste0(wd_marxan,"pu.csv"), row.names = F, quote = F) 
# Genus list
write.csv(spec, paste0(wd_marxan,"spec.csv"), row.names = F, quote = F)
# Community matrix
write.csv(puvsp, paste0(wd_marxan,"puvsp.csv"), row.names = F, quote = F)
