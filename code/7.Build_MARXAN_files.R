### Working directory and packages ####
wd <- paste0(getwd(), "/data/")
setwd(wd)
pacman::p_load(dplyr, stringi, stringr, tidyr, # data wrangling
               terra, tidyterra #gis
               )
options(dplyr.summarise.inform = FALSE)
#### Set parameters ####

# discard <- "no" # no or 10 or 50
# frag = "nofrag" # "" o "_frag"

file_name <- paste0(discard, "discard_nofrag")
file_name_frag <- paste0(discard, "discard_frag")
niche_folder <- prediction_folder <- paste0(discard, "discard/")

print(paste("___________________", file_name, "___________________" ))

#### Load data  ####

### Hydroatlas
if(!exists("hyriv_hybas", envir = globalenv())){
  hyriv_hybas <- read.csv(paste0(wd,"/environment/hydroatlas/hyriv_hybas_equivalence.csv"),
                          header=T, sep=",", stringsAsFactors = FALSE) %>%
    select(HYRIV_ID, HYBAS_L12)
  hybas_country <-terra::vect(paste0(wd, "environment/hydroatlas/hybas_eu_lev01-12_v1c/hybas12_europe.shp")) %>%
    as.data.frame() %>%
    select(HYBAS_ID, NAME_EN) %>%
    rename(country = NAME_EN)
  
  bbox_terra <- c(-10.5,32.702, 34.856, 71.31) # geographic extent of the study area
  rivers <-terra::vect(paste0(wd,"/environment/hydroatlas/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp"),
                       extent=bbox_terra) #already cropped to the right extent
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
  
  n <- length(nb_catch)
  remaining_target <- target_tot[1]
  remaining_cap <- nb_catch
  allocation <- rep(0, n)
  
  if (remaining_target[1] < n ) {
    p <- round(remaining_target[1], 0)
    rand_country <- sample.int(n, p)
    allocation[rand_country] <- 1
    
  } else {
    
    while (remaining_target >= 1) {
      active <- which(remaining_cap > 0)
      if (length(active) == 0) break
      
      L <- remaining_target / length(active)
      
      for (i in active) {
        give <- min(L, remaining_cap[i])
        allocation[i] <- allocation[i] + give
        remaining_target <- remaining_target - give
        remaining_cap[i] <- remaining_cap[i] - give
      }
    }
    
  }
  
  return(allocation)
}

#### Building or loading table of taxa occurrence by river for present and future ####

if (!file.exists(paste0(wd,"output/predictions/", prediction_folder,
                        "hydroriver_",  file_name, "_rich_present.csv"))){
  
  if(!exists("bio", envir = globalenv())){
    bio  <- read.csv(paste0(wd,"biological/bio_genus.csv"),
                     header=T, sep=",", stringsAsFactors = FALSE)
    rf_parameters <- read.csv(paste0(wd,"/biological/RF_parameters.csv"),
                              header=T, sep=",", stringsAsFactors = FALSE) %>%
      filter(!is.na(tss))
    
  }

  list_file_predictions <- list.files(paste0(wd,"/output/predictions/", prediction_folder, ""),
                                      pattern = paste0("_diff"), full.names = TRUE)
  list_file_predictions_selec <- list_file_predictions[str_remove(basename(list_file_predictions),
                                                                  paste0("_", discard,"discard_diff.csv")) %in% bio$taxon]

  hydroriver_taxa_present <- as.data.frame(matrix(nrow = nrow(rivers_df), ncol = 1)) %>%
    setNames("HYRIV_ID") %>%
    mutate(HYRIV_ID = rivers_df$HYRIV_ID)

  hydroriver_taxa_future <- as.data.frame(matrix(nrow = nrow(rivers_df), ncol = 1)) %>%
    setNames("HYRIV_ID") %>%
    mutate(HYRIV_ID = rivers_df$HYRIV_ID)

  for (i in 1:length(list_file_predictions_selec)){
    print(i)
    prediction <- read.csv(list_file_predictions_selec[i], header=T, sep=",", stringsAsFactors = FALSE) %>%
      rename(future = pred_occur) %>%
      mutate(present = future - diff)

    hydroriver_taxa_present <- hydroriver_taxa_present %>%
      left_join(prediction %>% select(HYRIV_ID, present) %>%
                  rename_with(~str_remove(basename(list_file_predictions_selec[i]), paste0("_", discard,"discard_diff.csv")), 2),
                by = "HYRIV_ID")

    hydroriver_taxa_future <- hydroriver_taxa_future %>%
      left_join(prediction %>% select(HYRIV_ID, future) %>%
                  rename_with(~str_remove(basename(list_file_predictions_selec[i]), paste0("_", discard,"discard_diff.csv")), 2),
                by = "HYRIV_ID")

  }

  # if fragmented scenario, take base file and change only future columns for aquatic taxa
  if (frag == "frag"){

    list_file_predictions_frag=list_file_predictions[grepl("frag", list_file_predictions)]
    for (i in 1:length(list_file_predictions_frag)){
      prediction <- read.csv(list_file_predictions_frag[i], header=T, sep=",", stringsAsFactors = FALSE)
      taxo_frag <- str_remove(basename(list_file_predictions_frag[i]),
                              paste0("_", discard,"discard_diff_frag.csv"))
      hydroriver_taxa_future <- hydroriver_taxa_future %>%
        replace(grep(taxo_frag, colnames(hydroriver_taxa_future)), as.vector(prediction$future))
    }
  }

  hydroriver_taxa_present <- hydroriver_taxa_present %>%
    mutate(HYRIV_ID = as.factor(HYRIV_ID)) %>%
    mutate_if(is.numeric, ~1 * (. != 0)) %>%
    select(any_of(c("HYRIV_ID", (rf_parameters %>% filter(tss >=0.4))$taxon))) %>%
    mutate(sum = rowSums(across(-HYRIV_ID)),
           HYRIV_ID = as.integer(as.character(HYRIV_ID)))

  write.csv2(hydroriver_taxa_present,
             paste0(wd,"output/predictions/", prediction_folder, "hydroriver_",  file_name, "_rich_present.csv"), row.names = FALSE)

  hydroriver_taxa_future <- hydroriver_taxa_future %>%
    mutate(HYRIV_ID = as.factor(HYRIV_ID)) %>%
    mutate_if(is.numeric, ~1 * (. != 0)) %>%
    select(any_of(c("HYRIV_ID", (rf_parameters %>% filter(tss >=0.4))$taxon))) %>%
    mutate(sum = rowSums(across(-HYRIV_ID)))

  write.csv2(hydroriver_taxa_future,
             paste0(wd,"output/predictions/", prediction_folder, "hydroriver_",  file_name, "_rich_future.csv"), row.names = FALSE)


  hydroriver_taxa_diff <- cbind(hydroriver_taxa_future$HYRIV_ID,
                                (hydroriver_taxa_future[,-1] - hydroriver_taxa_present[,-1])) %>%
    rename("HYRIV_ID" = 1) %>%
    mutate(count_ext =  rowSums(across(-c(HYRIV_ID, sum)) == -1),
           count_col = rowSums(across(-c(HYRIV_ID, sum)) == 1))

  pc_change <- round((hydroriver_taxa_future$sum - hydroriver_taxa_present$sum)/ hydroriver_taxa_present$sum * 100, 2)
  hydroriver_taxa_diff <- hydroriver_taxa_diff %>%
    mutate(pc_diff = pc_change)

  write.csv2(hydroriver_taxa_diff,
             paste0(wd,"output/predictions/", prediction_folder, "hydroriver_",  file_name, "_rich_diff.csv"), row.names = FALSE)


} else {

  hydroriver_taxa_present <- read.csv(paste0(wd,"output/predictions/", prediction_folder, "hydroriver_",  file_name, "_rich_present.csv"), header=T, sep=";", stringsAsFactors = FALSE)
  hydroriver_taxa_future <- read.csv(paste0(wd,"output/predictions/", prediction_folder, "hydroriver_",  file_name, "_rich_future.csv"), header=T, sep=";", stringsAsFactors = FALSE)
  hydroriver_taxa_diff <- read.csv(paste0(wd,"output/predictions/", prediction_folder, "hydroriver_",  file_name, "_rich_diff.csv"), header=T, sep=";", dec=",", stringsAsFactors = FALSE)

}
#### Start puvsp file - long community matrix ####
puvsp_raw <- hydroriver_taxa_present %>%
  left_join(hyriv_hybas, by = "HYRIV_ID") %>%
  select(- c(sum, HYRIV_ID)) %>%
  select(HYBAS_L12, everything()) %>%
  group_by(HYBAS_L12) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE) %>%
  ungroup() %>%
  mutate(across(- HYBAS_L12,  ~1 * (. != 0))) %>%
  pivot_longer(!HYBAS_L12, names_to = "species", values_to = "amount") %>%
  rename(pu = HYBAS_L12) %>%
  mutate(species = as.factor(species)) %>%
  filter(amount>0)

puvsp <- puvsp_raw %>%
  right_join(hybas_country, by = c("pu" = "HYBAS_ID")) %>%
  mutate(sp_country = paste0(species, country),
         sp_country = as.integer(as.factor(sp_country)))

equiv_spec_id <- puvsp %>%
  rename(species_id = sp_country) %>%
  select(-c(pu, amount)) %>%
  distinct()

#### Build pu file - site characteristics  ####
hybas_temp_turn <- temporal_turn %>%
  left_join(hyriv_hybas, by = "HYRIV_ID") %>%
  select(HYBAS_L12, everything()) %>%
  group_by(HYBAS_L12) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup()

max_turn <- max(hybas_temp_turn$temp_turn)
macroclim_diff <- hybas_temp_turn %>%
  select(HYBAS_L12, temp_turn) %>%
  mutate(cost_diff = max_turn - temp_turn)

pu <- puvsp %>%
  select(pu) %>%
  distinct() %>%
  mutate(pu = as.integer(as.character(pu))) %>%
  left_join(macroclim_diff %>%
              select(HYBAS_L12, temp_turn, cost_diff), by = c("pu" = "HYBAS_L12")) %>%
  mutate(
    cost = cost_diff,
    status = case_when(temp_turn >= 0.50 ~ 2, .default = 0)
    ) %>%
  filter(!is.na(cost)) %>%
  select(pu, cost, status) %>%
  rename(id = pu)

#### Build spec file - species characteristics - and complete puvsp file ####
taxa_weight <- read.csv(paste0(wd,"output/taxa_vuln_", file_name,".csv"),
                        header = TRUE, sep = ",", dec= ".")  %>%
  select(taxa, taxa_vuln) %>%
  rename(weight = taxa_vuln) %>%
  mutate(exp_weight = exp(weight))

min_score_exp = min(taxa_weight$exp_weight)
max_score_exp = max(taxa_weight$exp_weight)

taxa_weight <- taxa_weight %>%
  mutate(weight_calibrated =
           round((exp_weight - min_score_exp) / (max_score_exp - min_score_exp), 3))

spec <- equiv_spec_id %>%
  left_join(taxa_weight, by = c("species" = "taxa")) %>%
  mutate(spf = 100) %>%
  filter(!is.na(weight_calibrated)) %>%
  rename(name = species,
         id = species_id)

### count catchments per species
count_catch <- puvsp_raw %>%
  select(-c(amount, pu)) %>%
  group_by(species) %>%
  summarize(nb_catch_tot = length(species)) %>%
  left_join(taxa_weight, by = c("species" = "taxa")) %>%
  filter(!is.na(weight)) %>%
  rename(id = species)

puvsp <- puvsp %>%
  select(sp_country, pu, amount) %>%
  rename(species = sp_country)

count_catch_by_country <- puvsp %>%
  select(-c(amount, pu)) %>%
  group_by(species) %>%
  summarize(nb_catch = length(species)) %>%
  left_join(spec %>% select(name, id, spf), by = c("species" = "id")) %>%
  left_join(equiv_spec_id %>% select(species_id, country), by = c("species" = "species_id")) %>%
  left_join(count_catch, by = c("name" = "id")) %>%
  filter(!is.na(weight)) %>%
  mutate(target_tot = weight_calibrated*459) %>%
  group_by(name) %>%
  mutate(allocation = allocate_fair(nb_catch, first(target_tot)),
         prop_alloc = allocation/nb_catch) %>%
  ungroup()

spec <- count_catch_by_country %>%
  select(species, prop_alloc, spf, name) %>%
  rename(prop = prop_alloc,
         id = species)

pu_id <- pu$id

puvsp <- puvsp %>%
  filter(species %in% spec$id) %>%
  filter(pu %in% pu_id)


#### Export files ####
wd_marxan <- paste0(wd,"/output/marxan/", file_name, "/input/")
dir.create(wd_marxan, recursive = T)

write.csv(pu, paste0(wd_marxan,"pu.csv"), row.names = F, quote = F)
write.csv(spec, paste0(wd_marxan,"spec.csv"), row.names = F, quote = F)
write.csv(puvsp, paste0(wd_marxan,"puvsp.csv"), row.names = F, quote = F)
