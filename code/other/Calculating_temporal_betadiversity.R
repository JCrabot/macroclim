#### The files used in this script (e.g. hydroriver_taxa_present.csv) are the
# same as in the "7.Build_MARXAN_files.R" script and can be loaded using it.

#### Calculating temporal beta diversity ####

# Creating empty table
beta <- data.frame(matrix(ncol=2))

for (i in 1:nrow(hydroriver_taxa_present)){
  
  # Calculation of temporal beta diversity
  betapair <- vegan::vegdist(rbind(hydroriver_taxa_present[i,-1],
                                   hydroriver_taxa_future[i,-1]), method = "bray")
  beta <- rbind(beta,
                c(hydroriver_taxa_present$HYRIV_ID[i], betapair))
  
}
beta <- beta[-1,]

# Adding column names
names(beta) <- c("HYRIV_ID", "temp_turn")

# # Exporting resulting table
# write.csv(beta, paste0(wd,"biological/temporal_turnover.csv"), row.names = FALSE)

#### Exporting shapefile at river scale ####
if(!exists("rivers", envir = globalenv())){ # if river shapefile not loaded yet

  bbox_terra <- c(-10.5,32.702, 34.856, 71.31) # geographic extent of the study area
  rivers <-terra::vect(paste0(wd,"/environment/hydroatlas/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp"),
                       extent=bbox_terra) #already cropped to the right extent
}

# Merge result with river network of HydroRIVERS
hydroriver_beta_shp <- merge(rivers, beta, by = "HYRIV_ID")

# Export shapefile
writeVector(hydroriver_beta_shp,
            filename = paste0(wd, "gis/hydroriver_temporal_turnover.shp") ,
            overwrite=TRUE)

#### Exporting shapefile at catchment scale ####

# Getting the equivalence between ID of river reaches and ID of catchments
hyriv_hybas <- read.csv(paste0(wd,"environment/hydroatlas/hyriv_hybas_equivalence.csv"),
                        header=T, sep=",", stringsAsFactors = FALSE) %>%
  select(HYRIV_ID, HYBAS_L12)

# Loading the temporal beta diversity
temporal_turn <- read.csv(paste0(wd, "biological/temporal_turnover.csv"),
                          header=T, sep=";", dec =".")

# Averaging the temporal beta diversity at the scale of the catchment
hybas_temp_turn <- temporal_turn %>%
  left_join(hyriv_hybas, by = "HYRIV_ID") %>%
  select(HYBAS_L12, everything()) %>%
  group_by(HYBAS_L12) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  rename(HYBAS_ID = HYBAS_L12)

# Recalling spatial extent if needed
bbox_terra <- c(-10.5,32.702, 34.856, 71.31) 

# Loading catchments of HydroBasins
basins <-terra::vect(paste0(wd,"hydroatlas/hybas_eu_lev01-12_v1c/hybas_eu_lev12_v1c.shp"),
                     extent=bbox_terra) #already cropped to the right extent

# Merging results on beta diversity with the catchments of HydroBasins
hydrobasin_beta_shp <- merge(basins, hybas_temp_turn, by = "HYBAS_ID")

# Exporting shapefile
writeVector(hydrobasin_beta_shp,
            filename = paste0(wd, "output/hydrobasin_genus_temporal_turnover.shp") ,
            overwrite=TRUE)
