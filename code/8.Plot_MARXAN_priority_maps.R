#### Working environment ####
wd <- paste0(getwd(), "/data/")
setwd(wd)
pacman::p_load(dplyr, tidyterra, terra) # Data wrangling and map plots

# Create folder for the output shapefiles of MARXAN if needed
if (!dir.exists(paste0(wd, "output/marxan/gis/"))){
  dir.create(paste0(wd, "output/marxan/gis/"), recursive=T)
}

#### Set parameters ####

#whether to discard "polluted" sites, can be "no", 10 for 10%, 50 for 50%
discard <- "no"
# whether to consider the fragmented scenario ("_frag") or the base scenario
# without fragmentation ("nofrag")
frag <- "nofrag" # "" o "_frag"

### Setting associated names for files and folders
file_name <- paste0(discard, "discard_nofrag")
file_name_frag <- paste0(discard, "discard_frag")
niche_folder <- prediction_folder <- paste0(discard, "discard/")

print(paste("___________________", file_name, "___________________" ))


#### Run MARXAN ####
# Defining the output folder for the model considered
wd_marxan <- paste0(wd, "output/marxan/", file_name)
setwd(wd_marxan)

# Copy/ pasting the files necessary to run MARXAN from the root folder to the
# folder of the model considered
file.copy(paste0(wd, "/marxan/Marxan_x64.exe"), # execution command
          wd_marxan)
file.copy(paste0(wd, "/marxan/input.dat"), # indications called in the execution
          wd_marxan)

# Creating a folder to save the MARXAN solutions for this model
dir.create(paste0(wd_marxan, "/output"), recursive = T)

# Run MARXAN
system("Marxan_x64.exe", wait = T, invisible = T) 

#### Load HydroATLAS shapefile ####
if (! ("basins" %in% ls(envir = .GlobalEnv))) { # if catchments not already loaded
  bbox_terra <- c(-10.5,32.702, 34.856, 71.31) # spatial extent
  basins <- terra::vect(paste0(wd,
                               "/environment/hydroatlas/hybas_eu_lev01-12_v1c/hybas_eu_lev12_v1c.shp"),
                        extent = bbox_terra)
}

#### Export solution shapefile ####

# Among all the output files of MARXAN (one file for each solution), the
# "output_best.csv" is the best solution, as it minimized the cost. 
# 1 indicates that the associated planning unit is selected, 0 that it is not.
file_best <- read.csv(paste0(wd_marxan, "/output/output_best.csv"), 
                      header=T, sep=",", stringsAsFactors = FALSE) %>%
  rename(HYBAS_ID = PUID)

# Merge solutions with catchment shapefile
solution <- merge(basins %>% select(HYBAS_ID) , file_best, by = "HYBAS_ID")

# Set file name
solution_name <- paste0(wd, "/output/marxan/gis/", file_name, ".shp")

# Export solution shapefile
writeVector(solution,
            filename = solution_name,
            overwrite=TRUE)

# Print message
cat(paste(file_name, ": Marxan solution shapefile saved"))
