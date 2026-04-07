#### Working environment ####
wd <- paste0(getwd(), "/data/")
setwd(wd)
pacman::p_load(dplyr, tidyterra, terra) # Data wrangling and map plots

# Create folder for the output shapefiles of MARXAN if needed
if (!dir.exists(paste0(wd, "output/marxan/gis/"))){
  dir.create(paste0(wd, "output/marxan/gis/"), recursive=T)
}

#### Set parameters ####

discard <- "no" # no or 10 or 50
frag = "nofrag" # "" o "_frag"

file_name <- paste0(discard, "discard_nofrag")
file_name_frag <- paste0(discard, "discard_frag")
niche_folder <- prediction_folder <- paste0(discard, "discard/")

print(paste("___________________", file_name, "___________________" ))


#### Load hydroatlas data ####
if (! ("basins" %in% ls(envir = .GlobalEnv))) {
  bbox_terra <- c(-10.5,32.702, 34.856, 71.31)
  basins <- vect(paste0(wd,
                        "/environment/hydroatlas/hybas_eu_lev01-12_v1c/hybas_eu_lev12_v1c.shp"),
                 extent=bbox_terra)
}

#### Run MARXAN ####
wd_marxan <- paste0(wd,"output/marxan/", file_name)
setwd(wd_marxan)

file.copy(paste0(wd, "/marxan/Marxan_x64.exe"),
          wd_marxan)
file.copy(paste0(wd, "/marxan/input.dat"),
          wd_marxan)

dir.create(paste0(wd_marxan, "/output"), recursive=T)

system("Marxan_x64.exe", wait = T, invisible = T) # Call Marxan to execute

#### export  shapefile ####

file_best <- read.csv(paste0(wd_marxan, "/output/output_best.csv"), 
                      header=T, sep=",", stringsAsFactors = FALSE) %>%
  rename(HYBAS_ID = PUID)

solution <- merge(basins %>% select(HYBAS_ID) , file_best, by = "HYBAS_ID")

solution_name <- paste0(wd, "/output/marxan/gis/", file_name, ".shp")

writeVector(solution,
            filename = solution_name,
            overwrite=TRUE)

cat(paste(file_name, ": Marxan solution shapefile saved"))
