#################################################################
##       Download environmental data and crop them             ##
#################################################################
#### Working directory and packages ####
wd <- "C:/Users/julie/Documents/macroclim/data/"
setwd(wd) 

pacman::p_load(dplyr, ncdf4, fst, #data wrangling
               terra, tidyterra, # gis
               httr, XML # to work with online content (download data)
)

bbox_terra <- c(-10.5,32.702, 34.856, 71.31)

##--------------------------------------------------
##          Download environmental files        --
##--------------------------------------------------
#### Download CHELSA variables: present and future ####
# Create download directory
dir.create(paste0(wd, "/chelsa_bioclim"))

# Extend timeout to allow uninterrupted downloading
options(timeout = 1000000)

# Download
chelsa_var=c(paste0("bio0", 1:9), paste0("bio", 10:19), # "traditionnal" bioclimatic variables
             "scd" #, snow cover days
             # "cmi_max","cmi_mean","cmi_min","cmi_range" #, #moisture but no values for future, we take this variable in HydroRIVERS
             # "hurs_max","hurs_mean","hurs_min","hurs_range" #surface relative humidity
             ) # not integrating koppen geiger climate classification, vapor pressure, evapotranspiration, cloud cover, growing degree days

# Parameters chosen for future prediction
period = "2041-2070" # [ choices are 2011-2040, 2041-2070, 2071-2100]
gcm = "IPSL-CM6A-LR" # [ choices are GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL]
ssp = "126" # [choices are rcp 126, rcp 370, rcp 585]

url_chelsa <- "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/"

for (i in chelsa_var){
  if (!file.exists(paste0("/chelsa_bioclim/",i,"_1981-2010.tif"))){ # if not already downloaded
    #Present
    download.file(paste0(url_chelsa, i, "/1981-2010/CHELSA_",i,"_1981-2010_V.2.1.tif"),
                  destfile = paste0("chelsa_bioclim/",i,"_1981-2010.tif"), mode = "wb")
    
    # Future
    download.file(paste0(url_chelsa, i, "/", period, "/", gcm, "/ssp", ssp, "/CHELSA_",
                         tolower(gcm), "_ssp", ssp, "_", i,"_", period, "_V.2.1.tif"),
                  destfile = paste0("chelsa_bioclim/", i, "_", tolower(gcm), "_ssp", ssp, "_", period, ".tif"), mode = "wb")
    }
}
Sys.time()

#### Download FutureStreams variables: present and future ####
dir.create(paste0(wd, "futurestreams"))

url_FS <- "https://geo.public.data.uu.nl/vault-futurestreams/research-futurestreams%5B1633685642%5D/original/indicators/"

# select files to download
var_names_FS=c("waterTemp", "discharge")
scenario_FS=c("hist","rcp2p6","rcp4p5","rcp6p0","rcp8p5")
decade_FS=c("2021-01-31_to_2040-12-31","2041-01-31_to_2060-12-31","2061-01-31_to_2080-12-31","2081-01-31_to_2099-12-31")

for (v in var_names_FS){ # loop on variables (discharge or temperature)
  
  url_FS_var <- paste(url_FS,v,sep="")
  indicators_FS <- readHTMLTable(content(GET(url_FS_var), "text"))[[1]]$Name
  indicators_FS <- indicators_FS[!grepl("Store", indicators_FS)]
  
  for (i in indicators_FS){  # loop on indicators for this variable
    
    if(!(i %in% c("Q-bfi/", "Q-hvi/"))){
    url_FS_var_ind=paste(url_FS_var,"/",i,sep="")
    
      # for (s in scenario_FS) {  # loop on scenarios
      for (s in c("hist","rcp2p6")){
        # s = scenario_FS[1] # scenario = hist
        url_FS_var_ind_sce=paste(url_FS_var_ind, s, sep="")
        GCM_FS=readHTMLTable(content(GET(url_FS_var_ind_sce), "text"))[[1]]$Name # select a  General Circulation Model for meteo data
        
        # for (g in GCM_FS){  # loop on GCMs
          g = GCM_FS[4] # GCM = ipsl
          url_FS_var_ind_sce_GC=paste(url_FS_var_ind_sce,"/", g, sep="")
          
          if (s != "hist") {
            # for (d in decade_FS) { # loop on available decades
            d = decade_FS[2] # decade 2040-2060
                print(paste(v, i, s, g, d))
              
              url_FS_to_download=paste0(url_FS_var_ind_sce_GC,gsub('.{1}$', '', i),"_",
                                        gsub('.{1}$', '',g),"_",s,"_",
                                        d,".nc")
              destfile_FS <- paste0(wd, "/futurestreams/",gsub('.{1}$', '', i),"_",
                                    gsub('.{1}$', '',g),"_",s,"_",d,".nc")
              if (!file.exists(destfile_FS)){
                download.file(url_FS_to_download, destfile_FS, mode = "wb")
              }
  
              # }
            } else {
              print(paste(v, i, s, g))
              
              url_FS_to_download=paste0(url_FS_var_ind_sce_GC,gsub('.{1}$', '', i),"_",
                                        gsub('.{1}$', '',g),"_hist_1976-01-31_to_2005-12-31.nc")
              destfile_FS <- paste0(wd, "/futurestreams/",gsub('.{1}$', '', i),"_",
                                    gsub('.{1}$', '',g),"_hist_1976-01-31_to_2005-12-31.nc")
              if (!file.exists(destfile_FS)){
                download.file(url_FS_to_download, destfile_FS, mode = "wb")
              }
            }
            
        # } 
      }
    }
    
  }

  rm(list=c("destfile_FS","url_FS_to_download"))
}

#### Download LUCAS LUC land cover: present and future####

# the download requires to create an account and manually select the files


# https://www.wdc-climate.de/ui/entry?acronym=LUC_hist_EU_afts_v1.1
# https://www.wdc-climate.de/ui/entry?acronym=LUC_future_EU_ssp126_v1.1
# https://www.wdc-climate.de/ui/entry?acronym=LUC_future_EU_ssp585_v1.1


#### Download flow intermittence data ####

# Manually downloaded at https://figshare.com/articles/journal_contribution/Streamflow_intermittence_in_European_reaches/24591807

#### Download HydroATLAS ####
dir.create(paste0(wd, "/hydroatlas/"))

# manual downloads at https://www.hydrosheds.org/hydroatlas#download

##--------------------------------------------------
##   Crop the environmental data to study extent  --
##--------------------------------------------------

#### Define output directory for merged files #### 
layer_dir <- paste0(wd, "final_layers")
# Create the directory if it has not already been created before
if (!dir.exists(layer_dir)) {
  dir.create(layer_dir)
}

### CHELSA ------
# directory containing the layers to be cropped
dirs_chelsa <- paste0(wd, "/chelsa_bioclim")
# crop the files using the function crop_to_extent() in a loop
files_chelsa <- list.files(dirs_chelsa, pattern = ".tif", full.names = TRUE)

for(ifile in files_chelsa) {
  if (!file.exists(paste0(layer_dir,"/",basename(ifile)))){ # if crop file does not already exists
    # crop_to_extent(
    #   raster_layer = ifile,
    #   bounding_box = bbox_hydrographr,
    #   out_dir = layer_dir,
    #   file_name = basename(ifile),
    #   read = FALSE,
    #   quiet = TRUE)
    raster_ch<- raster::raster(ifile)
    terra::crop(x = raster_ch,
                y = bbox_terra,
                file.path(layer_dir, basename(ifile)),
                todisk=TRUE,
                overwrite=TRUE)
  }
}

### FutureStreams ------
# directory containing the layers to be cropped
dirs_fs <- paste0(wd, "futurestreams")
# crop the files using the function crop_to_extent() in a loop
files_fs <- list.files(dirs_fs, pattern = "\\.nc$", full.names = TRUE)

for(ifile in files_fs) {
  # print(ifile)
  if (!file.exists(paste0(layer_dir,"/",paste0(str_remove(basename(ifile), ".nc"), ".tif")))){ # if crop file does not already exists
    # print("cropping")
    raster_fs<- raster::raster(ifile)
    terra::crop(x = raster_fs,
                y = bbox_terra,
                file.path(layer_dir, paste0(str_remove(basename(ifile), ".nc"), ".tif")),
                todisk=TRUE,
                overwrite=TRUE)
    rm(raster_fs); gc()
  }
}

### HydroATLAS ------
atlas <-terra::vect(paste0(wd,"/hydroatlas/RiverATLAS_v10_shp/RiverATLAS_v10_eu.shp"),
                    extent=bbox_terra) #already cropped to the right extent

# names(atlas)[order(names(atlas))] # see RiverATLAS_Catalog_v10.pdf downloaded with RiverATLAS for details

### flow intermittence  ------
dirs_dryver <- paste0(wd, "dryver")

dryver_hyriv  <- read.csv(paste0(dirs_dryver,"/hydrorivers_dryver_link.csv"),
                 header=T, sep=",", stringsAsFactors = FALSE)

dryver_hist_name <- file.path(dirs_dryver, "reference1985-2014", "dryver_net_irclasses_ipsl_eu.fst")
dryver_hist <- read_fst(dryver_hist_name)%>%
  select(-c(from_node:upa)) %>%
  mutate(jan = rowMeans(select(., contains("-01-")), na.rm = T)) %>%
  mutate(feb = rowMeans(select(., contains("-02-")), na.rm = T)) %>%
  mutate(mar = rowMeans(select(., contains("-03-")), na.rm = T)) %>%
  mutate(apr = rowMeans(select(., contains("-04-")), na.rm = T)) %>%
  mutate(may = rowMeans(select(., contains("-05-")), na.rm = T)) %>%
  mutate(jun = rowMeans(select(., contains("-06-")), na.rm = T)) %>%
  mutate(jul = rowMeans(select(., contains("-07-")), na.rm = T)) %>%
  mutate(aug = rowMeans(select(., contains("-08-")), na.rm = T)) %>%
  mutate(sep = rowMeans(select(., contains("-09-")), na.rm = T)) %>%
  mutate(oct = rowMeans(select(., contains("-10-")), na.rm = T)) %>%
  mutate(nov = rowMeans(select(., contains("-11-")), na.rm = T)) %>%
  mutate(dec = rowMeans(select(., contains("-12-")), na.rm = T)) %>%
  select(-(matches("19|200|201"))) %>%
  mutate(yr_irclass_present = rowMeans(select(., -matches("DRYVER_RIV")), na.rm = T)) %>%
  right_join(dryver_hyriv, by = "DRYVER_RIV") %>%
  select(HYRIV_ID, yr_irclass_present) %>%
  arrange(desc(yr_irclass_present)) %>%
  group_by(HYRIV_ID) %>%
  slice(1) %>%
  ungroup()

dryver_future_name <- file.path(dirs_dryver, "future_periods", "RCP2.6", "2041-2070", "dryver_net_irclasses_ipsl_eu.fst")
dryver <- read_fst(dryver_future_name) %>%
  select(-c(from_node:upa)) %>%
  mutate(jan = rowMeans(select(., contains("-01-")), na.rm = T)) %>%
  mutate(feb = rowMeans(select(., contains("-02-")), na.rm = T)) %>%
  mutate(mar = rowMeans(select(., contains("-03-")), na.rm = T)) %>%
  mutate(apr = rowMeans(select(., contains("-04-")), na.rm = T)) %>%
  mutate(may = rowMeans(select(., contains("-05-")), na.rm = T)) %>%
  mutate(jun = rowMeans(select(., contains("-06-")), na.rm = T)) %>%
  mutate(jul = rowMeans(select(., contains("-07-")), na.rm = T)) %>%
  mutate(aug = rowMeans(select(., contains("-08-")), na.rm = T)) %>%
  mutate(sep = rowMeans(select(., contains("-09-")), na.rm = T)) %>%
  mutate(oct = rowMeans(select(., contains("-10-")), na.rm = T)) %>%
  mutate(nov = rowMeans(select(., contains("-11-")), na.rm = T)) %>%
  mutate(dec = rowMeans(select(., contains("-12-")), na.rm = T)) %>%
  select(-(matches("19|200|201"))) %>%
  mutate(yr_irclass_future = rowMeans(select(., -matches("DRYVER_RIV")), na.rm = T)) %>%
  right_join(dryver_hyriv, by = "DRYVER_RIV") %>%
  select(HYRIV_ID, yr_irclass_future) %>%
  arrange(desc(yr_irclass_future)) %>%
  group_by(HYRIV_ID) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(dryver_hist, by = "HYRIV_ID") %>%
  filter(HYRIV_ID %in% atlas$HYRIV_ID)

write.csv(dryver,
          paste0(dirs_dryver,"dryver_pres_fut.csv"), row.names = FALSE)

### LUCAS LUC ------
# directory containing the layers to be cropped
dirs_ll<- paste0(wd,"lucas_luc")

# crop the files using the function crop_to_extent() in a loop
files_ll <- c(paste0(dirs_ll,"/LUC_future_EU_ssp126_v1.1_1-9/LUCAS_LUC_v1.1_ssp126_Europe_0.1deg_2046_2055.nc"),
              paste0(dirs_ll,"/LUC_hist_EU_afts_v1.1_1-7/LUCAS_LUC_v1.1_historical_Europe_0.1deg_2010_2015.nc"))

for(ifile in files_ll) {
  raster_ll <- terra::rast(ifile)
  # there are separate layers for each land cover type and each year
  # we take the last year, there are 16 land cover types
  n <- dim(raster_ll)[3]
  raster_ll <- subset(raster_ll, (n-15):n)
  for (i in c(3:5,7:9,13:15)) { # selection of land cover types, see https://essd.copernicus.org/preprints/essd-2021-252/essd-2021-252.pdf
    raster_lli <- subset(raster_ll,i) # select one layer = one land cover type
    time(raster_lli) <- NULL # time not necessary because we keep only 1 year and it creates a .json file if we keep it
    writeRaster(terra::crop(x = raster_lli,
                                y = bbox_terra),
                paste0(layer_dir,"/",str_remove(basename(ifile), ".nc"),"_lc",i, ".tif"),
                overwrite=TRUE)
    rm(raster_lli)
  }
  rm(raster_ll) # to clean up memory
}

gc()