####  Selecting the site locations of the selected taxon ####

species <- bio %>%
  filter(taxon == taxa_selected & abundance > 0) %>%
  select(site_id, date, abundance) %>%
  mutate(id = paste(site_id, date, sep = "_"))

spdata <- species %>%
  left_join(sites %>% select(site_id, longitude, latitude), by = "site_id") %>%
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude))

# ### To plot the sites where the taxon was observed
# library(leaflet)
# bbox <- c(-10.5, 34.856,  32.702, 71.31)
# map <- leaflet() %>%
#   addProviderTiles('Esri.WorldShadedRelief') %>%
#   setMaxBounds(bbox[1], bbox[2], bbox[3], bbox[4]) %>%
#   addCircleMarkers(data = spdata,
#                    color = "#0066CC",
#                    radius = 0.5)
# map
