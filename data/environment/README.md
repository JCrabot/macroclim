Environmental data to test the code associated with the manuscript "A roadmap for
monitoring climate change effects on river biodiversity" from Crabot et al. 

amber-project-data/ folder:
   - atlas.csv: raw data from the AMBER project, with information on European
   barriers to dispersal, from https://amber.international/european-barrier-atlas/
   - intersect_amber_dams_hydroshed.csv: result of the snapping of dams from 
   AMBER to the HydroRIVERS network
   
hydroatlas/ folder:
   - The three subfolders contain files of HydroRIVERS, RiverATLAS and HydroBASINS,
   downloaded from https://www.hydrosheds.org/hydroatlas
   - hyriv_hybas_equivalence.csv: gives the equivalence between ID of HydroRIVERS
   reaches and ID of catchments of HydroBASINS

environment.csv: This file combines all environmental variables that stem from
RiverATLAS, CHELSA, FutureStreams, LUCAS LUC, and Abbasi et al. 2025. They were
associated to the HydroRIVERS river reaches as in the chunk "Associate all
environmental variables with HydroRIVERS network" of the R script in the 
code/other folder "Download_environmental_data_and_crop.R".

list_polluted_segments_10.csv: List of the IDs of the river reaches that should 
be discarded if "discard" is set to 10 (meaning that 10% of the most impacted
sites should be excluded when running the distribution lmodels). This is the
output of the chunk "Extract HYDROSHED ID of polluted river segments if
necessary" in the R script "2.2_Environmental_data_wrangling.R"

list_polluted_segments_50.csv: List of the IDs of the river reaches that should 
be discarded if "discard" is set to 50 (meaning that 50% of the most impacted
sites should be excluded when running the distribution lmodels). This is the
output of the chunk "Extract HYDROSHED ID of polluted river segments if
necessary" in the R script "2.2_Environmental_data_wrangling.R"

network_length.csv: pairwise distance between reach outlets along the river
network of HydroRIVERS as calculated in the R script in the code/other folder
"Calculating_network_distance.R"

