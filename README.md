Code associated with the manuscript "A roadmap for monitoring climate change effects on river biodiversity" from Crabot et al.
All scripts are available in the /code folder, files that were not too heavy for GitHub were made available in the /data folder.
As a large part of the biological data is private and not ours to share, and in order to have maller files that make the scripts
faster to run, only a very small subset of the biological data of Haase et al. (2023) is used in this repository.

Haase, P. et al. The recovery of European freshwater biodiversity has come to a halt. Nature 620, 582–588 (2023).

In the /code folder, scripts are numbered in the order in which they should be used to follow the procedure outlined in the paper.
All scripts are commented but here is the general structure of the procedure:

- 1.Running_All_Models: This script contains the code that wraps up the whole process, setting the working environment, loading data
and calling all the scripts from 2.0_SDM_and_distribution_descriptors to 6.Calculate_vulnerability_scores.

- 2.0_SDM_and_distribution_descriptors: This script is used to run the species distribution models (SDM). It is a loop running on all
taxa. The steps of the SDM have been broken down into different scripts from 2.1 to 2.4 (see below).
   - 2.1_Biological_data_wrangling: very short script to select the biological data related to the selected taxon. Available lines of
   code also allows plotting the distribution map of the observations.
   - 2.2_Environmental_data_wrangling: Setting up the environmental variables for the SDMs if not already loaded.
   - 2.3_SpeciesDistributionModels: Finalizing the formatting of the training dataset, running the random forest and projecting
   future habitat suitability. Possibility to export distribution maps.
   - 2.4_Extract_model_information: Extracting information on the model (variable importance, out-of-bag-error), on the distribution
   for present and future (latitudinal range, altitudinal range...) and compiling everything in a table ("distrib_summary").
  
- 3.Calculate_distribution_change: Calculating the distribution change (e.g. difference in maximal altitude between present and 
future) for each taxa and compiling everything in a table ("taxa_selection_parameters").
- 4.New_predictions_with_fragmentation: Scenario considering fragmentation by dams. For all non-flying taxa, it considers all 
"colonization" events (increasing habitat suitability between present and future) and removes new occurrences if they occur
upstream of a dam. Concretely, for a given genus, it alters the prediction file for future by setting the occurrence from 1 to 0.
- 5.Calculate_distribution_change_after_fragmentation: script similar to the script 3. but considering the changes brought by script 4.
- 6.Calculate_vulnerability_scores: for a given scenario (e.g. without discarding impacted sites, without considering fragmentation
by dams), selecting good indicator taxa and inputing a Macroclim score to them.

- 7.Build_MARXAN_files: MARXAN is the prioritization software. It requires a very specific formatting of the files and requires at
least three files: the community matrix in long format, a list of the planning units (here the catchments) with the cost of their
monitoring, and the list of the taxa with the proportion of their distribution area that should be targeted. This script builds
these files.
- 8.Plot_MARXAN_priority_maps: Running MARXAN and saving the map of the best solution.

/other/ folder :
   - Calculating_network_distance.R: Calculating pairwise distance between reach outlets along the river network of HydroRIVERS.
   - Calculating_temporal_betadiversity.R: Calculating temporal beta diversity using the vegdist function on the predicted
     occurrence for present and future scenario.
   - Download_environmental_data_and_crop.R: Indications on the download of all environmental variables that stem from RiverATLAS,
     CHELSA, FutureStreams, LUCAS LUC, and Abbasi et al. 2025. Compilation of all data at the river reach scale.
   - Tuning_RF_parameters.R: Running the random forest models as in 2.3_SpeciesDistributionModels.R but checking for different
     possible combination of parameters and assessing model performance.
