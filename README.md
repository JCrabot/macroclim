Code associated with the manuscript "A roadmap for monitoring climate change effects on river biodiversity" from Crabot et al.
All scripts are available in the /code folder, files that were not too heavy for GitHub were made available in the /data folder.
As a large part of the biological data is private and not ours to share, and in order to have maller files that make the scripts
faster to run, only a small subset of the biological data of Haase et al. (2023) is used in this repository.

Haase, P. et al. The recovery of European freshwater biodiversity has come to a halt. Nature 620, 582–588 (2023).

In the /code folder, scripts are numbered in the order in which they should be used to follow the procedure outlined in the paper.
All scripts are commented but here is the general structure of the procedure:

0_Download_environmental_data_and_crop.R: 

1.Running_All_Models:

2.0_SDM_and_distribution_descriptors:
  2.1_Biological_data_wrangling:
  2.2_Environmental_data_wrangling:
  2.3_SpeciesDistributionModels:
  2.4_Extract_model_information:
  
3.Calculate_distribution_change:
4.New_predictions_with_fragmentation:
5.Calculate_distribution_change_after_fragmentation
6.Calculate_vulnerability_scores

7.Build_MARXAN_files
8.Plot_MARXAN_priority_maps
