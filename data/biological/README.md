Biological data to test the code associated with the manuscript "A roadmap for
monitoring climate change effects on river biodiversity" from Crabot et al. Only
a small subset of the biological data of Haase et al. (2023) is used in this
repository.

Haase, P. et al. The recovery of European freshwater biodiversity has come to a
halt. Nature 620, 582–588 (2023).

bio_genus.csv: Very small subset of data from Haase et al. 2023. Provides the 
location of 9 invertebrate genera on each site and for one or several dates.
Columns are site_id, date, season, taxon, abundance, resolution, and then the 
classification of each genus (phylum to genus). Details can be checked directly
in Haase et al. 2023.

classification_macroclim.csv: Only the classification (phylum to genus) of each
taxon

datasets_by_genus.csv: Datasets to use in the random forest models for each given
taxon.

RF_parameters.csv: Parameters to use for the random forest, which stems from a 
parameter tuning script. The function used for the parameter tuning is available in
the script "Tuning_RF_parameters.R" in the code/other/ folder. The same function
was used than in the "2.3_SpeciesDistributionModels.R" script but this time
we checked for different values of "min.node.size", "min.bucket" and "max.depth"
and compared the model performance using different indicators but mostly the
True Skill Statistic (TSS).

sites.csv: Details on each sampling site as in Haase et al. 2023.
Columns "site_id" to "resolution" provide information on the country, associated 
dataset, details on sampling method, main taxonomic resolution reached and spatial
coordinates. Additional columns (from HYRIV_ID to HYBAS_L12) stem from the
snapping of the sampling sites to the HydroRIVER network and are the 
matching columns of the HydroRIVERS shapefile.

summary_occurences.csv: provides for each genus the total number of occurrences,
the frequency of occurrence across sites and the number of countries where it was
observed.

temporal_turnover.csv: temporal beta diversity calculated as in the script 
in the code/other/ folder "Calculating_temporal_betadiversity.R", using the 
vegdist function on the predicted occurrence for present and future scenario.
