# NumerosityTiming

These scripts were used to run the analysis for

Hendrikx, E., Paul, J.M., van Ackooij, M., van der Stoep, N., & Harvey, B.M. Cortical quantity representations of visual numerosity and timing increasingly overlap into superior cortices but remain distinct.

## Running the pipeline
### Creating the dataset
NumTim_run_all calls other scripts you need. You still need to provide paths where the subjects' information is stored. It initially creates datasets with the required parameters for all the analyses (NumTim_load_data and NumTim_load_time_series). Note that the datasets created here for our data are available at the following DOIs: map parameters (DOI) & time series (DOI). It then re-combines the data according to requested voxel selections (NumTim_restructure_info & NumTim_get_select_ids & NumTim_prepare_whichCibi); performs the analyses (except the ones concerning repeated measures) (NmTime_chance_anova_stats & NumTim_compare_proportions & NumTime_compare_correlations); and plots figures of the results (NumTim_plots_prop_corr & NumTim_plots_topo & NumTim_compare_proportions & NumTime_compare_correlations)

### Creating meshes
As the meshes require 3d information about the voxel coordinates, this information cannot be extracted for the previously created datasets. They are created separately (surface_plots or manually).

### Contact
For any questions or remarks, please send an email to e.h.h.hendrikx@uu.nl
