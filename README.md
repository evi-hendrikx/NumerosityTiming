# NumerosityTiming

These scripts were used to run the analysis for

Hendrikx, E., Paul, J.M., van Ackooij, M., van der Stoep, N., & Harvey, B.M. Cortical quantity representations of visual numerosity and timing overlap increasingly into superior cortices but remain distinct.

## Running the pipeline
### Creating the dataset
NumTim_run_all calls other scripts you need. You still need to provide paths where the subjects' information is stored. It initially creates datasets with the required parameters for all the analyses (NumTim_load_data and NumTim_load_time_series). Note that the datasets created here for our data are available at the following DOIs: map parameters (10.6084/m9.figshare.24126633) & time series (10.6084/m9.figshare.24126663). It then re-combines the data according to requested voxel selections (NumTim_restructure_info & NumTim_get_select_ids & NumTim_prepare_whichCombi); performs the analyses (except Bayesian topology stats, which are done in R) (NumTime_chance_anova_stats & NumTim_compare_proportions & NumTime_compare_correlations & NumTime_compare_topo); and plots figures of the results (NumTim_plots_prop_corr & NumTim_plots_topo & NumTim_compare_proportions & NumTime_compare_correlations). 

### Creating meshes
As the meshes require 3d information about the voxel coordinates, this information cannot be extracted for the previously created datasets. They are created separately (surface_plots or manually).

### Bayesian topology statistics
Performed with R_bayesian_github.R

### Other required folders
PyColormap4Matlab (can be found here: https://nl.mathworks.com/matlabcentral/fileexchange/68239-pycolormap4matlab)
circstat-matlab-master (can be found here: https://github.com/circstat/circstat-matlab)
For analyses that calculate relations with the coordinates of the ROIs, files containing the N27 coordinates can be found here: 10.6084/m9.figshare.24147252

## Licence
Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## Contact
For any questions or remarks, please send an email to e.h.h.hendrikx@uu.nl
