This repository contains code and data created in support of the project "Molecular and micro-architectural mapping of abnormal gray matter developmental trajectories in psychosis". All code was written in Matlab and R. Folders, files, and first steps are described below.

**Data**

The Data folder contains all the data required for running the analyses. Here are the files that need to be downloaded and stored in a specific location. The remaining files will be automatically generated:

-	The code used to compute centiles (BrainCharts folder) is available at https://github.com/brainchart/Lifespan

-	The all_microsc_DesikanKilliany68.csv file is available at https://github.com/netneurolab/netneurotools

-	The code used for the PCA-CCA analyses (cca_pls_toolkit_dev-grcca folder) is available at https://github.com/anaston/cca_pls_toolkit

First steps

1.	Download Code folder, which contains the scripts and functions used for the analyses.

2.	Download Data folder, which contains data used to run the analyses, and will include data derived from them.

3.	Centiles folder in Data must contain 4 subfolders: PE, FEP, SCZ_SAD, and groups. 
-	PE folder will contain data corresponding to Psychotic Experiences. Please store BrainCharts folder in it along with your .csv file containing volume data from these patients.

-	FEP folder will contain data corresponding to First Episode of Psychosis. Please store in it your .csv file containing volume data from these patients.

-	SCZ_SAD folder will contain data corresponding to relatives and chronic SCZ and SAD individuals. Please store in it your .csv file containing volume data from these patients.

-	groups folder will contain data from patients classified under different groups according to their clinical profile. Please store in it your .csv file containing volume data from these patients.


4.	Please store all_microsc_DesikanKilliany68.csv file and cca_pls_toolkit_dev-grcca folder in PCA_CCA folder.


Code
The Code folder contains all the code required for running the analyses and generate the figures. All scripts are designed to be sequentially executed. Don't forget to change the location variable regularly. 

-	NCP_1_means_PARC_aparc_volume.R – calculates the average of the different brain volume regions between both hemispheres, left and right. It takes a .csv file with volume data from both hemispheres as input, and returns a mean_volumes.csv file.

-	NCP_2_centile_estimation.R – calculates the centiles of brain volume regions from the mean_volumes.csv file. It returns a centiles.csv file.

-	NCP_3_centile_PE.m – classifies centiles by PE diagnosis (suspected, definite, or clinical). It storages the results in (1) a PE.mat variable that will be used for the rest of analyses, (2) a mean_PE-suspected.csv, mean_PE-definite.csv, and mean_PE-clinical.csv files that will be used for regional brain mapping (NCP_13_centile_mapping.R).

-	NCP_4_centile_SCZ_SAD.m – classifies centiles by SCZ or SAD diagnosis, or by relatives of SCZ or SAD patients (SCZ, SAD, SCZ-relatives, or SAD-relatives, respectively). It storages the results in (1) a SCZ.mat and Schizoaffective_Disorder.mat variables that will be used for the rest of analyses, (2) a mean_SCZ.csv, mean_SAD.csv, mean_SAD-relatives.csv, and mean_SCZ-relatives.csv that will be used for regional brain mapping (NCP_13_centile_mapping.R).

-	NCP_5_centile_FEP.m – selects centiles from the first session of a First Episode of Psychosis. It storages the results in (1) a FEP.mat variable that will be used for the rest of analyses, (2) a mean_FEP.csv file that will be used for regional brain mapping (NCP_13_centile_mapping.R).

-	NCP_6_centile_and_effect_sizes_groups.m – classifies centiles by groups according to their clinical profile: SCZ and SAD-relatives, PE, FEP or SCZ and SAD-chronic. It storages the results in (1) a G0.mat, G1.mat, G1_5.mat, and G2.mat variables, respectively, that will be used for the rest of analyses, (2) a mean_SCZ_and_SAD-relatives.csv, mean_PE.csv, mean_FEP.csv, and mean_SCZ_and_SAD-chronic.csv files that will be used for regional brain mapping (NCP_13_centile_mapping.R). Additionally, it computes the effect sizes between groups and saves the corresponding results in Data > Effect sizes.

-	NCP_7_region_permutation_test.m – it tests the regions that exhibit significant centile differences from healthy controls (Wilcoxon rank sum test), or that the effect sizes are significantly different from 0 (permutation test) with FDR correction. It takes as input all the .mat variables previously generated, and returns a significant_values.csv file that will be used for regional brain mapping (NCP_13_centile_mapping.R).

-	NCP_8_SSD.m – calculates the Sum of Squared Diferences (SSD) between diagnoses and groups (relatives, PE, FEP and chronic). It takes as input all the .mat variables previously generated, and displays two figures representing the SSD between diagnoses, and between groups.

-	NCP_9_parfor_groups_CCA.m – runs in parallel multiple PCA-CCA of psychosis-related groups classified according to their clinical profile. See NCP_10_CCA_cent_var.m and NCP_11_run_CCA.m scripts before execution.
 
-	NCP_10_CCA_cent_var.m – it is not a script to be executed, but for setting parameters by default for NCP_9_parfor_groups_CCA.m script. Default parameters include (see https://mlnl.github.io/cca_pls_toolkit/cfg/ for detailed information):

% Machine settings
o	cfg.machine.name = 'cca';
o	cfg.machine.metric = {'correl' 'trexvarx' 'trexvary'}; 
o	cfg.machine.param.name = {'VARx', 'VARy'}; % explained variance by the PCA components
o	cfg.machine.param.VARx = 0.6:0.1:0.9; % variance of data kept in the principal components during the SVD step of PCA-CCA  
o	cfg.machine.param.VARy = 1;   
o	cfg.machine.svd.varx = 1; % variance of X kept during the SVD step of PCA-CCA 
o	cfg.machine.svd.vary = 1; % variance of Y kept during the SVD step of PCA-CCA
o	cfg.machine.alignw = 'wX';

   	% Framework settings
o	cfg.frwork.name = 'permutation';     
o	cfg.frwork.split.nout % number of outer splits/folds
o	cfg.frwork.nlevel = 1;
    
% Deflation settings
o	cfg.defl.name = 'generalized'; 
    
% Environment settings
o	cfg.env.comp = 'local'; %  ['local', 'cluster']
o	cfg.env.save.tableHeading = {'set' 'varx' 'correl' 'pval' 'npcax'};
    
% Number of permutations
o	cfg.stat.nperm = 1000;
o	cfg.stat.nboot = 1000;

-	NCP_11_run_CCA.m – it is not a script to be executed, but for setting parameters by default for NCP_9_parfor_groups_CCA.m script. Default parameters include (see https://mlnl.github.io/cca_pls_toolkit/res/ for detailed information):

% Weight or loading visualization

o	res.gen.weight.type = 'correlation'; % for plotting Loadings instead of weights
o	%     res.behav.weight.norm = 'minmax'; % Normalize weight {'std','zscore'}
o	res.gen.weight.subset.type = 'significant'; % Filter weights by keeping only significant ones 
% {'top' (Filter weights by keeping top ones based on absolute value) 'minmax' (Filter weights by keeping top positive and negative ones)}
o	res.gen.weight.stat.ci = 'sd'; % {'pi'} Confidence interval
o	res.behav.weight.errorbar.side = 'one';

-	NCP_12_parfor_dx_CCA.m – is equivalent to NCP_9_parfor_groups_CCA.m script, but for individual diagnoses. It runs in parallel multiple PCA-CCA of psychosis-related diagnoses. See NCP_10_CCA_cent_var.m and NCP_11_run_CCA.m scripts before execution.

-	NCP_13_centile_mapping.R – generates the regional brain maps of: (1) predicted centiles/effect sizes, (2) empirical centiles/effect sizes, and (3) molecular and micro-architectural features from the .csv files previously generated.

-	NCP_14_join_images.m – generates the combined regional brain maps of empirical centile/effect sizes with each regional brain map of their respective (1) significant regions and (2) predicted centiles/effect sizes.

-	NCP_15_join_3_images.m – generates the combined regional brain maps of empirical centile/effect sizes with each regional brain map of their respective (1) significant regions and (2) predicted centiles/effect sizes into a single figure.

-	NCP_16_figures_CCA_loadings.m – generates the figures corresponding to the significant loadings and weights of PCA-CCA modeling.

-	NCP_17_figures_correlation.m – generates the figures corresponding to the (1) neurobiological similarity matrix, (2) co-vulnerability to psychosis matrix, and (3) their respective correlation.


Function calls
This section contains the functions that are essential for running the scripts but must not be executed.

-	computeCohen_d.m – computes the Cohen’s distance between two vectors. It is called by NCP_6_centile_and_effect_sizes_groups.m script.

-	mix_dx.m – creates randomized groups by mixing patients with different diagnoses or group membership. It is called by NCP_7_region_permutation_test.m and NCP_8_SSD.m scripts. It was also called by NCP_9_parfor_groups_CCA.m and NCP_12_parfor_dx_CCA.m scripts when doing the leave-one-study-out cross-validation.

