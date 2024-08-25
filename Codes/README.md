# Human2 manuscript data and scripts
This directory contains the scripts and data files used to perform analyses and generate models and figures that were used in the Human2 manuscript.

### Required Software:
* A functional Matlab installation (MATLAB 9.6.0 (R2019a) and higher)
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN) (version 2.8.2)
* The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox) (version 3.4)
* The [GECKO toolbox for MATLAB](https://github.com/opencobra/cobratoolbox) (version 3.2.1)
* The Human-GEM GitHub repository (version 2.0.0)
    * Ensure that the directories of this repository have been added to the MATLAB path
* The [Sensitivity analysis for matlab] (https://gitlab.com/csb.ethz/functionalcomparisonmetabnetworks/-/tree/main)
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for MATLAB (version 10.0.2)
* For generating some figures, a functional [R installation](https://www.r-project.org/) (version 3.6.1)


## Generation of Figure 2 (b,c,e)
* Run the `Figure2b_geneCompartment.py`, `Figure2c_MainContents.py` and `Figure2e_GPRRatio.py` script (in the `figures` subdirectory) to generate the different panels in Figure 2 of the main text.

## Validation of generatic models
* Run the `IEM_simulation.py` script (in the `figures` subdirectory) to get the results of inborn-error metabolism simulation and Run the `Figure3c_IEM.py` to generate the panel in Figure 3 of the main text.
* Run the `Figure3a_GeneEssential.py` script (in the `figures` subdirectory) to obtain the plot in Figure3 (in the `figures` subdirectory) after getting the results of "Generation of all ftINIT models" and "Gene essentiality analysis for different Human-GEM versions" below.
* Run the `runFluxDistribution.m` script  to get the different fluxes distribution and run `Figure3b_FluxesDistribution.py` generate the panel in Figure 3 of the main text.

## Generation of all ftINIT models
* The `gen_all_ftINIT_models.m` script contains the code used to generate the ftINIT GEMs for cell lines and organs presented in the manuscript.
* **NOTE:** This function should be run on a compute cluster as it will take a long time to run.

## Gene essentiality analysis for different Human-GEM versions
* The first part of `gen_models_get_essentials.m` function will generate the 5 cell-line GEMs corresponding to the Hart 2015 study and evaluate gene essentiality for any given version(s) of Human-GEM.
* The second part of `gen_models_get_essentials.m` function will generate the 625 cell-line GEMs corresponding to the Depmap database and evaluate gene essentiality for any given version(s) of Human-GEM.
* The third part of `gen_models_get_essentials.m` function will generate different age and gender organ GEMs as base models for whole-body models(WBM) 

## Comparation of organ-specific GEMs
* Run the `compareMultipleModels.m` script to compare the subsystem coverage, similarity, function and sentivity of all organ-specific GEMs.
* Compared results could be plotted after runing 'SupportFigure6_SubsystemCover.py', 'SupportFigure7_FullTasks.py', 'SupportFigure8_Sentivity.R' script (in the `figures` subdirectory) to generate the different panels in Support Figures.

## Comparation of organ-specific GEMs
* Run the `runMergeModels.m` function will generate WBM for different age and gender group.
* Run the `runBMR.m` function will get the simulated results of basic metabolic rate and it could be plotted with running 'Figure4f_BMR.py'.
* Run the `PrepareSHAPvalue.m` function will get the simulated BMR with different parameters and it could be analysised and plotted with running 'SHAP_BMR.py'.

## Generation and simulation of CoreWBM
* Run the `runCoreWBM.m` function will generate coreWBM and different food simualtion results and it could be plotted with running 'Figure5a_FoodCorrelation.py'.
* Run the `dFBA_Fasting.m`  function will generate the different panels in Figure 6 of the main text.





