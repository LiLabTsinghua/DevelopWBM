# Human2 manuscript data and scripts
## This directory contains the scripts and data files used to perform analyses and generate models and figures that were used in the Human2 manuscript.

### Required Software:
* A functional Matlab installation (MATLAB 9.6.0 (R2019a) and higher)
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN) (version 2.9.0)
* The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox) (version 3.4)
* The [GECKO toolbox for MATLAB](https://github.com/opencobra/cobratoolbox) (version 3.2.2)
* The Human-GEM GitHub repository (version 2.0.0)
    * Ensure that the directories of this repository have been added to the MATLAB path
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for MATLAB (version 10.0.2)
* For generating some figures, a functional [R installation](https://www.r-project.org/) (version 3.6.1)


## Generation of Figure 2 (b,c,d,e)
* Run the `Figure2b_geneCompartment.py`, `Figure2c_MainContents.py` and `Figure2e_GPRRatio.py` script (in the `figures` subdirectory) to generate the different panels in Figure 2 of the main text.

## Validation of Human2
* Validation of essential genes. Run the `Run_GeneEssentiality.m` script to get results of essential gene analysis.
* Running the `IEM_simulation.py` script to get the results of inborn-error metabolism simulation. Then run the `Figure3c_IEM.py` script to generate the panel in Figure 3 of the main text.
* Running the `FVA_FluxDistribution.m` function to get the different fluxes distributions and the panel in Support-figure 3 of the main text.
* To use the ecGEMs and non-ecGEMs to predict growth rates and metabolite exchange rates with increasing levels of constraints (as shown in XXX), run the `predict_cellLines_gRates.m`.

## Generation and comparation of Organ-specific models
* The `Gen_Organ_models.m` script is used to generate organ-specific models. There are five types of individuals available for selection, including AdultMale, AdultFemale, ElderlyMale, ElderlyFemale and Fetus.
* Running the `Run_Comparation_OrganTasks.m` can check and compare the specific organ metabolic tasks that the generated organ-specific models can perform. Then run the `Fig_4d_GEM_tasks.py` script to generate the panel in Figure 4d of the main text.
* Running the `ModelComparation.m` script to compare the subsystem coverage, similarity, function and sentivity of all organ-specific GEMs. Then run the `Sup_Fig_6a_SubsystemCoverage.py` and `Sup_Fig_6b_metabolicTasks.py` script to generate the panel in Support Figure 6.
* The Sensitivity analysis of organ-specific GEMS can refer to (https://gitlab.com/csb.ethz/functionalcomparisonmetabnetworks/-/tree/main). The results were stored in './results/Sensitivity/', and run `Sup_Fig_7_OrganModel_sensitivity.py` script to generate the panel in Support Figure 7.

## Generation and simulation of whole-body models
* Running the `runMergeModels.m` function can generate WBM for different age and gender groups. Then run the `Fig_4e_WBM_contents.py` script to generate the panel in Figure 4e of the main text.
* Running the `WBMSimilarity.m` function can get similarity among WBMs. Then run the `Fig_4f_WBM_similarity.py` script to generate the panel in Figure 4f of the main text.
* Running the `runBMR.m` function can get the simulated results of basic metabolic rate and it could be plotted with running 'Sup_Fig_8a-d.py'. Further more, Running the `BMR_deleteGenes.m` could get all essential genes which could effect the BMR.
* Running the `PrepareSHAPvalue.m` can obtain the BMR of adult males under different physiological parameters. Furthermore, running `Sup_Fig_8e_BMR_SHAPValue.py` can conduct a SHAP analysis and obtain Support figure 8e.
* Running the `run_AgingMets.m` can obtain the changes of blood and urine metabolites as age increases. Then run the `Sup_Fig_9b-d_Aging_mets.py` script to generate the panel in Figure 9.

## Generation and simulation of CoreWBM
* Running the `runCoreWBM.m` function will generate coreWBM and different food simualtion results and it could be plotted with running `Fig_5a_Food.py`.
* Running the `dFBA.m`  function will generate the different panels in `Figure 6` of the main text.





