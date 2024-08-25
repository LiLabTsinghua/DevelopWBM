% This script is used to merge organ models into whole-body model.
% The supported options are 'male', 'female', and 'fetus'. Please
% define those using the variable 'sex' (e.g., sex = 'male').

Age = 'Adult';
Sex = 'Female';
modelPath = ['../Models/models', '_', Age, Sex];
load(modelPath);
% models = load(modelPath);

% 1. Remove unused genes in each organ model
for i=1:length(models)
    model = models{i,1};
    % "removeUnusedGenes" function comes from CobraTool
    [model, ~] = removeUnusedGenes(model); 
    models{i,1} = model;
end

% models_1 = models;

% 2. Add prefix for the metabolites and reactions of each organ model to 
% distinguish their contents.
modelID_rep = [Age, Sex, '-'];
for i=1:length(models)
    model = models{i,1};
    OrganName = [strrep(model.id, modelID_rep, ''), '_'];
    NewRxnID = cellfun(@(x) [OrganName x], model.rxns, 'UniformOutput', false);
    NewMetID = cellfun(@(x) [OrganName x], model.mets, 'UniformOutput', false);
    model.rxns = NewRxnID;
    model.mets = NewMetID;

    models{i,1} = model;
end

% 3. merge every organ models into whole-body model
n = 1;
while n<=numel(models)-1
    if n == 1
        MergeModel = mergeTwoModels(models{n,1},models{n+1,1});
    else
        MergeModel = mergeTwoModels(MergeModel,models{n+1,1});
    end
    n = n+1;
end

check_rxns = 0;
check_mets = 0;
for i=1:length(models)
    check_rxns = check_rxns+numel(models{i,1}.rxns);
    check_mets = check_mets+numel(models{i,1}.mets);
end

missing_rxns = {};
for i=1:length(models)
    model = models{i,1};
    Diff = setdiff(model.rxns, MergeModel.rxns);
    missing_rxns{end+1} = Diff;
end
    
    
% 4. Remove exchange reactions 
cMergeModel = MergeModel;
[selExc, ~] = findExcRxns(cMergeModel);
[cMergeModel, ~, ~] = removeRxns(cMergeModel, cMergeModel.rxns(selExc), 'metFlag', false);


% 5. Remove tasks that cannot be completed in human body.
removeTask = 'No_tasks.tsv'; % Some metabolic tasks should not be included in WBM.
removeData = readtable(['../Data/', removeTask], 'Delimiter', '\t', 'FileType', 'text');
removeRxns_id = removeData.WBMID;
removeRxns_id = intersect(removeRxns_id, cMergeModel.rxns);
[cMergeModel, ~, ~] = removeRxns(cMergeModel, removeRxns_id, 'metFlag', false);


% 6. Add biofluid metabolites and reactions.
% To establish a more physiologically meaningful association between organs,
% we referred to Thiele's work on the construction of WBMs, which intercon-
% nects all organs through biofluid compartments36. In total 13 biofluid c-
% ompartments were incorporated into the WBM reconstructions, including Di-
% et ([d]), Lumen ([lu]), Lumen of small intestine ([luSI]), Lumen of larg-
% e intestine ([luLI]), Feces ([fe]), Blood of circulation ([bc]), Blood 
% of portal vein ([bp]), Bile duct ([bd]), Cerebrospinal fluid ([csf]), Ur-
% ine ([u]), Sweat ([sw]), Breast milk (for females only, [mi]), and Air 
% ([a]). These compartments facilitate the exchange of metabolites with re-
% lated compartments and corresponding organs, thereby establishing connec-
% tivity among all organs in the WBMs. The established table of needed
% metabolites and reactions were loaded and added these contents into WBM.

 % set the bounds of rxns in WBM to 1000000.
 WBMLb = cMergeModel.lb;
 WBMLb(WBMLb == -1000) = -1000000;
 WBMLb(WBMLb == 1000) = 1000000;
 cMergeModel.lb = WBMLb;

 WBMUb = cMergeModel.ub;
 WBMUb(WBMUb == -1000) = -1000000;
 WBMUb(WBMUb == 1000) = 1000000;
 cMergeModel.ub = WBMUb;
 
 % Add compartments
 if Sex == 'Male'
     BiofluidComp = { 'luI'; 'bpI'; 'bp'; 'luSI'; 'bc'; 'luC'; 'bpC'; 'luLI';...
                      'lu'; 'fe'; 'bdL'; 'bpL'; 'bd'; 'luP'; 'bpP'; 'bcK';...
                      'u'; 'bpS'; 'csf'; 'a'; 'd'; 'aL'; 'swS'; 'sw'};
     BiofluidCompName = {'luSI-sIEC'; 'sIEC-bp'; 'Blood, portal vein'; 'Lumen, small intestine';...
                         'Blood, circulation'; 'luLI-Colon'; 'Colon-bp'; 'Lumen, large intestine';...
                         'Lumen'; 'Feces'; 'bd-Liver'; 'bp-Liver'; 'Bile duct'; 'lu-Pancreas';...
                         'Pancreas-bp'; 'Kidney-bc'; 'Urine'; 'Scord-bp'; 'Cerebrospinal fluid';...
                         'Air'; 'Diet'; 'air-Lung'; 'sw-Skin'; 'Sweat'};
     WBMComp = [cMergeModel.comps; BiofluidComp];
     WBMCompName = [cMergeModel.compNames; BiofluidCompName];
     cMergeModel.comps = WBMComp;
     cMergeModel.compNames = WBMCompName;
 elseif  Sex == 'Female'
     BiofluidComp = { 'luI'; 'bpI'; 'bp'; 'luSI'; 'bc'; 'luC'; 'bpC'; 'luLI';...
                      'lu'; 'fe'; 'bdL'; 'bpL'; 'bd'; 'luP'; 'bpP'; 'bcK';...
                      'u'; 'bpS'; 'csf'; 'a'; 'd'; 'aL'; 'swS'; 'sw'; 'miB'; 'mi'};
     BiofluidCompName = {'luSI-sIEC'; 'sIEC-bp'; 'Blood, portal vein'; 'Lumen, small intestine';...
                         'Blood, circulation'; 'luLI-Colon'; 'Colon-bp'; 'Lumen, large intestine';...
                         'Lumen'; 'Feces'; 'bd-Liver'; 'bp-Liver'; 'Bile duct'; 'lu-Pancreas';...
                         'Pancreas-bp'; 'Kidney-bc'; 'Urine'; 'Scord-bp'; 'Cerebrospinal fluid';...
                         'Air'; 'Diet'; 'air-Lung'; 'sw-Skin'; 'Sweat'; 'mi-Breast'; 'Breast milk'};
     WBMComp = [cMergeModel.comps; BiofluidComp];
     WBMCompName = [cMergeModel.compNames; BiofluidCompName];
     cMergeModel.comps = WBMComp;
     cMergeModel.compNames = WBMCompName;
 
 end
 
 % Add metabolites
 AddMetsData = readtable(['../Data/', Sex, '_Biofluid_mets.tsv'], 'Delimiter', '\t', 'FileType', 'text');
 metsToAdd = struct();
 metsToAdd.mets = AddMetsData.ID;
 metsToAdd.metNames = AddMetsData.Name;
 metsToAdd.compartments = AddMetsData.Compartment;
 metsToAdd.metFormulas = AddMetsData.Formula;
 metsToAdd.metCharges = AddMetsData.Charge;
 
 newModel=addMets(cMergeModel,metsToAdd); % A function in RAVEN.
 
 % Add reactions
 AddRxnsData = readtable(['../Data/', Sex, '_Biofluid_rxns.tsv'], 'Delimiter', '\t', 'FileType', 'text');
 rxnsToAdd = struct();
 rxnsToAdd.rxns = AddRxnsData.ID;
 rxnsToAdd.equations = AddRxnsData.Reaction;
 rxnsToAdd.mets = AddRxnsData.Metabolites;
 rxnsToAdd.stoichCoeffs = AddRxnsData.Coefficient;
 rxnsToAdd.lb = AddRxnsData.Lb;
 rxnsToAdd.ub = AddRxnsData.Ub;
 rxnsToAdd.subSystems = AddRxnsData.Subsystem;
 
 newModel1 = addRxns(newModel,rxnsToAdd); % A function in RAVEN.

% 7. set physiologically and stoichiometrically constrained modeling (PSCM)
% and Diet constraints into WBM. There has no physiological data for fetus.
% Thus, this step is no need to run for WBM of fetus.
if Age == 'Adult';
    sex =  lower(Sex);
    gender = lower(Sex);
    standardPhysiolDefaultParameters;
    IndividualParameters_ = IndividualParameters
    
    % set PSCM to WBM. This function comes from CobraTool.
    % There would get some warning about reactions not in the model, 
    % including 'Gall_H2Ot[bdG]','Brain_DM_atp_c_' and 'Heart_DM_atp_c_'. 
    % Gall is not included by this stduy. As for 'Brain_DM_atp_c_' and 
    % 'Heart_DM_atp_c_' is no need to demand in this research. 
    % Thus, we ignore them.
    WBM = physiologicalConstraintsHMDBbased(newModel1,IndividualParameters_);
    
    % set Diet constraints to WBM. This function comes from CobraTool
    EUAverageDietNew
    WBM = setDietConstraints(WBM, Diet); 
    
    % set some more constraints
    WBM = setSimulationConstraints(WBM);
    WBM.lb(strmatch('BBB_KYNATE[CSF]upt',WBM.rxns)) = -1000000;
    WBM.lb(strmatch('BBB_LKYNR[CSF]upt',WBM.rxns)) = -1000000;
    WBM.lb(strmatch('BBB_TRP_L[CSF]upt',WBM.rxns)) = -1000000;
    WBM.ub(strmatch('Brain_EX_glc_D(',WBM.rxns)) = -100; 
    
elseif Age == 'Elderly';
    sex =  lower(Sex);
    gender = lower(Sex);
    standardPhysiolDefaultParameters;
    % change some physiological parameter for elderly group.(PMID:15544194)
    if sex == 'male'
        IndividualParameters.bodyWight = 78.58;
        IndividualParameters.Height = 172.57;
        IndividualParameters.HeartRate = 65;
    elseif sex == 'female'
        IndividualParameters.bodyWight = 68.51;
        IndividualParameters.Height = 159.01;
        IndividualParameters.HeartRate = 68;
    else
        EM='Something wrong with the gender!';
        dispEM(EM);
    end
    
    IndividualParameters_ = IndividualParameters
    
    % set PSCM to WBM. This function comes from CobraTool.
    % There would get some warning about reactions not in the model, 
    % including 'Gall_H2Ot[bdG]','Brain_DM_atp_c_' and 'Heart_DM_atp_c_'. 
    % Gall is not included by this stduy. As for 'Brain_DM_atp_c_' and 
    % 'Heart_DM_atp_c_' is no need to demand in this research. 
    % Thus, we ignore them.
    WBM = physiologicalConstraintsHMDBbased(newModel1,IndividualParameters_);
    
    % set Diet constraints to WBM. This function comes from CobraTool
    EUAverageDietNew
    WBM = setDietConstraints(WBM, Diet); 
    
    % set some more constraints
    WBM = setSimulationConstraints(WBM);
    WBM.lb(strmatch('BBB_KYNATE[CSF]upt',WBM.rxns)) = -1000000;
    WBM.lb(strmatch('BBB_LKYNR[CSF]upt',WBM.rxns)) = -1000000;
    WBM.lb(strmatch('BBB_TRP_L[CSF]upt',WBM.rxns)) = -1000000;
    WBM.ub(strmatch('Brain_EX_glc_D(',WBM.rxns)) = -100; 
else
    EM="Something wrong with the Age group! Only allowed for 'Adult' or 'Elderly'."
    dispEM(EM);
    
% 8. Simulation for biomass maintain of WBM
newModel1 = changeObjective(newModel1, 'Whole_body_objective_rxn');
sol = optimizeCbModel(newModel1, 'max', 'one');

WBM = changeObjective(WBM, 'Whole_body_objective_rxn');
sol = optimizeCbModel(WBM, 'max', 'one');

WBM = releaseWBMConstraints(WBM, Sex, Age);

BMR_value = HumanBMR(WBM, sol);
save('Adult_Male_WBM.mat', 'WBM');


WBMc = releaseWBMConstraints(WBMc, Sex, Age);
WBMc = changeObjective(WBMc, 'Whole_body_objective_rxn');
[GeneClasses, RxnClasses, modelIrrevFM, MinimizedFlux] = pFBA(WBMc, 'geneoption',0, 'tol',1e-6, 'skipclass', 1);
BMR_value = HumanBMR(modelIrrevFM, MinimizedFlux);
    
