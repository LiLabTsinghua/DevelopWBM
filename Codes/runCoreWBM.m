% This script is used to generate core WBM which included Brain, Heart, 
% Liver, Lung, Kidney, Muscle, Skin.
% The supported options are 'male', 'female', and 'fetus'. Please
% define those using the variable 'sex' (e.g., sex = 'male').

Age = 'Adult';
Sex = 'Male';
modelPath = ['../models/', Age, Sex, '_models.mat'];
load(modelPath);
% models = load(modelPath);

% 1. Remove unused genes in each organ model
for i=1:length(models)
    model = models{i,1};
    % "removeUnusedGenes" function from CobraTool
    [model, ~] = removeUnusedGenes(model); 
    models{i,1} = model;
end

% 2. Select core organs
CoreOrgans = {'Brain'; 'Heart'; 'Liver'; 'Lung'; 'Muscle'; 'Kidney'; 'Skin'};
CoreModels = {};
for i=1:length(models)
    model = models{i,1};
    namesA = strrep(model.id, [Age Sex '-'], '');
    if ismember(namesA, CoreOrgans)
        CoreModels{end+1} = model;
    end
end
CoreModels = CoreModels.';

% 3. Add prefix for the metabolites and reactions of each organ model to 
% distinguish their contents.
modelID_rep = [Age, Sex, '-'];
for i=1:length(CoreModels)
    model = CoreModels{i,1};
    OrganName = [strrep(model.id, modelID_rep, ''), '_'];
    NewRxnID = cellfun(@(x) [OrganName x], model.rxns, 'UniformOutput', false);
    NewMetID = cellfun(@(x) [OrganName x], model.mets, 'UniformOutput', false);
    NewMetName = cellfun(@(x) [OrganName x], model.metNames, 'UniformOutput', false);
    model.rxns = NewRxnID;
    model.mets = NewMetID;
    model.metNames = NewMetName;

    CoreModels{i,1} = model;
end

% 4. merge every organ models into core-WBM
n = 1;
while n<=numel(CoreModels)-1
    if n == 1
        MergeModel = mergeTwoModels(CoreModels{n,1},CoreModels{n+1,1});
    else
        MergeModel = mergeTwoModels(MergeModel,CoreModels{n+1,1});
    end
    n = n+1;
end


% 5. Remove exchange reactions 
cMergeModel = MergeModel;
[selExc, ~] = findExcRxns(cMergeModel);
[cMergeModel, ~, ~] = removeRxns(cMergeModel, cMergeModel.rxns(selExc), 'metFlag', false);

% update comp of mets
new_comp = cell(size(cMergeModel.mets));
for i = 1:numel(cMergeModel.mets)
    currentString = cMergeModel.mets{i};
    lastChar = currentString(end);
    new_comp{i} = lastChar;
end

new_comp_num = [];
for i=1:length(new_comp)
    new_comp_char = new_comp{i,1};
    for j=1:length(cMergeModel.comps)
        model_comp = cMergeModel.comps{j,1};
        if strcmp(new_comp_char, model_comp)
            new_comp_num(end+1) = j;
        end
    end
end

cMergeModel.metComps = new_comp_num.';
        



% 6. Remove tasks that cannot be completed in human body.
removeTask = 'No_tasks.tsv'; % Some metabolic tasks should not be included in WBM.
removeData = readtable(['../data/metabolicTasks/', removeTask], 'Delimiter', '\t', 'FileType', 'text');
removeRxns_id = removeData.WBMID;
removeRxns_id = intersect(removeRxns_id, cMergeModel.rxns);
[cMergeModel, ~, ~] = removeRxns(cMergeModel, removeRxns_id, 'metFlag', false);

if strcmp(Sex, 'Male')
    AddMetsData = readtable('../data/metabolicTasks/Male_PresentTasks_Mets.tsv', 'Delimiter', '\t', 'FileType', 'text');
    AddRxnsData = readtable('../data/metabolicTasks/Male_PresentTasks_Rxns.tsv', 'Delimiter', '\t', 'FileType', 'text');
elseif strcmp(Sex, 'Female')
    AddMetsData = readtable('../data/metabolicTasks/Female_PresentTasks_Mets.tsv', 'Delimiter', '\t', 'FileType', 'text');
    AddRxnsData = readtable('../data/metabolicTasks/Female_PresentTasks_Rxns.tsv', 'Delimiter', '\t', 'FileType', 'text');
end

Met_index = [];
All_need_mets = AddMetsData.Met_id;
All_need_met_names = AddMetsData.Met_name;
for i=1:length(All_need_mets)
    met_id = All_need_mets{i,1};
    met_name = All_need_met_names{i,1};
    met_logical = 0;
    for j=1:length(CoreOrgans)
        organ = CoreOrgans{j,1};
        if startsWith(met_id, organ)
            if ~ismember(met_id, cMergeModel.mets)
%                 if ~ismember(met_name, cMergeModel.metNames)
                met_logical = 1;
%                 end
            end
        end
    end
    Met_index(end+1) = met_logical;
end
Met_index = find(Met_index.');


metsToAdd_task = struct();
metsToAdd_task.mets = AddMetsData.Met_id(Met_index);
metsToAdd_task.metNames = AddMetsData.Met_name(Met_index);
metsToAdd_task.compartments = AddMetsData.Met_compartment(Met_index);
metsToAdd_task.metFormulas = AddMetsData.Met_formula(Met_index);
metsToAdd_task.metCharges = AddMetsData.Met_charge(Met_index);

cMergeModel=addMets(cMergeModel,metsToAdd_task); % A function in RAVEN.



Rxn_index = [];
All_need_rxns = AddRxnsData.Rxn_id;
All_need_Rxn_names = AddRxnsData.Rxn_name;
for i=1:length(All_need_rxns)
    rxn_id = All_need_rxns{i,1};
    rxn_name = All_need_Rxn_names{i,1};
    rxn_logical = 0;
    for j=1:length(CoreOrgans)
        organ = CoreOrgans{j,1};
        if startsWith(rxn_id, organ)
            if ~ismember(rxn_id, cMergeModel.rxns)
%                 if ~ismember(met_name, cMergeModel.metNames)
                rxn_logical = 1;
%                 end
            end
        end
    end
    Rxn_index(end+1) = rxn_logical;
end
Rxn_index = find(Rxn_index.');

rxnsToAdd_task = struct();
rxnsToAdd_task.rxns = AddRxnsData.Rxn_id(Rxn_index);
rxnsToAdd_task.equations = AddRxnsData.Rxn_eq(Rxn_index);
rxnsToAdd_task.mets = AddRxnsData.Mets_id(Rxn_index);
rxnsToAdd_task.stoichCoeffs = AddRxnsData.Mets_coff(Rxn_index);
rxnsToAdd_task.lb = AddRxnsData.Lb(Rxn_index);
rxnsToAdd_task.ub = AddRxnsData.Ub(Rxn_index);
rxnsToAdd_task.subSystems = AddRxnsData.Rxn_sub(Rxn_index);

cMergeModel = addRxns(cMergeModel,rxnsToAdd_task); % A function in RAVEN.

cMergeModel_1 = cMergeModel;
cMergeModel = cMergeModel_1;


% 7. Add metabolites for for transfer in various organs.
AddMetsData = readtable(['../data/metabolicTasks/', 'CoreTrans_Mets.tsv'], 'Delimiter', '\t', 'FileType', 'text');
[isCommon, idx] = ismember(AddMetsData.ID, cMergeModel.mets);

metsToAdd = struct();
metsToAdd.mets = AddMetsData.ID(~isCommon);
metsToAdd.metNames = AddMetsData.Name(~isCommon);
metsToAdd.compartments = AddMetsData.Compartment(~isCommon);
metsToAdd.metFormulas = AddMetsData.Formula(~isCommon);
metsToAdd.metCharges = AddMetsData.Charge(~isCommon);

newModel=addMets(cMergeModel,metsToAdd); % A function in RAVEN.

% 8. Add reactions for for transfer in various organs.
AddRxnsData = readtable(['../data/metabolicTasks/', 'CoreTrans_Rxns.tsv'], 'Delimiter', '\t', 'FileType', 'text');
[isCommon, idx] = ismember(AddRxnsData.ID, newModel.rxns);

rxnsToAdd = struct();
rxnsToAdd.rxns = AddRxnsData.ID(~isCommon);
rxnsToAdd.equations = AddRxnsData.Reaction(~isCommon);
rxnsToAdd.mets = AddRxnsData.Metabolites(~isCommon);
rxnsToAdd.stoichCoeffs = AddRxnsData.Coefficient(~isCommon);
rxnsToAdd.lb = AddRxnsData.Lb(~isCommon);
rxnsToAdd.ub = AddRxnsData.Ub(~isCommon);
rxnsToAdd.subSystems = AddRxnsData.Subsystem(~isCommon);

newModel1 = addRxns(newModel,rxnsToAdd); % A function in RAVEN.

save('../models/AdultMale_coreWBM_0324.mat', 'newModel1')

% 9. 
% The uptake reactions were classified.
CarbohydrateRxns = {'Diet_EX_xyl_D[d]'; 'Diet_EX_gal[d]'; 'Diet_EX_xylt[d]';...
                    'Diet_EX_rib_D[d]'; 'Diet_EX_strch1[d]'; 'Diet_EX_malt[d]';...
                    'Diet_EX_sucr[d]'; 'Diet_EX_lcts[d]'; 'Diet_EX_fru[d]';...
                    'Diet_EX_sbt_D[d]'; 'Diet_EX_glc_D[d]'; 'Diet_EX_acgam[d]'};
                
FattyRxns = {'Diet_EX_ocdca[d]'; 'Diet_EX_lnlnca[d]'; 'Diet_EX_lnlc[d]';...
             'Diet_EX_lgnc[d]'; 'Diet_EX_lanost[d]'; 'Diet_EX_ddca[d]';...
             'Diet_EX_octa[d]'; 'Diet_EX_tdchola[d]'; 'Diet_EX_CE2510[d]';...
             'Diet_EX_arachd[d]'; 'Diet_EX_vitd3[d]'; 'Diet_EX_avite1[d]';...
             'Diet_EX_ttdca[d]'; 'Diet_EX_clpnd[d]'; 'Diet_EX_caro[d]';...
             'Diet_EX_dgchol[d]'; 'Diet_EX_2obut[d]'; 'Diet_EX_hdca[d]';...
             'Diet_EX_hpdca[d]'; 'Diet_EX_retinol[d]'; 'Diet_EX_phyQ[d]';...
             'Diet_EX_dlnlcg[d]'; 'Diet_EX_arach[d]'; 'Diet_EX_dca[d]';...
             'Diet_EX_chsterol[d]'; 'Diet_EX_ptrc[d]'; 'Diet_EX_docosac[d]';...
             'Diet_EX_but[d]'; 'Diet_EX_crvnc[d]'; 'Diet_EX_doco13ac[d]';...
             'Diet_EX_M01989[d]'; 'Diet_EX_ptdca[d]'; 'Diet_EX_tchola[d]';...
             'Diet_EX_hdcea[d]'; 'Diet_EX_ac[d]'};
         
ProteinRxns = {'Diet_EX_cys_L[d]'; 'Diet_EX_glu_L[d]'; 'Diet_EX_thr_L[d]';...
               'Diet_EX_arg_L[d]'; 'Diet_EX_asp_L[d]'; 'Diet_EX_pro_L[d]';...
               'Diet_EX_ile_L[d]'; 'Diet_EX_trp_L[d]'; 'Diet_EX_his_L[d]';...
               'Diet_EX_lys_L[d]'; 'Diet_EX_phe_L[d]'; 'Diet_EX_ala_L[d]';...
               'Diet_EX_ser_L[d]'; 'Diet_EX_gly[d]'; 'Diet_EX_Lcystin[d]';...
               'Diet_EX_leu_L[d]'; 'Diet_EX_met_L[d]'; 'Diet_EX_asn_L[d]';...
               'Diet_EX_gln_L[d]'; 'Diet_EX_ala_D[d]'; 'Diet_EX_tyr_L[d]';...
               'Diet_EX_val_L[d]'};
           
NucleicAcidRxns = {'Diet_EX_thm[d]'; 'Diet_EX_gsn[d]'; 'Diet_EX_gua[d]';...
                   'Diet_EX_urate[d]'; 'Diet_EX_csn[d]'; 'Diet_EX_cytd[d]';...
                   'Diet_EX_amet[d]'; 'Diet_EX_dad_2[d]'; 'Diet_EX_hxan[d]';...
                   'Diet_EX_adn[d]'; 'Diet_EX_amp[d]'; 'Diet_EX_ade[d]';...
                   'Diet_EX_thymd[d]'; 'Diet_EX_dgsn[d]'; 'Diet_EX_ura[d]';...
                   'Diet_EX_xan[d]'; 'Diet_EX_dcyt[d]'};

% close all uptake reactions of nutrient metabolites
coreWBM = newModel1;
CarbohydrateRxnsIdx = findIndex(coreWBM.rxns, CarbohydrateRxns);
FattyRxnsIdx = findIndex(coreWBM.rxns, FattyRxns);
ProteinRxnsIdx = findIndex(coreWBM.rxns, ProteinRxns);
NucleicAcidRxnsIdx = findIndex(coreWBM.rxns, NucleicAcidRxns);

NutritionRxnsIdx = [CarbohydrateRxnsIdx; FattyRxnsIdx;...
                    ProteinRxnsIdx; NucleicAcidRxnsIdx];

coreWBM.lb(NutritionRxnsIdx) = 0;

% Add reaction of WBM-ATPM and set it as objective
ATPM = 'MAR03964';
OrganATPMidx = cellfun(@(x) contains(x, 'MAR03964'), coreWBM.rxns);
OrganATPM = coreWBM.rxns(OrganATPMidx);

coreWBM_1 = coreWBM;
coreWBM = coreWBM_1;

Pseudo_Mets = {};
for i=1:length(OrganATPM)
    pseudoMet = strrep(OrganATPM{i,1}, ATPM, 'pseudo_ATP');
    metsToAdd = struct();
    metsToAdd.mets = pseudoMet;
    metsToAdd.metNames = pseudoMet;
    metsToAdd.compartments = 'e';
    metsToAdd.metFormulas = 'X';
    metsToAdd.metCharges = 0;
    coreWBM=addMets(coreWBM,metsToAdd);
    
    % add pseudo met into the reaction of organ ATPM as a product
    m = findIndex(coreWBM.mets, pseudoMet);
    r = findIndex(coreWBM.rxns, OrganATPM{i,1});
    coreWBM.S(m, r) = 1;
    
    Pseudo_Mets{end+1} = pseudoMet;
end

% add WBM-ATPM reaction as objective
% set stoichiometric coefficient of objective rxn(Bionumbers: 109725).
% In order to prevent the simulation result from having a flux value excee-
% ding 1000, the original BMR ratio data of each organ was divided by 10.
% The lungs and skin split the remainder of the data equally according to 
% weight.
CoeffPseudoMets = {'Brain', 1.9;
                  'Heart', 0.7;
                  'Liver', 2.7;
                  'Lung', 0.325;
                  'Kidney', 1;
                  'Muscle', 1.8;
                  'Skin', 1.575};
              
Pseudo_Mets = Pseudo_Mets.';
for i=1:length(Pseudo_Mets)
    for j=1:length(CoeffPseudoMets)
        if contains(Pseudo_Mets{i,1}, CoeffPseudoMets{j,1})
            Pseudo_Mets{i,2} = num2str(-CoeffPseudoMets{j,2});
        end
    end
end


rxnsToAdd = struct();
rxnsToAdd.rxns = 'WBM-ATPM';
% rxnsToAdd.equations = '=> ';
rxnsToAdd.mets = Pseudo_Mets(:,1);
rxnsToAdd.stoichCoeffs = cell2mat(cellfun(@str2double, Pseudo_Mets(:, 2), 'UniformOutput', false));
rxnsToAdd.lb = 0;
rxnsToAdd.ub = 1000;
% rxnsToAdd.subSystems = AddRxnsData.Subsystem;
coreWBM = addRxns(coreWBM,rxnsToAdd);

% check the rxn
printRxnFormula(coreWBM, 'WBM-ATPM');
printRxnFormula(coreWBM, 'Heart_MAR03964')

% check simulation
coreWBM.lb(findIndex(coreWBM.rxns, 'EX_Liver_MAM03161e')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'EX_Liver_MAM03161e')) = 0;
coreWBM.lb(findIndex(coreWBM.rxns, 'EX_MAM02674e')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'EX_MAM02674e')) = 0;

coreWBM.lb(findIndex(coreWBM.rxns, 'Diet_EX_etoh[d]')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'Diet_EX_etoh[d]')) = 0;
coreWBM.lb(findIndex(coreWBM.rxns, 'Diet_EX_for[d]')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'Diet_EX_for[d]')) = 0;

coreWBM.lb(findIndex(coreWBM.rxns, CarbohydrateRxns)) = -10;
coreWBM.lb(findIndex(coreWBM.rxns, FattyRxns)) = -10;
coreWBM.lb(findIndex(coreWBM.rxns, ProteinRxns)) = -10;

coreWBM = changeObjective(coreWBM, 'WBM-ATPM');
sol = optimizeCbModel(coreWBM, 'max', 'one');

save('../models/AdultMale_coreWBM.mat', 'coreWBM')

%% Simulation for food uptake

% You can load the model directly without running the above code.
Age = 'Adult';
Sex = 'Male';
modelPath = ['../models/', Age, Sex, '_coreWBM.mat'];
load(modelPath);

FoodData = readtable(['../data/Food_data/', 'FNDDS_Food_Data.tsv'], 'Delimiter', '\t', 'FileType', 'text');
ProteinData = FoodData.Protein;
LipidData = FoodData.Lipid;
CarbohydrateData = FoodData.Carbohydrate;
EnergyData = FoodData.Energy;


% Close this exchange reactions which were used for fasting simulation
coreWBM.lb(findIndex(coreWBM.rxns, 'Brain_Exchange_glycogen')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'Brain_Exchange_glycogen')) = 0;
coreWBM.lb(findIndex(coreWBM.rxns, 'Heart_Exchange_glycogen')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'Heart_Exchange_glycogen')) = 0;
coreWBM.lb(findIndex(coreWBM.rxns, 'Kidney_Exchange_glycogen')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'Kidney_Exchange_glycogen')) = 0;
coreWBM.lb(findIndex(coreWBM.rxns, 'Lung_Exchange_glycogen')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'Lung_Exchange_glycogen')) = 0;
coreWBM.lb(findIndex(coreWBM.rxns, 'Muscle_Exchange_glycogen')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'Muscle_Exchange_glycogen')) = 0;
coreWBM.lb(findIndex(coreWBM.rxns, 'Skin_Exchange_glycogen')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'Skin_Exchange_glycogen')) = 0;

coreWBM.lb(findIndex(coreWBM.rxns, 'Liver_Trans_MAM01307e')) = -1000;
coreWBM.ub(findIndex(coreWBM.rxns, 'Liver_Trans_MAM01307e')) = 1000;
coreWBM.lb(findIndex(coreWBM.rxns, 'Muscle_Trans_MAM01307e')) = -1000;
coreWBM.ub(findIndex(coreWBM.rxns, 'Muscle_Trans_MAM01307e')) = 1000;

coreWBM.lb(findIndex(coreWBM.rxns, 'Liver_Exchange_Lac')) = -1000;
coreWBM.ub(findIndex(coreWBM.rxns, 'Liver_Exchange_Lac')) = 0;
coreWBM.lb(findIndex(coreWBM.rxns, 'Muscle_Exchange_Lac')) = 0;
coreWBM.ub(findIndex(coreWBM.rxns, 'Muscle_Exchange_Lac')) = 1000;


PrdictEnergy = zeros(length(EnergyData), 1);
% It will take a while.
parfor i=1:length(ProteinData)
    coreWBMc = coreWBM;
    
    cab = -CarbohydrateData(i,1);
    lip = -LipidData(i,1);
    pro = -ProteinData(i,1);
    
    % 
    coreWBMc.lb(findIndex(coreWBMc.rxns, CarbohydrateRxns)) = cab*0.2;
    coreWBMc.lb(findIndex(coreWBMc.rxns, FattyRxns)) = lip*0.2;
    coreWBMc.lb(findIndex(coreWBMc.rxns, ProteinRxns)) = pro*0.2;

    coreWBMc = changeObjective(coreWBMc, 'WBM-ATPM');
    %sol = optimizeCbModel(coreWBMc, 'max', 'one');
    sol = solveLP(coreWBMc, '1');
    
    if sol.stat == 1
        PrdictEnergy(i,1) = sol.f;
    end
end

%% output results
FoodData = addvars(FoodData, PrdictEnergy);
FoodData.Properties.VariableNames{end} = 'PrdictEnergy';
outputFileName = '../results/Result_FoodEnergy_Prediction.tsv';
writetable(FoodData, outputFileName, 'Delimiter', '\t', 'FileType', 'text');
