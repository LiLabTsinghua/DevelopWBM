% This script is used to merge organ models into whole-body model.
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


% check_rxns = 0;
% check_mets = 0;
% for i=1:length(models)
%     check_rxns = check_rxns+numel(models{i,1}.rxns);
%     check_mets = check_mets+numel(models{i,1}.mets);
% end

% missing_rxns = {};
% for i=1:length(models)
%     model = models{i,1};
%     Diff = setdiff(model.rxns, MergeModel.rxns);
%     missing_rxns{end+1} = Diff;
% end
    

% 4. Remove exchange reactions 
cMergeModel = MergeModel;
[selExc, ~] = findExcRxns(cMergeModel);
[cMergeModel, ~, ~] = removeRxns(cMergeModel, cMergeModel.rxns(selExc), 'metFlag', false);


% 5. Remove tasks that cannot be completed in human body.
cMergeModel = WBMTaskRxns(cMergeModel, Sex);


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
 
 % Add biofluid compartments
 if strcmp(Sex, 'Male')
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
 elseif  strcmp(Sex, 'Female')
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
 else
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
 end
 
 % Add metabolites
 if strcmp(Sex, 'Male')
     AddMetsData = readtable('../data/metabolicTasks/Biofluids/Harvey_biofluid_mets.tsv', 'Delimiter', '\t', 'FileType', 'text');
 elseif strcmp(Sex, 'Female')
     AddMetsData = readtable('../data/metabolicTasks/Biofluids/Harvetta_biofluid_mets.tsv', 'Delimiter', '\t', 'FileType', 'text');
 else
     AddMetsData = readtable('../data/metabolicTasks/Biofluids/Fetus_biofluid_mets.tsv', 'Delimiter', '\t', 'FileType', 'text');
 end

 metsToAdd = struct();
 metsToAdd.mets = AddMetsData.Met_id;
 metsToAdd.metNames = AddMetsData.Met_name;
 metsToAdd.compartments = AddMetsData.Met_compartment;
 metsToAdd.metFormulas = AddMetsData.Met_formula;
 metsToAdd.metCharges = AddMetsData.Met_charge;
 
 newModel=addMets(cMergeModel,metsToAdd); % A function in RAVEN.
 
 % Add reactions
 if strcmp(Sex, 'Male')
     AddRxnsData = readtable('../data/metabolicTasks/Biofluids/Harvey_biofluid_rxns.tsv', 'Delimiter', '\t', 'FileType', 'text');
 elseif strcmp(Sex, 'Female')
     AddRxnsData = readtable('../data/metabolicTasks/Biofluids/Harvetta_biofluid_rxns.tsv', 'Delimiter', '\t', 'FileType', 'text');
 else
     AddRxnsData = readtable('../data/metabolicTasks/Biofluids/Fetus_biofluid_rxns.tsv', 'Delimiter', '\t', 'FileType', 'text');
 end


 % AddRxnsData = readtable(['../Data/', Sex, '_Biofluid_rxns.tsv'], 'Delimiter', '\t', 'FileType', 'text');
 % AddRxnsData = readtable(['../Data/', 'Male_Biofluid_rxns_0910.tsv'], 'Delimiter', '\t', 'FileType', 'text');
 % AddRxnsData = readtable(['../Data/', 'Female_Biofluid_rxns_0924_end.tsv'], 'Delimiter', '\t', 'FileType', 'text');
 
 rxnsToAdd = struct();
 rxnsToAdd.rxns = AddRxnsData.ID;
 rxnsToAdd.equations = AddRxnsData.Reaction;
 % rxnsToAdd.mets = AddRxnsData.Metabolites;
 rxnsToAdd.mets = AddRxnsData.Metabolites_new;
 % rxnsToAdd.stoichCoeffs = AddRxnsData.Coefficient;
 rxnsToAdd.stoichCoeffs = AddRxnsData.Mets_Coefficient_new;
 rxnsToAdd.lb = AddRxnsData.Lb;
 rxnsToAdd.ub = AddRxnsData.Ub;
 rxnsToAdd.subSystems = AddRxnsData.Subsystem;
 rxnsToAdd.grRules = AddRxnsData.GPR;
 
  parfor i=1:length(rxnsToAdd.mets)
     mets = eval(rxnsToAdd.mets{i,1}).';
     if all(ismember(mets, newModel.mets))
         missingMet(i,1) = 1;
     else
         missingMet(i,1) = 0;
     end
 end
 
 index_of_zeros = find(missingMet == 0);
 let_go_rxns = rxnsToAdd.rxns(index_of_zeros);
 
 [TF, loc] = ismember(newModel.rxns, rxnsToAdd.rxns);
 % commonrxnsIdx = find(TF);
 commonrxnsIdx = findIndex(rxnsToAdd.rxns, newModel.rxns(find(TF)));
 
 index_of_zeros = [index_of_zeros; commonrxnsIdx];
 
 cellLength = length(rxnsToAdd.rxns);
 keepRows = true(cellLength, 1);
 keepRows(index_of_zeros) = false;
 fieldNames = fieldnames(rxnsToAdd);
for i = 1:length(fieldNames)
    fieldName = fieldNames{i};
    if iscell(rxnsToAdd.(fieldName))
        rxnsToAdd.(fieldName) = rxnsToAdd.(fieldName)(keepRows);
    elseif isnumeric(rxnsToAdd.(fieldName)) && ismatrix(rxnsToAdd.(fieldName)
        if isvector(rxnsToAdd.(fieldName))
            rxnsToAdd.(fieldName) = rxnsToAdd.(fieldName)(:);
            rxnsToAdd.(fieldName) = rxnsToAdd.(fieldName)(keepRows);
        else
            rxnsToAdd.(fieldName) = rxnsToAdd.(fieldName)(keepRows, :);
        end
    end
end
 
 
 eqnType=1;
 compartment = [];
 allowNewMets=false;
 allowNewGenes=true;
 newModel1=addRxns(newModel,rxnsToAdd,eqnType,compartment,allowNewMets,allowNewGenes);
 
 
% 7. set physiologically and stoichiometrically constrained modeling (PSCM)
% and Diet constraints into WBM. There has no physiological data for fetus.
% Thus, this step is no need to run for WBM of fetus.
if strcmp(Age, 'Adult')
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
    
elseif strcmp(Age, 'Elderly')
    sex =  lower(Sex);
    gender = lower(Sex);
    standardPhysiolDefaultParameters;
    % change some physiological parameter for elderly group.(PMID:15544194)
    if strcmp(sex, 'male')
        % IndividualParameters.bodyWeight = 78.58;
        IndividualParameters.bodyWeight = 65;
        % IndividualParameters.bodyWeight = 70;
        IndividualParameters.Height = 172.57;
        IndividualParameters.HeartRate = 65;
        
        IndividualParameters.StrokeVolume = 70;
        IndividualParameters.Hematocrit = 0.35;
        
    elseif strcmp(sex, 'female')
        % IndividualParameters.bodyWeight = 68.51;
        IndividualParameters.bodyWeight = 54;
        IndividualParameters.Height = 159.01;
        IndividualParameters.HeartRate = 65;
        
        IndividualParameters.StrokeVolume = 70;
        IndividualParameters.Hematocrit = 0.35;
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
    EM="Something wrong with the Age group! Only allowed for 'Adult' or 'Elderly'.";
    dispEM(EM);
end
    
% 8. Simulation for biomass maintain of WBM
newModel1 = changeObjective(newModel1, 'Whole_body_objective_rxn');
sol = optimizeCbModel(newModel1, 'max', 'one');
save(['../models/WBMs/',Age,Sex,'_WBM_withoutPSCM.mat'], 'newModel1')

WBM = releaseWBMConstraints(WBM, Sex, Age);
WBM = changeObjective(WBM, 'Whole_body_objective_rxn');
%sol = optimizeCbModel(WBM, 'max', 'one');
sol = solveLP(WBM);
sol.f
save(['../models/WBMs/',Age,Sex,'_WBM_withPSCM.mat'], 'WBM')
