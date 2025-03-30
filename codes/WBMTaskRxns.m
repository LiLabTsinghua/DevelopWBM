function newModel = WBMTaskRxns(model, Gender)

if nargin < 2 
    Gender = 'Male';
end

newModel = model;

% Remove tasks that cannot present in WBM.
removeTask = 'No_tasks.tsv'; % Some metabolic tasks should not be included in WBM.
removeData = readtable(['../data/metabolicTasks/', removeTask], 'Delimiter', '\t', 'FileType', 'text');
removeRxns_id = removeData.WBMID;
removeRxns_id = intersect(removeRxns_id, newModel.rxns);
[newModel, ~, ~] = removeRxns(newModel, removeRxns_id, 'metFlag', false);

% Add tasks that should present in WBM.
if strcmp(Gender, 'Male')
    AddMetsData = readtable('../data/metabolicTasks/Male_PresentTasks_Mets.tsv', 'Delimiter', '\t', 'FileType', 'text');
    AddRxnsData = readtable('../data/metabolicTasks/Male_PresentTasks_Rxns.tsv', 'Delimiter', '\t', 'FileType', 'text');
else
    AddMetsData = readtable('../data/metabolicTasks/Female_PresentTasks_Mets.tsv', 'Delimiter', '\t', 'FileType', 'text');
    AddRxnsData = readtable('../data/metabolicTasks/Female_PresentTasks_Rxns.tsv', 'Delimiter', '\t', 'FileType', 'text');
end

% add tasks mets
AddMets_id = AddMetsData.Met_id;
[tf, ~] = ismember(AddMets_id, newModel.mets);
idx = find(~tf);

metsToAdd = struct();
metsToAdd.mets = AddMetsData.Met_id(idx);
metsToAdd.metNames = AddMetsData.Met_name(idx);
metsToAdd.compartments = AddMetsData.Met_compartment(idx);
metsToAdd.metFormulas = AddMetsData.Met_formula(idx);
metsToAdd.metCharges = AddMetsData.Met_charge(idx);

newModel=addMets(newModel,metsToAdd);

% add tasks rxns
AddRxns_id = AddRxnsData.Rxn_id;
[tf, ~] = ismember(AddRxns_id, newModel.rxns);
idx = find(~tf);

rxnsToAdd = struct();
rxnsToAdd.rxns = AddRxnsData.Rxn_id(idx);
rxnsToAdd.equations = AddRxnsData.Rxn_eq(idx);
rxnsToAdd.mets = AddRxnsData.Mets_id(idx);
rxnsToAdd.stoichCoeffs = AddRxnsData.Mets_coff(idx);
rxnsToAdd.lb = AddRxnsData.Lb(idx);
rxnsToAdd.ub = AddRxnsData.Ub(idx);
rxnsToAdd.subSystems = AddRxnsData.Rxn_sub(idx);
rxnsToAdd.grRules = AddRxnsData.Rxn_gpr(idx);

eqnType=1;
compartment = [];
allowNewMets=false;
allowNewGenes=true;
newModel=addRxns(newModel,rxnsToAdd,eqnType,compartment,allowNewMets,allowNewGenes);

end