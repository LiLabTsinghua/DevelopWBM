%% This script is used to compare the synthesis capacity 
% of blood and urine metabolites in different age groups, serving as a
% marker of metabolites for aging.We obtained metabolites with 
% significantly and stably altered concentrations from the metaboAge database as 
% controls to evaluate the accuracy of WBMs. Afterwards, we simulated 
% alterations in all blood and urine metabolites. 

Sex = 'Male'; % or 'Female'
age = 'Adult'; % or 'Elderly'

%% Get results 
% load data
% DataPath = ['../data/', Sex, '_AgingMets.tsv'];
% AgeData = readtable(DataPath, 'Delimiter', '\t', 'FileType', 'text');
% Mets = AgeData.Met_id;
Mets_comp = 'bc'; % or 'u'

% age = Age{i,1};
Modelpath = ['../models/WBMs/', age, Sex, '_WBM_withPSCM.mat'];
load(Modelpath)
WBMc = WBM;
WBMc = releaseWBMConstraints(WBMc, Sex, age);
% [met_test, ~, ~] = intersect(Mets, WBMc.mets);
idx = findIndex(WBMc.comps, Mets_comp);
idx_all = find(WBMc.metComps == idx);
Mets = WBMc.mets(idx_all);
met_test = Mets;

PredValues = {};
for j=1:length(met_test)
    model = WBMc;
    if ismember(met_test{j,1}, model.mets)
        [model,rxnNames] = addDemandReaction(model, met_test(j,1));
        model = changeRxnBounds(model, rxnNames, 1000000, 'u');
        model = changeObjective(model, rxnNames);
        try
            % sol = optimizeCbModel(model, 'max', 'one');
            sol = solveLP(model);
            met_bc_value = sol.f;
        catch
            met_bc_value = [];
        end
%           PredValues{i,j} = met_bc_value;
    else
        met_bc_value = [];
    end
    PredValues{j,1} = met_bc_value;
    %disp(met_bc_value)
    disp(j)
end


Output = {};
Output(:,1) = Mets;
Output(:,2) = PredValues;
tbl = array2table(Output, 'VariableNames',{'Mets', 'DemandValue'});
outputFileName = ['../results/', age, Sex, '_', Mets_comp,'_MetsDemand.tsv'];
writetable(tbl, outputFileName, 'Delimiter', '\t', 'WriteRowNames', true, 'FileType', 'text');