%Run Growth rate prediction analysis on both models and ecModels for the
%analysed NCI60 cell-Lines subset subject to cummulative levels of exchange 
%fluxes constraints.
%
% Ivan Domenzain.      Last edited: 2019-12-03
% Jiahao Luo           Last edited: 2024-12-23

%% Model preparation
ihuman = readYAMLmodel('../models/Human2.yml');
taskStruct = '../data/metabolicTasks/metabolicTasks_Essential.txt';
prepData = prepHumanModelForftINIT(ihuman, false, taskStruct, '../data/reactions.tsv');

%% load NCI-60 data
load('../data/eGenesData/DepMap_RNAseq_data.mat')

cellNames   = {'HS_578T', 'ACH_000148';...
               'RPMI_8226', 'ACH_000817';...
               'HT29', 'ACH_000552';...
               'MALME_3M', 'ACH_000477';...
               'SR', 'ACH_000338';...
               'UO_31', 'ACH_000428';...
               'MDMAMB_231', 'ACH_000768';...
               'HOP62', 'ACH_000861';...
               'NCI_H226', 'ACH_000367';...
               'HOP92', 'ACH_000825';...
               'O_786', 'ACH_000649'};
cellidx = findIndex(depmap.cellID, cellNames(:,2));
% [~,cellidx] = ismember(depmap.cellID, cellNames(:,2));

data_struct.genes = depmap.genes;
data_struct.tissues = depmap.cellID(cellidx);
data_struct.levels = depmap.tpm(:,cellidx);
data_struct.threshold = 1;
data_struct

numSamp = length(data_struct.tissues);
models = cell(numSamp, 1);

%% generate context-specific models
mkdir('../models/NCI60')
for i = 1:numSamp
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
    % use '1+1' mode to reduce number of reactions without GPR
    m = ftINIT(prepData, data_struct.tissues{i}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), true, true);
    m.description = ['Automatically generated model for ', data_struct.tissues{i,1}];
    idx = findIndex(cellNames(:,2), data_struct.tissues{i,1});
    cellNames_ = cellNames(:,1);
    m.id = cellNames_{idx,1};
    mkdir(['../models/NCI60/',m.id]);
    model = m;
    save(['../models/NCI60/',m.id,'/model.mat'], 'model')
    models{i} = m;
    % disp(data_struct.tissues{i,1});
end

%% Generate ec-context-specific models
% More details can be seen in the tutorials of GECKO3(https://github.com/SysBioChalmers/GECKO/tree/main/tutorials/light_ecModel)
%% Simulation

constraints = [{'media'} {'glucose'} {'L-lactate'} {'threonine'}];

pred_Grates    = [];
pred_EC_Grates = [];
meanError      = [];
mean_EC_error  = [];
for constLevel = [0 1 2 3]
    %Run growth simulations for ecModels
    [meanErr,pred,meas_GRates] = ExchFluxesComparison_NCI60(constLevel,true,false);
    pred_EC_Grates   = [pred_EC_Grates, pred];
    mean_EC_error    = [mean_EC_error; meanErr]; 
    %Run growth simulations for standard models
    [meanErr,pred,~] = ExchFluxesComparison_NCI60(constLevel,false,false);
    pred_Grates      = [pred_Grates, pred];
    meanError        = [meanError; meanErr]; 
end
experimental = meas_GRates';
cellNames    = cellNames(:,1);
%Write results ecModels
results    = table(cellNames,experimental);
errorTable = table(cellNames);
for i=1:length(constraints)
    results    = [results table(pred_EC_Grates(:,i),'VariableNames',constraints(i))];
    errorVec   = abs(experimental-pred_EC_Grates(:,i))./experimental;
    errorTable = [errorTable table(errorVec,'VariableNames',constraints(i))];
end


% save results
folderpath = '../results/11_cellLines_NCI60/';
if ~exist(folderpath, 'dir')
    mkdir(folderpath)
end

fileName = '../results/11_cellLines_NCI60/ecModels_gRates.txt';
writetable(results,fileName,'Delimiter','\t','QuoteStrings',false);
fileName = '../results/11_cellLines_NCI60/ecModels_error_gRates.txt';
writetable(errorTable,fileName,'Delimiter','\t','QuoteStrings',false);
%Write results Models
results  = table(cellNames,experimental);
errorTable = table(cellNames);
for i=1:length(constraints)
    results    = [results table(pred_Grates(:,i),'VariableNames',constraints(i))];
    errorVec   = abs(experimental-pred_Grates(:,i))./experimental;
    errorTable = [errorTable table(errorVec,'VariableNames',constraints(i))];
end
fileName = '../results/11_cellLines_NCI60/models_gRates.txt';
writetable(results,fileName,'Delimiter','\t','QuoteStrings',false);
fileName = '../results/11_cellLines_NCI60/models_error_gRates.txt';
writetable(errorTable,fileName,'Delimiter','\t','QuoteStrings',false);
