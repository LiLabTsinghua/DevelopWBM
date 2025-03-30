function [models] = Gen_Organ_models(model, Age, Sex)

if nargin < 3
    Sex = 'Male';
end

if nargin < 3
    Age = 'Adult';
end

if strcmp(Age, 'Fetus')
    Sex = '';
end

% Generate organ-specific models
taskStruct = '../data/metabolicTasks/metabolicTasks_Essential.txt';
prepData = prepHumanModelForftINIT(model, false, taskStruct, '../data/reactions.tsv');

% load data
load('../models/RBC.mat');
% modelPath = ['../Models/models', '_', Age, Sex];
gtex_data = readtable(['../data/Organ_Exp/', Age, Sex, '_expression_data.txt']);
[~, n] = size(gtex_data);
numSamp = n-1; %the first two columns are the genes in ENSEMBL

gtex_data(1:5, 1:5)

data_struct.genes = gtex_data{:, 1}; % gene names
data_struct.tissues = gtex_data.Properties.VariableNames(2:n); % sample (tissue) names
data_struct.levels = gtex_data{:, 2:n}; % gene TPM values

% set 25th percentile of metabolic gene expression values as the threshold
if strcmp(Age, 'Fetus')
    % use the expression data from adults' organ to define threshold
    organ_adults = {'Adipocytes';'Agland';'Colon';'Scord';'Skin';'Thyroidgland';'Urinarybladder'};
    idx_organs = findIndex(data_struct.tissues, organ_adults);
    ref_TPM = data_struct.levels(:,idx_organs);
    interGenes = intersect(model.genes, data_struct.genes);
    M_ExpValues = ref_TPM(findIndex(data_struct.genes, interGenes), :);
    percentile25 = prctile(M_ExpValues(:), 25);
    data_struct.threshold = percentile25;
else
    interGenes = intersect(model.genes, data_struct.genes);
    M_ExpValues = data_struct.levels(findIndex(data_struct.genes, interGenes), :);
    percentile25 = prctile(M_ExpValues(:), 25);
    data_struct.threshold = percentile25; % gene threshold
end

models = cell(numSamp, 1);
for i = 1:numSamp
    disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
    models{i} = ftINIT(prepData, data_struct.tissues{i}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), false, true);
    models{i}.id = data_struct.tissues{1,i};
end

models{end+1,1} = RBC;
save(['../models/', Age, Sex, '_models.mat'], 'models')

end

