% This script describes the process of essential gene validation. This code
% can support the essential gene validation for Human2, Human1, and 
% Recon3D. Before running the code, it is necessary to ensure that the 
% model has been placed in the "./models/" folder, and the relevant dataset
% is in the "./data/eGenesData/" folder.

% Generic GEM
modelName = 'Human2'; % or 'Human1', 'Recon3D'
% DataSet
GeneExpData = 'Hart2015'; % or 'DepMap'

% Generate cell-specific models and predict essential genes.
[eGenes, models] = Gen_ftINIT_models(modelName,GeneExpData);

% Validation of essential genes
results = evaluateEssentialGenes(eGenes, GeneExpData);

% Output results
variableNames = results(1, 2:end);
data = results(2:end, :);
rowNames = data(:, 1);
tbl = cell2table(data(:, 2:end), 'VariableNames', variableNames);
tbl.Properties.RowNames = rowNames;
fileNames = ['../results/eGenes/', modelName, '_', GeneExpData, '_essentialGenes_ftINIT.tsv'];
writetable(tbl, fileNames, 'Delimiter', '\t', 'FileType', 'text', 'WriteRowNames', true);