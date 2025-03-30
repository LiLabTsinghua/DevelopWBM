%% Comparation of organ-specific GEMs

% load AdultMale
modelPath = '../models/AdultMale_models.mat';
load(modelPath);
AdultMale = models;

% remove GEM of RBC organ GEM as is identifical in all WBMs.
for i=1:length(AdultMale)
    model = AdultMale{i,1};
    if contains(model.id, 'RBC')
        AdultMale(i) = [];
    else
        model.id = ['AdultMale_', model.id];
        AdultMale{i,1} = model;
    end
end
AdultMale = RemoveSexOrgans(AdultMale,'Male');

% load ElderlyMale
modelPath = '../models/ElderlyMale_models.mat';
load(modelPath);
ElderlyMale = models;

% remove GEM of RBC as it are identifical in all WBMs.
for i=1:length(ElderlyMale)
    model = ElderlyMale{i,1};
    if contains(model.id, 'RBC')
        ElderlyMale(i) = [];
    else
        model.id = ['ElderlyMale_', model.id];
        ElderlyMale{i,1} = model;
    end
end
ElderlyMale = RemoveSexOrgans(ElderlyMale,'Male');

% load AdultFemale
modelPath = '../models/AdultFemale_models.mat';
load(modelPath);
AdultFemale = models;

% remove GEM of RBC as it are identifical in all WBMs.
for i=1:length(AdultFemale)
    model = AdultFemale{i,1};
    if contains(model.id, 'RBC')
        AdultFemale(i) = [];
    else
        model.id = ['AdultFemale_', model.id];
        AdultFemale{i,1} = model;
    end
end
AdultFemale = RemoveSexOrgans(AdultFemale,'Female');

% load ElderlyFemale
modelPath = '../models/ElderlyFemale_models.mat';
load(modelPath);
ElderlyFemale = models;

% remove GEM of RBC as it are identifical in all WBMs.
for i=1:length(ElderlyFemale)
    model = ElderlyFemale{i,1};
    if contains(model.id, 'RBC')
        ElderlyFemale(i) = [];
    else
        model.id = ['ElderlyFemale_', model.id];
        ElderlyFemale{i,1} = model;
    end
end
ElderlyFemale = RemoveSexOrgans(ElderlyFemale,'Female');

% Assemble all organ-specific models
All_models = [AdultMale; ElderlyMale; AdultFemale; ElderlyFemale];

OrganName = {'Adipocytes', 'Agland', 'Brain', 'Colon', 'Heart', 'Kidney',...
             'Liver', 'Lung', 'Muscle', 'Pancreas', 'Scord', 'sIEC',...
             'Skin', 'Spleen', 'Stomach', 'Thyroidgland', 'Urinarybladder'};
% rank models with oragn name
RankModels = {};
for i=1:length(OrganName)
    organ = OrganName{1,i};
    for j=1:length(All_models)
        model = All_models{j,1};
        if contains(model.id, organ)
            RankModels{end+1} = model;
        end
    end
end
RankModels = RankModels.';


%% Compare all organ-specific models
compStruct = compareMultipleModels(RankModels); 

% HeatMap of similarity of models
clustergram(compStruct.structComp, 'Symmetric', false, 'Colormap', 'bone',...
    'RowLabels', compStruct.modelIDs, 'ColumnLabels', compStruct.modelIDs);


% tsne of similarity of models
rxn2Dmap = tsne(compStruct.reactions.matrix', 'Distance', 'hamming',...
    'NumDimensions', 2, 'Perplexity', 5);
% plot and label the GEMs in tSNE space
scatter(rxn2Dmap(:,1), rxn2Dmap(:,2));
hold on
text(rxn2Dmap(:,1), rxn2Dmap(:,2), compStruct.modelIDs);

tSNE_result = {};
tSNE_result(:,1) = compStruct.modelIDs;
tSNE_result(:,2) = num2cell(rxn2Dmap(:,1));
tSNE_result(:,3) = num2cell(rxn2Dmap(:,2));

tbl = array2table(tSNE_result, 'VariableNames', {'ModelNames', 'X', 'Y'});
%tbl.Properties.RowNames = Rows;
outputFileName = '../results/OrganModels_tsne.tsv';
writetable(tbl, outputFileName, 'Delimiter', '\t', 'WriteRowNames', true, 'FileType', 'text');


%% Subsystems coverage
useModels = compStruct.modelIDs;
keep = ismember(compStruct.modelIDs, useModels);
subMat = compStruct.subsystems.matrix(:, keep);
subCoverage = (subMat - mean(subMat, 2)) ./ mean(subMat, 2) * 100;
% select subsystems to include in plot
inclSub = any(abs(subCoverage) > 25, 2);
subNames = compStruct.subsystems.ID(inclSub);

% generate clustergram
cg = clustergram(subCoverage(inclSub,:), 'Colormap', redbluecmap,...
    'DisplayRange', 100, 'rowLabels', subNames, 'columnLabels',...
    useModels, 'Cluster', 'Column', 'ShowDendrogram', 'ON');

% Output the subsystems coverage result to "./Data"
SubOutput = subCoverage(inclSub,:);
Columns = useModels;
Columns = strrep(Columns, '-', '_');
Rows = subNames;

tbl = array2table(SubOutput, 'VariableNames', Columns);
tbl.Properties.RowNames = Rows;
outputFileName = '../results/SubCoverage.tsv';
writetable(tbl, outputFileName, 'Delimiter', '\t', 'WriteRowNames', true, 'FileType', 'text');

%% Comparation of metabolic function using full metabolic tasks
% check the metabolic tasks file whether in "./Data" folder.
MetaTasksPath = '../data/metabolicTasks/metabolicTasks_Full.txt';
% MetaTasks = readtable(DataPath, 'Delimiter', '\t', 'FileType', 'text');

% add boundary metabolites to each model for task checking.
for i=1:length(RankModels)
    model = RankModels{i,1};
    model = addBoundaryMets(model); % A function from RAVEN.
    RankModels{i,1} = model;
end

% Check the full metabolic tasks. It will take a while.
res_func = compareMultipleModels(RankModels, false, false, [], true, MetaTasksPath);

res_func.funcComp

% Identify which tasks differed among the GEMs
isDiff = ~all(res_func.funcComp.matrix == 0, 2) & ~all(res_func.funcComp.matrix == 1, 2);
diffTasks = res_func.funcComp.tasks(isDiff)

% visualize the matrix
spy(res_func.funcComp.matrix(isDiff,:), 10);

% apply some formatting changes
set(gca, 'XTick', 1:numel(useModels), 'XTickLabel', useModels, 'XTickLabelRotation', 90, ...
    'YTick', 1:numel(diffTasks), 'YTickLabel', diffTasks, 'YAxisLocation', 'right');
xlabel(gca, '');

% Outputhe Function tasks result to "./Data"
FunctionOutput = res_func.funcComp.matrix(isDiff,:);
Columns = useModels;
Columns = strrep(Columns, '-', '_');
Rows = diffTasks;

tbl = array2table(FunctionOutput, 'VariableNames', Columns);
tbl.Properties.RowNames = Rows;
outputFileName = '../results/FullTasks_Results.tsv';
writetable(tbl, outputFileName, 'Delimiter', '\t', 'WriteRowNames', true, 'FileType', 'text');
