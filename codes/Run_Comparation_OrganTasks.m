% This script describes the process of checking organ metabolic tasks. 
% It can support the metabolic task checks for two models, Human2 and
% Recon3D, in four organs: adipocytes, liver, heart, and kidneys. 
% Before running the code, you need to ensure that the models are stored 
% in the "./models/" folder, and the relevant metabolic tasks are placed 
% in the "./data/eGenesData/metabolicTasks/OrganTasks/" folder.

% You need to define Age and Sex group
Age = 'Adult'; % Only Adult, Elderly GEMs have no model to compare.
Sex = 'Male'; % or Female

ModelNames = {'Human2'; 'Recon3D'};
Organs = {'Liver'; 'Heart'; 'Kidney'; 'Adipocytes'};

TaskRatio = {};
for i=1:length(ModelNames)
    modelName = ModelNames{i,1};
    if strcmp(modelName, 'Human2')
        % load models
        load(['../models/', Age, Sex, '_models.mat']);
        % check tasks
        for j=1:length(Organs)
            organ = Organs{j,1};
            for m =1:length(models)
                if strcmp(models{m,1}.id, organ)
                    model = models{m,1};
                    bModel = addBoundaryMets(model);
                    [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(bModel, organ);
                    count_ones = sum(taskReport.ok == 1);
                    SuccessTasks_ratio = count_ones/numel(taskReport.ok);
                    fprintf([organ,': Ratio of successful tasks: %d\n'], SuccessTasks_ratio);
                end
            end
            TaskRatio{1,j} = round(SuccessTasks_ratio, 4); 
        end
    
    elseif strcmp(modelName, 'Recon3D')
        % load models
        if strcmp(Sex, 'Male')
            modelPath = '../models/Harvey_Harvetta/OrganAtlas_Harvey.mat';
        elseif strcmp(Sex, 'Female')
            modelPath = '../models/Harvey_Harvetta/OrganAtlas_Harvetta.mat';
        end
        load(modelPath);
        
        varNames = who;
        for var = 1:length(varNames)
            if startsWith(varNames{var}, 'OrganCompendium')
                models = eval(varNames{var});
            end
        end
        % check tasks
        for j=1:length(Organs)
            organ = Organs{j,1};
            for m =1:length(models)
                if isfield(models, organ)
                    model = models.(organ).modelAllComp;
                    type = 'Recon';
                    [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, organ, type);
                    count_ones = sum(taskReport.ok == 1);
                    SuccessTasks_ratio = count_ones/numel(taskReport.ok);
                    fprintf([organ,': Ratio of successful tasks: %d\n'], SuccessTasks_ratio);
                end
            end
            TaskRatio{2,j} = round(SuccessTasks_ratio, 4); 
        end
    end
end

tbl = cell2table(TaskRatio, 'VariableNames', Organs);
tbl.Properties.RowNames = ModelNames;
fileNames = ['../results/metabolicTasks_comparation/',Age,Sex, '_MetabolicTasks.tsv'];
writetable(tbl, fileNames, 'Delimiter', '\t', 'FileType', 'text', 'WriteRowNames', true);
        
    
