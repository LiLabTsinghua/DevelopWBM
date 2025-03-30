% check metabolic tasks for organ-specific GEMs

%% Liver tasks
% load metabolic tasks
essentialTasksFilePath = '../Data/Tasks/LiverTasks_Human.txt';
% essentialTasksFilePath = '../Data/Tasks/Liver_tasks_1.txt';
taskStruct = parseTaskList(essentialTasksFilePath);

% load models
modelPath = '../Models/AdultMale_models_(threshold=50th)_1004.mat';
load(modelPath);

% taskStruct_test = parseTaskList('./data/metabolicTasks/metabolicTasks_Essential.txt');
% bModel = closeModel(model);
% [taskReport, essentialRxnMat, ~, essentialFluxes] = checkTasks(bModel,[],true,true,false,taskStruct);
Organ = 'Liver';

% Human Liver
for i =1:length(models)
    if strcmp(models{i,1}.id, 'Liver')
        model = models{i,1};
        % bModel = closeModel(model);
        [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, Organ);
        count_ones = sum(taskReport.ok == 1);
        SuccessTasks_ratio = count_ones/numel(taskReport.ok);
        fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
    end
end



% Recon Liver
% modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvey.mat';
modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvetta.mat';
load(modelPath);

Organ = 'Liver';
type = 'Recon';

varNames = who;
for i = 1:length(varNames)
    if startsWith(varNames{i}, 'OrganCompendium')
        models = eval(varNames{i});
    end
end

for i =1:length(models)
    if isfield(models, 'Liver')
        model = models.Liver.modelAllComp;
        % bModel = closeModel(model);
        [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, Organ, type);
        count_ones = sum(taskReport.ok == 1);
        SuccessTasks_ratio = count_ones/numel(taskReport.ok);
        fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
    end
end

%% Heart tasks
% load metabolic tasks
essentialTasksFilePath = '../Data/Tasks/HeartTasks.txt';
% essentialTasksFilePath = '../Data/Tasks/Liver_tasks_1.txt';
taskStruct = parseTaskList(essentialTasksFilePath);

% load models
modelPath = '../Models/AdultMale_models_(threshold=50th)_1004.mat';
load(modelPath);

% taskStruct_test = parseTaskList('./data/metabolicTasks/metabolicTasks_Essential.txt');
% bModel = closeModel(model);
% [taskReport, essentialRxnMat, ~, essentialFluxes] = checkTasks(bModel,[],true,true,false,taskStruct);
Organ = 'Heart';

% Human Liver
for i =1:length(models)
    if strcmp(models{i,1}.id, Organ)
        model = models{i,1};
        % bModel = closeModel(model);
        [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, Organ);
        count_ones = sum(taskReport.ok == 1);
        SuccessTasks_ratio = count_ones/numel(taskReport.ok);
        fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
    end
end


% Recon Liver
% modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvey.mat';
modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvetta.mat';
load(modelPath);

Organ = 'Heart';
type = 'Recon';

varNames = who;
for i = 1:length(varNames)
    if startsWith(varNames{i}, 'OrganCompendium')
        models = eval(varNames{i});
    end
end

for i =1:length(models)
    if isfield(models, 'Heart')
        model = models.Heart.modelAllComp;
        % bModel = closeModel(model);
        [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, Organ, type);
        count_ones = sum(taskReport.ok == 1);
        SuccessTasks_ratio = count_ones/numel(taskReport.ok);
        fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
    end
end

%% Kidney tasks
% load metabolic tasks
essentialTasksFilePath = '../Data/Tasks/KidneyTasks.txt';
% essentialTasksFilePath = '../Data/Tasks/Liver_tasks_1.txt';
taskStruct = parseTaskList(essentialTasksFilePath);

% load models
modelPath = '../Models/AdultMale_models_(threshold=50th)_1004.mat';
load(modelPath);

% taskStruct_test = parseTaskList('./data/metabolicTasks/metabolicTasks_Essential.txt');
% bModel = closeModel(model);
% [taskReport, essentialRxnMat, ~, essentialFluxes] = checkTasks(bModel,[],true,true,false,taskStruct);
Organ = 'Kidney';

% Human Liver
for i =1:length(models)
    if strcmp(models{i,1}.id, Organ)
        model = models{i,1};
        % bModel = closeModel(model);
        [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, Organ);
        count_ones = sum(taskReport.ok == 1);
        SuccessTasks_ratio = count_ones/numel(taskReport.ok);
        fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
    end
end


% Recon Kidney
% modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvey.mat';
modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvetta.mat';
load(modelPath);

Organ = 'Kidney';
type = 'Recon';

varNames = who;
for i = 1:length(varNames)
    if startsWith(varNames{i}, 'OrganCompendium')
        models = eval(varNames{i});
    end
end

for i =1:length(models)
    if isfield(models, 'Kidney')
        model = models.Kidney.modelAllComp;
        % bModel = closeModel(model);
        [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, Organ, type);
        count_ones = sum(taskReport.ok == 1);
        SuccessTasks_ratio = count_ones/numel(taskReport.ok);
        fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
    end
end


        
%% Adipocytes tasks
% load metabolic tasks
essentialTasksFilePath = '../Data/Tasks/AdipocytesTasks.txt';
% essentialTasksFilePath = '../Data/Tasks/Liver_tasks_1.txt';
taskStruct = parseTaskList(essentialTasksFilePath);

% load models
modelPath = '../Models/AdultMale_models_(threshold=50th)_1004.mat';
load(modelPath);

% taskStruct_test = parseTaskList('./data/metabolicTasks/metabolicTasks_Essential.txt');
% bModel = closeModel(model);
% [taskReport, essentialRxnMat, ~, essentialFluxes] = checkTasks(bModel,[],true,true,false,taskStruct);
Organ = 'Adipocytes';

% Human Liver
for i =1:length(models)
    if strcmp(models{i,1}.id, Organ)
        model = models{i,1};
        % bModel = closeModel(model);
        [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, Organ);
        count_ones = sum(taskReport.ok == 1);
        SuccessTasks_ratio = count_ones/numel(taskReport.ok);
        fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
    end
end



% Recon Adipocytes
% modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvey.mat';
modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvetta.mat';
load(modelPath);

Organ = 'Adipocytes';
type = 'Recon';

varNames = who;
for i = 1:length(varNames)
    if startsWith(varNames{i}, 'OrganCompendium')
        models = eval(varNames{i});
    end
end

for i =1:length(models)
    if isfield(models, 'Adipocytes')
        model = models.Adipocytes.modelAllComp;
        % bModel = closeModel(model);
        [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, Organ, type);
        count_ones = sum(taskReport.ok == 1);
        SuccessTasks_ratio = count_ones/numel(taskReport.ok);
        fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
    end
end


%% All tasks
% load metabolic tasks
essentialTasksFilePath = '../Data/Tasks/All_tasks.txt';
% essentialTasksFilePath = '../Data/Tasks/All_tasks_human1.txt';
% essentialTasksFilePath = '../Data/Tasks/Liver_tasks_1.txt';
taskStruct = parseTaskList(essentialTasksFilePath);

bModel = closeModel(model);
[taskReport, essentialRxnMat, ~, essentialFluxes] = checkTasks(bModel,[],true,true,false,taskStruct);
count_ones = sum(taskReport.ok == 1);
SuccessTasks_ratio = count_ones/numel(taskReport.ok);
fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);

OrgantaskReport = struct();
OrgantaskReport.organs = cell(numel(models),1);
OrgantaskReport.tasks = cell(numel(taskStruct),1);
for i=1:length(taskStruct)
    OrgantaskReport.tasks{i,1} = taskStruct(i).description;
end
OrgantaskReport.checks = zeros(numel(taskStruct), numel(models));


for i =1:length(models)
    
    model = models{i,1};
    bModel = closeModel(model);
    [taskReport, essentialRxnMat, ~, essentialFluxes] = checkTasks(bModel,[],false,false,false,taskStruct);
    OrgantaskReport.checks(:,i) = taskReport.ok;
    OrgantaskReport.organs{i,1} = model.id;
    
%     count_ones = sum(taskReport.ok == 1);
%     SuccessTasks_ratio = count_ones/numel(taskReport.ok);
%     fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
    
end



%%
% load models
modelPath = '../Models/ElderlyFemale_models_(threshold=25th)_1004.mat';
load(modelPath);

% taskStruct_test = parseTaskList('./data/metabolicTasks/metabolicTasks_Essential.txt');
% bModel = closeModel(model);
% [taskReport, essentialRxnMat, ~, essentialFluxes] = checkTasks(bModel,[],true,true,false,taskStruct);
Organs = {'Liver'; 'Heart'; 'Kidney'; 'Adipocytes'};
TaskRatio = {};

% Human Liver
for i=1:length(Organs)
    organ = Organs{i,1};
    for j =1:length(models)
        if strcmp(models{j,1}.id, organ)
            model = models{j,1};
            % bModel = closeModel(model);
            [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, organ);
            count_ones = sum(taskReport.ok == 1);
            SuccessTasks_ratio = count_ones/numel(taskReport.ok);
            fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
        end
    end
    TaskRatio{i,1} = round(SuccessTasks_ratio, 4); 
end


%%
% load models
% Recon Adipocytes
modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvey.mat';
% modelPath = '../Models/Harvey_Harvetta/OrganAtlas_Harvetta.mat';
load(modelPath);

% taskStruct_test = parseTaskList('./data/metabolicTasks/metabolicTasks_Essential.txt');
% bModel = closeModel(model);
% [taskReport, essentialRxnMat, ~, essentialFluxes] = checkTasks(bModel,[],true,true,false,taskStruct);
Organs = {'Liver'; 'Heart'; 'Kidney'; 'Adipocytes'};
TaskRatio = {};

type = 'Recon';

varNames = who;
for i = 1:length(varNames)
    if startsWith(varNames{i}, 'OrganCompendium')
        models = eval(varNames{i});
    end
end

for i=1:length(Organs)
    organ = Organs{i,1};
    for j =1:length(models)
        if isfield(models, organ)
            model = models.(organ).modelAllComp;
            [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, organ, type);
            count_ones = sum(taskReport.ok == 1);
            SuccessTasks_ratio = count_ones/numel(taskReport.ok);
            fprintf('Ratio of successful tasks: %d\n', SuccessTasks_ratio);
        end
    end
    TaskRatio{i,1} = round(SuccessTasks_ratio, 4); 
end