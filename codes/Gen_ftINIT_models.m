function [eGenes, models] = Gen_ftINIT_models(modelName,GeneExpData)

if nargin < 2
    GeneExpData = 'Hart2015';
end

%% load data
if strcmp(GeneExpData, 'Hart2015')
    GeneExpDataFloder = '../data/eGenesData/Hart2015_RNAseq.txt';
    gtex_data = readtable(GeneExpDataFloder);
    [~, n] = size(gtex_data);
    numSamp = n-1;

    gtex_data(1:5, 1:5)

    data_struct.genes = gtex_data{:, 1};
    data_struct.tissues = gtex_data.Properties.VariableNames(2:n);
    data_struct.levels = gtex_data{:, 2:n};
    data_struct.threshold = 1;
    data_struct;

elseif strcmp(GeneExpData, 'DepMap')
    GeneExpDataFloder = '../data/eGenesData/DepMap_RNAseq_data.mat';
    load(GeneExpDataFloder);

    data_struct.genes = depmap.genes;
    data_struct.tissues = depmap.cellID;
    data_struct.levels = depmap.tpm;
    data_struct.threshold = 1;
    
else
    error('Please select either Hart2015 or DepMap as the dataset!')
    
end

numSamp = length(data_struct.tissues);

%% Generate cell-specific GEMs
if strcmp(modelName, 'Human2')
    % load Human2
    load('../models/Human2.mat');
    % model = ihuman;

    % prepare Human-GEM
    taskStruct = '../data/metabolicTasks/metabolicTasks_Essential.txt';
    
    prepData = prepHumanModelForftINIT(model, false, taskStruct, '../data/reactions.tsv');

    % generate context-specific models
    models = cell(numSamp, 1);
    for i = 1:numSamp
        disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
        % use '1+1' mode to reduce number of reactions without GPR
        models{i} = ftINIT(prepData, data_struct.tissues{i}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), true, true);
        models{i}.id = data_struct.tissues{1,i};
    end
    
elseif strcmp(modelName, 'Human1')
    % Human1 ftINIT-Hart2015
    load('../models/Human1.mat')
    % model = ihumanGEM;

    %load tasks
    essentialTasksFilePath = '../data/metabolicTasks/metabolicTasks_essential_human1.txt';
    taskStruct = parseTaskList(essentialTasksFilePath);

    %Spontaneous reactions:
    spontRxnNames = {...
        'HMR_4740';...
        'HMR_4250';...
        'HMR_6875';...
        'HMR_6876';...
        'HMR_4840';...
        'HMR_4771';...
        'HMR_6997';...
        'HMR_7008';...
        'HMR_7011';...
        'HMR_7015';...
        'HMR_7016';...
        'HMR_5127'};

    %remove some reactions often not used from the model to speed up calculations

    %drug reactions:
    m = removeDrugReactions(model);
    %remove all AA triplet rxns, i.e. rxns of the type
    %2 H2O[c] + Tryptophanyl-Glycyl-Aspartate[c] <=> aspartate[c] + glycine[c] + tryptophan[c]
    %This is only used for import of such compounds, and if you don't have that in your "medium",
    %which we usually don't, these are pretty pointless.
    AATriplets = getAATripletReactions(m,false);
    %length(AATriplets)%735
    m = removeReactions(m, AATriplets);

    %Remove the duplicate complex I reaction with ROS - it only causes problems to include it
    m = removeReactions(m, {'CYOOm3i'});

    %prepare some reactions that can always be on:

    %protein reactions creation/degradation
    proteinRxns = { ...
        'HMR_5155'; ...
        'HMR_5156'; ...
        'HMR_5161'; ...
        'HMR_5167'; ...
        'HMR_5168'; ...
        'HMR_5169'; ...
        'HMR_5170'; ...
        'HMR_5171'; ...
        'HMR_5172'; ...
        'HMR_5174'; ...
        'HMR_5260'; ...
        'HMR_5262'; ...
        'HMR_5264'; ...
        'HMR_5266'; ...
        'HMR_5267'; ...
        'HMR_5268'; ...
        'HMR_5269'; ...
        'HMR_5270'; ...
        'HMR_5271'; ...
        'HMR_5273'; ...
        'HMR_5275'; ...
        'HMR_5277'; ...
        'HMR_5279'; ...
        'HMR_5281'; ...
        'HMR_5283'; ...
        'HMR_5291'; ...
        'HMR_9817'; ...
        'HMR_9818'};

    %reactions that just pool metabolites
    poolRxns = { ...
        'HMR_0011'; ...
        'HMR_0012'; ...
        'HMR_0477'; ...
        'HMR_5233'; ...
        'HMR_5234'; ...
        'HMR_5238'; ...
        'HMR_5239'; ...
        'HMR_5243'; ...
        'HMR_5244'; ...
        'HMR_5247'; ...
        'HMR_9022'; ...
        'HMR_0015'; ...
        'HMR_0016'; ...
        'HMR_0017'; ...
        'HMR_0033'; ...
        'HMR_0035'; ...
        'HMR_0036'; ...
        'HMR_0037'; ...
        'HMR_0038'; ...
        'HMR_0062'; ...
        'HMR_0063'; ...
        'HMR_0064'; ...
        'HMR_0065'; ...
        'HMR_3082'};

    %Radical reactions.
    %I have assumed that these are spontaneous, but I don't actually have any proof, so we leave them
    %out for now. These could be investigated more - if found spontaneous - set the spontaneous flag instead.
    customRxnsToIgnore = unique([proteinRxns;poolRxns]);

    % m = simplifyModel(m);
    convertGenes = false;
    prepDataHumanGEM = prepINITModel(m, taskStruct, spontRxnNames, convertGenes, customRxnsToIgnore, 's');
    prepData = prepDataHumanGEM;

    models = cell(numSamp, 1);
    for i = 1:numSamp
        disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
        % use '1+1' mode to reduce number of reactions without GPR
        models{i} = ftINIT(prepData, data_struct.tissues{i}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), false, true);
        models{i}.id = data_struct.tissues{1,i};
    end

elseif strcmp(modelName, 'Recon3D')
    % Recon3D ftINIT-Hart2015
    load('../models/Recon3D_301.mat')
    model = Recon3D;

    % remove some unnecessary fields in Recon3D
    model = rmfield(model,{'rules','rxnECNumbers','rxnKEGGID','rxnCOG', ...
        'rxnKeggOrthology','metSmiles','metHMDBID','metInChIString', ...
        'metKEGGID','metPubChemID','rxnNotes','rxnConfidenceScores', ...
        'metCHEBIID','metPdMap','metReconMap','rxnReconMap','csense','osense'});

    % add "Biomass[c]" metabolite as a product in the "biomass_reaction" rxn
    [~,rxn_ind] = ismember('biomass_reaction',model.rxns);
    [~,met_ind] = ismember('Temp001[c]',model.mets);  % Temp001[c] is the biomass metabolite
    model.S(met_ind,rxn_ind) = 1;

    % constrain the "HMR_biomass_Renalcancer" rxn so it cannot be used instead
    model = setParam(model,'eq','HMR_biomass_Renalcancer',0);

    % rename biomass sink reaction
    model.rxns(ismember(model.rxns,'DM_HMR_biomass_renalcancer')) = {'EX_biomass'};

    % several grRules are ambiguous due to lack of parentheses. Update
    % these rules, assuming that the order of operations begins with ANDs,
    % followed by ORs. Note: these ambiguous rules were identified using
    % the "cleanModelGeneRules" function.
    badRules = {'(1629.1) and (594.1 and 593.1) and (1738.1) or (1629.1) and (1738.1) and (594.2 and 593.1)';
        '(1738.1 and 8050.1) and (5161.1 and 5162.1) and (1737.1) or (1738.1 and 8050.1) and (5160.1 and 5162.1) and (1737.1)';
        '(2194.1) or (2194.1 and 79071.1) and (60481.1) and (54898.1) ';
        '(5476.1) and (2720.1) and (2588.1) and (4758.1) and (5660.1) or (5476.1) and (2720.1) and (2588.1) and (4758.1) and (2760.1) or (5476.1) and (2720.1) and (5660.1) or (5476.1) and (2720.1) and (2760.1)';
        '130.1 or 125.1 or 125.1 and 124.1 or 128.1 or 126.1 or 284273.2 or 137872.1 or 127.1 or 124.1 or (126.1 and 124.1) or 284273.1 or (125.1 and 126.1) or 131.1'};

    newRules = {'(1629.1 and 594.1 and 593.1 and 1738.1) or (1629.1 and 1738.1 and 594.2 and 593.1)';
        '(1738.1 and 8050.1 and 5161.1 and 5162.1 and 1737.1) or (1738.1 and 8050.1 and 5160.1 and 5162.1 and 1737.1)';
        '2194.1 or (2194.1 and 79071.1 and 60481.1 and 54898.1)';
        '(5476.1 and 2720.1 and 2588.1 and 4758.1 and 5660.1) or (5476.1 and 2720.1 and 2588.1 and 4758.1 and 2760.1) or (5476.1 and 2720.1 and 5660.1) or (5476.1 and 2720.1 and 2760.1)';
        '130.1 or 125.1 or (125.1 and 124.1) or 128.1 or 126.1 or 284273.2 or 137872.1 or 127.1 or 124.1 or (126.1 and 124.1) or 284273.1 or (125.1 and 126.1) or 131.1'};

    [isBad,ruleInd] = ismember(model.grRules,badRules);
    model.grRules(isBad) = newRules(ruleInd(isBad));

    % add some necessary fields
    model.id = 'Recon3D';

    % change peroxisome compartment abbreviation from [x] to [p]
    model.mets = regexprep(model.mets,'\[x\]$','[p]');
    if isfield(model,'comps')
        model.comps(ismember(model.comps,'x')) = {'p'};
    end

    % change extracellular compartment abbreviation from [e] to [s]
    model.mets = regexprep(model.mets,'\[e\]$','[s]');
    if isfield(model,'comps')
        model.comps(ismember(model.comps,'e')) = {'s'};
    end

    % add "metComps", "compNames", and "rev" fields
    compAbbrevName = {'c' 'Cytosol'
                      's' 'Extracellular'
                      'g' 'Golgi apparatus'
                      'i' 'Inner mitochondria';
                      'l' 'Lysosome';
                      'm' 'Mitochondria';
                      'n' 'Nucleus';
                      'r' 'Endoplasmic reticulum';
                      'p' 'Peroxisome'};
    model = addMetCompsField(model);
    [~,ind] = ismember(model.comps,compAbbrevName(:,1));
    model.compNames = compAbbrevName(ind,2);
    model.rev = double(model.lb < 0 & model.ub > 0);

    % add boundary metabolites to the model
    % NOTE: the second input to this function is FALSE, meaning that boundary
    % metabolites will be added even to reactions that involve metabolites in
    % compartments other than extracellular (e.g., sink and demand reactions).
    model = addBoundaryMets(model,false);

    % DELETE all sink and demand reactions in the model. If not, then these can
    % erroneously be added back during the tINIT model generation process.
    sink_dm_ind = find(startsWith(model.rxns,{'sink_','DM_'}));
    model = removeReactionsFull(model,sink_dm_ind);

    % translate model genes from Entrez (NCBI) to Ensembl IDs (ENSG)
    [grRules,genes,rxnGeneMat] = translateGrRules(model.grRules,'ENSG','Entrez');
    model.grRules = grRules;
    model.genes = genes;
    model.rxnGeneMat = rxnGeneMat;

    % Run some preliminary steps that will allow us to skip some pre-processing
    % steps in the tINIT algorithm, greatly reducing the overall run time.
    [~,deletedDeadEndRxns] = simplifyModel(model,true,false,true,true,true);
    cModel = removeReactions(model,deletedDeadEndRxns,false,true);

    model1 = simplifyModel(cModel);
    taskStruct = parseTaskList('../data/metabolicTasks/Recon3D_essential_tasks.txt');

    spontRxnNames = {};
    convertGenes = false; 
    customRxnsToIgnore = {}; 
    extComp = 's';

    prepDataRecon3D = prepINITModel(model1, taskStruct, spontRxnNames, convertGenes, customRxnsToIgnore, extComp);
    prepData = prepDataRecon3D;

    % run ftINIT
    models = cell(numSamp, 1);
    for i = 1:numSamp
        disp(['Model: ' num2str(i) ' of ' num2str(numSamp)])
        % use '1+1' mode to reduce number of reactions without GPR
        models{i} = ftINIT(prepData, data_struct.tissues{i}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), false, true);
        models{i}.id = data_struct.tissues{1,i};
    end
end

%% Predict essential genes
if strcmp(modelName, 'Human2')
    model = addBoundaryMets(model);
    taskStruct = parseTaskList(taskStruct);
end

% iterate through the different models
parfor i = 1:length(models)
    m = models{i};
    m = addBoundaryMets(m);
    
    % determine essential genes for each task
    [~,essentialGenes] = checkTasksGenes(m,[],false,false,true,taskStruct);
    
    % collect results
    tissues{i,1} = m.id;
    geneList{i,1} = m.genes;
    allEssentials{i,1} = essentialGenes;
end

% gather results into eGenes structure
eGenes = {};
eGenes.taskList = {taskStruct(:).description}';
eGenes.tissues = tissues;
eGenes.geneList = geneList;
eGenes.essentialGenes = allEssentials;
eGenes.refModel = model;

end


