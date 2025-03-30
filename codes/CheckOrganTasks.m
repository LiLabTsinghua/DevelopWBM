function [taskReport, essentialRxns, taskStructure, essentialFluxes] = CheckOrganTasks(model, Organ, type)


if nargin < 3
    type = 'Human';
end

if nargin < 2 
    EM = 'Models and organs need to be provided£¡';
    dispEM(EM);
end

essentialTasksFilePath = ['../data/metabolicTasks/OrganTasks/', Organ, 'Tasks_', type, '.txt'];
taskStruct = parseTaskList(essentialTasksFilePath);

if strcmp(type, 'Human')
    bModel = closeModel(model);
    [taskReport, essentialRxns, taskStructure, essentialFluxes] = checkTasks(bModel,[],false,false,false,taskStruct);
elseif strcmp(type, 'Recon')
    model.compNames = {'Extracellular'; 'Peroxisome'; 'Mitochondria';...
                       'Cytosol'; 'Lysosome'; 'Endoplasmic reticulum';...
                       'Golgi apparatus'; 'Nucleus'; 'Inner mitochondria';...
                       'luSI-sIEC'; 'sIEC-bp'; 'Blood, portal vein'; ...
                       'Lumen, small intestine';'Blood, circulation';...
                       'luLI-Colon'; 'Colon-bp'; 'Lumen, large intestine';...
                       'Lumen'; 'Feces'; 'bd-Liver'; 'bp-Liver'; 'Bile duct'; 'lu-Pancreas';...
                       'Pancreas-bp'; 'Kidney-bc'; 'Urine'; 'Scord-bp'; 'Cerebrospinal fluid';...
                       'Air'; 'Diet'; 'air-Lung'; 'sw-Skin'; 'Sweat'};
    model.comps = {'e'; 'x'; 'm'; 'c'; 'l'; 'r'; 'g'; 'n'; 'i';  'luI';...
                   'bpI'; 'bp'; 'luSI'; 'bc'; 'luC'; 'bpC'; 'luLI';...
                   'lu'; 'fe'; 'bdL'; 'bpL'; 'bd'; 'luP'; 'bpP'; 'bcK';...
                   'u'; 'bpS'; 'csf'; 'a'; 'd'; 'aL'; 'swS'; 'sw'};
    CopyComps = {'[e]'; '[x]'; '[m]'; '[c]'; '[l]'; '[r]'; '[g]'; '[n]'; '[i]';  '[luI]';...
                   '[bpI]'; '[bp]'; '[luSI]'; '[bc]'; '[luC]'; '[bpC]'; '[luLI]';...
                   '[lu]'; '[fe]'; '[bdL]'; '[bpL]'; '[bd]'; '[luP]'; '[bpP]'; '[bcK]';...
                   '[u]'; '[bpS]'; '[csf]'; '[a]'; '[d]'; '[aL]'; '[swS]'; '[sw]'};
    metComps = zeros(size(model.mets, 1), 1);
    % Newmets = cell(size(model.mets));
    for i = 1:numel(model.mets)
        met = model.mets{i};
            
        for j = 1:numel(CopyComps)
            if contains(met, CopyComps{j})
                metComps(i) = j;
                % CurrentMet = strrep(met, CopyComps{j}, '');
                break; 
            end
        end
        % Newmets{i} = CurrentMet;
    end
    
    if any(metComps == 0)
        EM = 'Some metabolites do not belong to any compartment£¡';
        dispEM(EM);
    else
        model.metComps = metComps;
        % model.mets = Newmets;
        bModel = closeModel(model);
        [taskReport, essentialRxns, taskStructure, essentialFluxes] = checkTasks(bModel,[],false,false,false,taskStruct);
    end
    
end

