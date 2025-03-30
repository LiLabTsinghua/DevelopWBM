function [results] = FVA_FluxDistribution(modelName)

% Run flux distribution for Human-GEM and ecHuman-GEM. The name of model
% and ecmodel should be defined.

if strcmp(modelName, 'Human1')
    % load models
    load(['../models/', modelName, '.mat']);
    load(['../models/ec', modelName, '.mat']);
    
    % set media constrains
    model   = setHamsMedium(model,false);
    ecModel = setHamsMedium(ecModel,false);
    
    % set prot_pool_exchange bounds for ecHuman
    ecModel = changeRxnBounds(ecModel, 'prot_pool_exchange', -20.83, 'l');
    
    CsourceUptk = 'HMR_9034';

elseif strcmp(modelName, 'Human2')
    % load models
    load(['../models/', modelName, '.mat']);
    load(['../models/ec', modelName, '.mat']);
    
    % set media constrains
    model   = setHamsMedium(model,false);
    ecModel = setHamsMedium(ecModel,false);
    
    % set prot_pool_exchange bounds for ecHuman
    ecModel = changeRxnBounds(ecModel, 'prot_pool_exchange', -20.83, 'l');
    
    CsourceUptk = 'MAR09034';
else
    disp("The model is neither Human1 nor Human2. You can upload the model and define the medium constraints yourself.")
end

[FVA_Dists,indexes,~,~] = comparativeFVA(model,ecModel,CsourceUptk,false,1E-8);

% save results
mkdir('../results/FVA')
variables = {'rxns' 'formulas' 'model_ranges' 'ecModel_ranges' 'subSystems'};
formulas  = constructEquations(model,indexes); % a function from RAVEN to construct equations;
results   = table(model.rxns(indexes),formulas,FVA_Dists{1},FVA_Dists{2},model.subSystems(indexes),'VariableNames',variables);
writetable(results,['../results/FVA/FVA_comp_' modelName],'Delimiter','\t','QuoteStrings',false)
% writetable(results)

end

