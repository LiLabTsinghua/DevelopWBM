% Run flux distribution for Human-GEM and ecHuman-GEM. The name of model
% and ecmodel should be defined.
model_id = "Human1";
% model = "Human2";

if model_id == 'Human1'
    model_name = 'HumanGEM_1_3_FVA.mat';
    ecmodel_name = 'ec_HumanGEM_1_3_FVA.mat';
    
    load(['../Models/', model_name])
    load(['../Models/', ecmodel_name]);
    
    % set media constrains
    model   = setHamsMedium(model,false);
    ecModel = setHamsMedium(ecModel,false);
    
    % set prot_pool_exchange bounds for ecHuman
    ecModel = changeRxnBounds(ecModel, 'prot_pool_exchange', -20.83, 'l');
    
    CsourceUptk = 'HMR_9034';

elseif model_id == 'Human2'
    model_name = 'Human-GEM_FVA.yml';
    ecmodel_name = 'Human-GEM_EC_FVA.mat';
    
    model = importYaml(['../Models/', model_name]);
    load(['../Models/', ecmodel_name]);
    
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


% model_name = 'Human-GEM_FVA.yml'; % Human2
model_name = 'HumanGEM_1_3_FVA.mat'; % Human1
model = importYaml(['../Models/', model_name]);
load(['../Models/', model_name])

% ecmodel_name = 'Human-GEM_EC_FVA.mat'; % Human2
ecmodel_name = 'ec_HumanGEM_1_3_FVA.mat'; % Human2
load(['../Models/', ecmodel_name]);

% set media constrains
model   = setHamsMedium(model,false);
ecModel = setHamsMedium(ecModel,false);

% set prot_pool_exchange bounds for ecHuman
ecModel = changeRxnBounds(ecModel, 'prot_pool_exchange', -20.83, 'l');

% run comparativeFVA to get feasible flux distribution
CsourceUptk = 'MAR09034';
[FVA_Dists,indexes,~,~] = comparativeFVA(model,ecModel,CsourceUptk,false,1E-8);

% save results
mkdir('../../Results/FVA')
variables = {'rxns' 'formulas' 'model_ranges' 'ecModel_ranges' 'subSystems'};
formulas  = constructEquations(model,indexes); % a function from RAVEN to construct equations;
results   = table(model.rxns(indexes),formulas,FVA_Dists{1},FVA_Dists{2},model.subSystems(indexes),'VariableNames',variables);
% writetable(results,['../../Results/FVA/FVA_comp_' cellLine],'Delimiter','\t','QuoteStrings',false)
writetable(results)

