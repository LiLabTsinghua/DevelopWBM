% Prepare SHAP data

Sex = 'Male';
Age = 'Adult';

% load WBM
modelPath = ['../Models/', Age, '_', Sex, '_WBM_withoutPSCM.mat'];
load(modelPath);

%% Parameter preparation
sex =  lower(Sex);
gender = lower(Sex);

% Allows fluctuation in the range of 0.8-1.2 for default parameters
min_fluc = 0.8; 
max_fluc = 1.2;

standardPhysiolDefaultParameters
IndividualParameters;

% bodyWeight
min_val = IndividualParameters.bodyWeight*min_fluc;
max_val = IndividualParameters.bodyWeight*max_fluc;
% Take a random sample of 100
bodyWeightRand = (max_val - min_val) .* rand(1, 100) + min_val;
bodyWeightRand = bodyWeightRand.';

% Height
min_val = IndividualParameters.Height*min_fluc;
max_val = IndividualParameters.Height*max_fluc;
% Take a random sample of 100
HeightRand = (max_val - min_val) .* rand(1, 100) + min_val;
HeightRand = HeightRand.';

% Hematocrit
min_val = IndividualParameters.Hematocrit*min_fluc;
max_val = IndividualParameters.Hematocrit*max_fluc;
% Take a random sample of 100
HematocritRand = (max_val - min_val) .* rand(1, 100) + min_val;
HematocritRand = HematocritRand.';

% HeartRate
min_val = IndividualParameters.HeartRate*min_fluc;
max_val = IndividualParameters.HeartRate*max_fluc;
% Take a random sample of 100
HeartRateRand = (max_val - min_val) .* rand(1, 100) + min_val;
HeartRateRand = HeartRateRand.';

% StrokeVolume
min_val = IndividualParameters.StrokeVolume*min_fluc;
max_val = IndividualParameters.StrokeVolume*max_fluc;
% Take a random sample of 100
StrokeVolumeRand = (max_val - min_val) .* rand(1, 100) + min_val;
StrokeVolumeRand = StrokeVolumeRand.';

% CardiacOutput
min_val = IndividualParameters.CardiacOutput*min_fluc;
max_val = IndividualParameters.CardiacOutput*max_fluc;
% Take a random sample of 100
CardiacOutputRand = (max_val - min_val) .* rand(1, 100) + min_val;
CardiacOutputRand = CardiacOutputRand.';

% CSFFlowRate
min_val = IndividualParameters.CSFFlowRate*min_fluc;
max_val = IndividualParameters.CSFFlowRate*max_fluc;
% Take a random sample of 100
CSFFlowRateRand = (max_val - min_val) .* rand(1, 100) + min_val;
CSFFlowRateRand = CSFFlowRateRand.';

% CSFBloodFlowRate
min_val = IndividualParameters.CSFBloodFlowRate*min_fluc;
max_val = IndividualParameters.CSFBloodFlowRate*max_fluc;
% Take a random sample of 100
CSFBloodFlowRateRand = (max_val - min_val) .* rand(1, 100) + min_val;
CSFBloodFlowRateRand = CSFBloodFlowRateRand.';

% UrFlowRate
min_val = IndividualParameters.UrFlowRate*min_fluc;
max_val = IndividualParameters.UrFlowRate*max_fluc;
% Take a random sample of 100
UrFlowRateRand = (max_val - min_val) .* rand(1, 100) + min_val;
UrFlowRateRand = UrFlowRateRand.';

% GlomerularFiltrationRate
min_val = IndividualParameters.GlomerularFiltrationRate*min_fluc;
max_val = IndividualParameters.GlomerularFiltrationRate*max_fluc;
% Take a random sample of 100
GlomerularFiltrationRateRand = (max_val - min_val) .* rand(1, 100) + min_val;
GlomerularFiltrationRateRand = GlomerularFiltrationRateRand.';
%% Calculation

for j=1:100
    WBMc = WBM;

    standardPhysiolDefaultParameters;
    IndividualParameters.bodyWeight = bodyWeightRand(j,1);
    IndividualParameters.Height = HeightRand(j,1);
    IndividualParameters.Hematocrit = HematocritRand(j,1);
    IndividualParameters.HeartRate = HeartRateRand(j,1);
    IndividualParameters.StrokeVolume = StrokeVolumeRand(j,1);
    IndividualParameters.CardiacOutput = CardiacOutputRand(j,1);
    IndividualParameters.CSFFlowRate = CSFFlowRateRand(j,1);
    IndividualParameters.CSFBloodFlowRate = CSFBloodFlowRateRand(j,1);
    IndividualParameters.UrFlowRate = UrFlowRateRand(j,1);
    IndividualParameters.GlomerularFiltrationRate = GlomerularFiltrationRateRand(j,1);

    IndividualParameters_ = IndividualParameters;
    % set PSCM to WBM. This function comes from CobraTool.
    % There would get some warning about reactions not in the model, 
    % including 'Gall_H2Ot[bdG]','Brain_DM_atp_c_' and 'Heart_DM_atp_c_'. 
    % Gall is not included by this stduy. As for 'Brain_DM_atp_c_' and 
    % 'Heart_DM_atp_c_' is no need to demand in this research. 
    % Thus, we ignore them.
    WBMc = physiologicalConstraintsHMDBbased(WBMc,IndividualParameters_);
    
    % set Diet constraints to WBM. This function comes from CobraTool
    EUAverageDietNew
    WBMc = setDietConstraints(WBMc, Diet);

    % set some more constraints
    WBMc = setSimulationConstraints(WBMc);
    WBMc.lb(strmatch('BBB_KYNATE[CSF]upt',WBMc.rxns)) = -1000000; %constrained uptake
    WBMc.lb(strmatch('BBB_LKYNR[CSF]upt',WBMc.rxns)) = -1000000; %constrained uptake
    WBMc.lb(strmatch('BBB_TRP_L[CSF]upt',WBMc.rxns)) = -1000000; %constrained uptake
    WBMc.ub(strmatch('Brain_EX_glc_D(',WBMc.rxns)) = -100; % currently -400 rendering many of the models to be infeasible in germfree state

    % get BMR value
    WBMc = releaseWBMConstraints(WBMc, Sex, Age);
    WBMc = changeObjective(WBMc, 'Whole_body_objective_rxn');
    try 
        [~, ~, modelIrrevFM, MinimizedFlux] = pFBA(WBMc, 'geneoption',0, 'tol',1e-6, 'skipclass', 1);
        if MinimizedFlux.stat == 1
            BMR_value = HumanBMR(modelIrrevFM, MinimizedFlux);
            Results_BMR(j,1) = BMR_value(1,1);
        else
            Results_BMR(j,1) = 0;
        end
    catch
        % disp('')
        Results_BMR(j,1) = 0;
    end

end

%% Output

