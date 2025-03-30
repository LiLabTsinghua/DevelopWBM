% Run BMR simulation

Sex = 'Male';
Age = 'Adult';

% load WBM
modelPath = ['../models/WBMs/', Age, Sex, '_WBM_withoutPSCM.mat'];
load(modelPath);
WBM = newModel1;

% load BMR data
DataPath = ['../data/', Age,  Sex, '_BMR.tsv'];
BMRData = readtable(DataPath, 'Delimiter', '\t', 'FileType', 'text');
Height = BMRData.Height_cm_;
Weight = BMRData.Fat_freeBodyMass_kg_;
% Weight = BMRData.Weight_kg_;

QPSolver = 'ILOGcomplex';

sex =  lower(Sex);
gender = lower(Sex);

Results_BMR = zeros(length(Height), 1);
for j=1:length(Height)
    % disp(j);
    H = Height(j,1);
    W = Weight(j,1);
    WBMc = WBM;

    standardPhysiolDefaultParameters;
    
    IndividualParameters_ = IndividualParameters;
    % set PSCM to WBM. This function comes from CobraTool.
    % There would get some warning about reactions not in the model, 
    % including 'Gall_H2Ot[bdG]','Brain_DM_atp_c_' and 'Heart_DM_atp_c_'. 
    % Gall is not included by this stduy. As for 'Brain_DM_atp_c_' and 
    % 'Heart_DM_atp_c_' is no need to demand in this research. 
    % Thus, we ignore them.
    IndividualParameters_.bodyWeight = W;
    IndividualParameters_.Height = H;
    WBMc = physiologicalConstraintsHMDBbased(WBMc,IndividualParameters_);

    % set Diet constraints to WBM. This function comes from CobraTool
    EUAverageDietNew
    WBMc = setDietConstraints(WBMc, Diet);

    % set some more constraints
    WBMc = setSimulationConstraints(WBMc);
    WBMc.lb(strmatch('BBB_KYNATE[CSF]upt',WBMc.rxns)) = -1000000;
    WBMc.lb(strmatch('BBB_LKYNR[CSF]upt',WBMc.rxns)) = -1000000;
    WBMc.lb(strmatch('BBB_TRP_L[CSF]upt',WBMc.rxns)) = -1000000;
    WBMc.ub(strmatch('Brain_EX_glc_D(',WBMc.rxns)) = -100; 
    
    %WBMc = addReaction(WBMc, 'Brain_ATPM', 'Brain_MAM01371c + Brain_MAM02040c -> Brain_MAM01285c + Brain_MAM02039c + Brain_MAM02751c');
    %WBMc = changeRxnBounds(WBMc, 'Brain_ATPM', 1000, 'l');

    %WBMc = addReaction(WBMc, 'Liver_ATPM', 'Liver_MAM01371c + Liver_MAM02040c -> Liver_MAM01285c + Liver_MAM02039c + Liver_MAM02751c');
    %WBMc = changeRxnBounds(WBMc, 'Liver_ATPM', 200, 'l');
    
    % WBMc = releaseWBMConstraints_old(WBMc, Sex, Age);
    WBMc = releaseWBMConstraints(WBMc, Sex, Age);
    
    % WBMc = changeRxnBounds(WBMc, 'Whole_body_objective_rxn', 1, 'b');
    WBMc = changeObjective(WBMc, 'Whole_body_objective_rxn');
    WBMc = changeRxnBounds(WBMc, 'Whole_body_objective_rxn', 1, 'b');
    [GeneClasses, RxnClasses, modelIrrevFM, MinimizedFlux] = pFBA(WBMc, 'geneoption',0, 'tol',1e-6, 'skipclass', 1);
    BMR_value = HumanBMR(modelIrrevFM, MinimizedFlux);
    
%     WBMc = changeObjective(WBMc,'Whole_body_objective_rxn');
%     [solution_male] = computeMin2Norm_HH(WBMc, QPSolver);
%     solution_male.x = solution_male.v;
%     BMR_value = HumanBMR(WBMc, solution_male);

    disp(BMR_value(1,1))
    
    Results_BMR(j,1) = BMR_value(1,1);
end

% output results
BMRdata = addvars(BMRData, Results_BMR);
BMRdata.Properties.VariableNames{end} = 'PredictedBMR';
outputFileName = ['../results/BMR/Reslut_' Age Sex '_BMR.tsv'];
writetable(BMRdata, outputFileName, 'Delimiter', '\t', 'FileType', 'text');


