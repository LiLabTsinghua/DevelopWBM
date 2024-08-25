% Run BMR simulation

Sex = 'Male';
Age = 'Adult';

% load WBM
modelPath = ['../Models/', Age, '_', Sex, '_WBM_withoutPSCM.mat'];
load(modelPath);

% load BMR data
DataPath = ['../Data/', Age,  Sex, '_BMR.tsv'];
BMRData = readtable(DataPath, 'Delimiter', '\t', 'FileType', 'text');
Height = BMRData.Height_cm_;
% Weight = BMRData.Fat_freeBodyMass_kg_;
Weight = BMRData.Weight_kg_;


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

    WBMc = releaseWBMConstraints(WBMc, Sex, Age);
    WBMc = changeObjective(WBMc, 'Whole_body_objective_rxn');
    [GeneClasses, RxnClasses, modelIrrevFM, MinimizedFlux] = pFBA(WBMc, 'geneoption',0, 'tol',1e-6, 'skipclass', 1);
    BMR_value = HumanBMR(modelIrrevFM, MinimizedFlux);
    Results_BMR(j,1) = BMR_value(1,1);
end



