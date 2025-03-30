%% 

Age = 'Adult';
Sex = 'Male';

% load WBM
modelPath = ['../models/WBMs/', Age, Sex, '_WBM_withoutPSCM.mat'];
load(modelPath);

WBM = releaseWBMConstraints(WBM, Sex, Age);
geneList = WBM.genes;

% geneList = WBM.genes;
All_solutionDel = {};
All_pfba_energy = {};

environment = getEnvironment();
parfor i=1:length(geneList)
    restoreEnvironment(environment);
    WBMc = WBM;
    WBMc = changeObjective(WBMc, 'Whole_body_objective_rxn');
    WBMc = changeRxnBounds(WBMc, 'Whole_body_objective_rxn', 1, 'b');

    del_gene = geneList{i,1};
    [modelDel, ~, constrRxnNames] = deleteHumanGenes(WBMc, {del_gene});
    
    try
        [GeneClasses, RxnClasses, modelIrrevFM, MinimizedFlux] = pFBA(modelDel, 'geneoption',0, 'tol',1e-6, 'skipclass', 1);
        BMR_value = HumanBMR(modelIrrevFM, MinimizedFlux);
        All_pfba_energy{i,1} = BMR_value(1,1);
    catch
        All_pfba_energy{i,1} = 'NA';
    end
    disp(i)
end

T = table(geneList(:), All_pfba_energy(:), 'VariableNames', {'DelGene', 'BMR'});
writetable(T, ['../results/', Age, Sex, '_delGene_BMR.tsv'], 'FileType', 'text', 'Delimiter', '\t');

