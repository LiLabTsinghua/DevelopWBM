% Time cost of FBA calculation

Sex = 'Male';
Age = 'Adult';

% load WBM
modelPath = ['../Models/', Age, '_', Sex, '_WBM.mat'];
load(modelPath);

WBM = releaseWBMConstraints(WBM, Sex, Age);

All_time = {};
for i=1:10
    t_start_1 = tic;
    WBMc = WBM;
    WBMc = changeObjective(WBMc,'Whole_body_objective_rxn');
    sol = optimizeCbModel(WBMc, 'max', 'one');
    disp(sol.f)
    elapsed_time_1 = toc(t_start_1);
    
    All_time{end+1} = elapsed_time_1;
end



