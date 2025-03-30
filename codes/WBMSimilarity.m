% Use the Jaccard Index to compare the WBM similarities among Fetus, 
% Adult Male, Elderly Male, Adult Female, and Elderly Female.
WBM_names = {'Fetus';'AdultMale';'AdultFemale';'ElderlyMale';'ElderlyFemale'};
for i=1:length(WBM_names)
    name = WBM_names{i,1};
    model = load(['../models/WBMs/',name, '_WBM_withoutPSCM.mat']);
    WBMs{i,1} = model.newModel1;
end

includeGene = false;
Similarity_values = {};
for i=1:length(WBMs)
    model1 = WBMs{i,1};
    for j=1:length(WBMs)
        model2 = WBMs{j,1};
        [JI, ~, ~] = JaccardIndex(model1,model2,includeGene);
        Similarity_values{i,j} = JI;
    end
end


tbl = cell2table(Similarity_values, 'VariableNames', WBM_names);
tbl.Properties.RowNames = WBM_names;
fileNames = '../results/WBM_Similarity.tsv';
writetable(tbl, fileNames, 'Delimiter', '\t', 'FileType', 'text', 'WriteRowNames', true);