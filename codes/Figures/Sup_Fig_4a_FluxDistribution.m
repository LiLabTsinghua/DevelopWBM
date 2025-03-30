% load FVA max/min flux results
df_Human1 = readtable('../../results/FVA/FVA_comp_Human1.txt');
df_Human2 = readtable('../../results/FVA/FVA_comp_Human2.txt');

rxns_Human1 = df_Human1.rxns;
flux_Human1 = df_Human1.model_ranges;
flux_ecHuman1 = df_Human1.ecModel_ranges;

rxns_Human2 = df_Human2.rxns;
flux_Human2 = df_Human2.model_ranges;
flux_ecHuman2 = df_Human2.ecModel_ranges;

% get same reactions
df_rxns = readtable('../../data/reactions.tsv', 'Delimiter', '\t', 'FileType', 'text');
All_MARs = df_rxns.rxns;
All_HMRs = df_rxns.rxnRetired;

inter_Rxns = {};
inter_flux_Human1 = {};
inter_flux_ecHuman1 = {};
inter_flux_Human2 = {};
inter_flux_ecHuman2 = {};
for i=1:length(All_MARs)
    r_1 = All_MARs{i,1};
    r_2 = All_HMRs{i,1};
    if ~isempty(r_1) && ~isempty(r_2)
        if any(strcmp(r_1, rxns_Human2)) && any(strcmp(r_2, rxns_Human1))
            idx1 = findIndex(rxns_Human2, r_1);
            idx2 = findIndex(rxns_Human1, r_2);

            inter_Rxns{end+1,1} = r_1;
            inter_flux_Human1{end+1,1} = flux_Human1(idx2);
            inter_flux_ecHuman1{end+1,1} = flux_ecHuman1(idx2);
            inter_flux_Human2{end+1,1} = flux_Human2(idx1);
            inter_flux_ecHuman2{end+1,1} = flux_ecHuman2(idx1);
        end
    end
end

% parameter setting
Human_FVA_Dists = {inter_flux_Human1; inter_flux_ecHuman1; inter_flux_Human2; inter_flux_ecHuman2};
legends = {'Human1', 'ecHuman1', 'Human2', 'ecHuman2', 'Box', 'off'};
titleStr   = 'Flux variability cumulative distribution';
colors = {'#d7191c', '#fdae61', '#2c7bb6', '#abd9e9'};
%[~, stats] = plotCumDist_v1(Human_FVA_Dists,legends,titleStr,colors);
[y_param, stats] = plotCumDist_v2(Human_FVA_Dists, legends, colors);
%[~, stats] = plotCumDist_v2(Human_FVA_Dists,legends,titleStr,colors);

