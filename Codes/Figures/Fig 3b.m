% load FVA max/min flux results
load('FVA1_pool20.7.mat')
load('FVA_Human1_pool20.7(1).mat')

[~, rxn_id] = xlsread('D:\All_Human_GTEx\EC_models\Human1_2_inter_flux_distribution.xlsx', 'Sheet1', 'A2:A9469');
[Human1, ~] = xlsread('D:\All_Human_GTEx\EC_models\Human1_2_inter_flux_distribution.xlsx', 'Sheet1', 'B2:B9469');
[ecHuman1, ~] = xlsread('D:\All_Human_GTEx\EC_models\Human1_2_inter_flux_distribution.xlsx', 'Sheet1', 'C2:C9469');
[Human2, ~] = xlsread('D:\All_Human_GTEx\EC_models\Human1_2_inter_flux_distribution.xlsx', 'Sheet1', 'D2:D9469');
[ecHuman2, ~] = xlsread('D:\All_Human_GTEx\EC_models\Human1_2_inter_flux_distribution.xlsx', 'Sheet1', 'E2:E9469');

% parameter setting
% Human_FVA_Dists = [Human1_FVA_Dists, Human2_FVA_Dists];
Human_FVA_Dists = {Human1, ecHuman1, Human2, ecHuman2};
legends = {'Human1-model', 'Human1-ecModel', 'Human2-model', 'Human2-ecModel', 'Box', 'off'};
% titleStr   = 'Flux variability cumulative distribution';
colors = {'#2c7bb6', '#abd9e9', '#d7191c', '#fdae61'};
[~, stats] = plotCumDist_v1(Human_FVA_Dists,legends,titleStr,colors);

