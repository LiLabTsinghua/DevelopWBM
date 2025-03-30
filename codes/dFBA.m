load('../models/AdultMale_coreWBM.mat')
model = coreWBM;

% Load metabolites related to the centralmetabolism, HMDB and exchange
removeTask = 'KeepMets.tsv'; % Some metabolic tasks should not be included in WBM.
KeepMets = readtable(['../data/', removeTask], 'Delimiter', '\t', 'FileType', 'text');
KeepMets_id = KeepMets.ID;


% Find reactions that all of these metabolites involve in the metabolite listed above
[metList, ~] = findMetsFromRxns(model, model.rxns); % a function from Cobra

KeepRxns = zeros(length(model.rxns), 1);
identifiers = cell(length(metList), 1);
for i=1:length(metList)
    mlist = metList{i,1};
    identifiers = cell(length(mlist), 1);
    for j=1:length(mlist)
        if contains(mlist{j,1}, '_')
            underscorePos = strfind(mlist{j,1}, '_');
            identifiers{j,1} = mlist{j,1}(underscorePos(end)+1:end-1);
        else
            identifiers{j,1} = mlist{j,1}(1:end-1);
        end
    end
    allInB = all(cellfun(@(x) any(strcmp(x, KeepMets_id)), identifiers));
    result = double(allInB);
    if result
        KeepRxns(i,1) = 1;
    end
end

KeepRxns_Idx = find(KeepRxns == 1);
% keep the NGAM reactions
NGAM = {'Brain_MAR03964'; 'Heart_MAR03964'; 'Liver_MAR03964';...
        'Lung_MAR03964'; 'Kidney_MAR03964'; 'Muscle_MAR03964';...
        'Skin_MAR03964'; 'WBM-ATPM'};
KeepRxns_id = [model.rxns(KeepRxns_Idx); NGAM];
NoKeepRxns = setdiff(model.rxns, KeepRxns_id);

% remove unrelated reactions
reducedModel=removeReactions(model,NoKeepRxns,true,true,true);
reducedModel = rmfield(reducedModel, 'rules');

% In order to achieve ketone exchange, some transport reactions need 
% to adjust the boundary.
KetoneTrans = {'Liver_MAR05018'; 'Brain_MAR04976';...
               'Brain_MAR04977'; 'Liver_MAR00174';...
               'Kidney_MAR00174'};
reducedModel = changeRxnBounds(reducedModel, KetoneTrans, -1000, 'l');
reducedModel = changeRxnBounds(reducedModel, KetoneTrans, 1000, 'u');

% uniqueToCellB = setdiff(model_before.rxns, KeepRxns_id);

%% Generate ecModel

% remove prefix of metabolite name before running GECKO3
MetNames = reducedModel.metNames;
MetNames_prefix = {'Brain_'; 'Heart_'; 'Liver_'; 'Lung_'; 'Muscle_'; 'Kidney_'; 'Skin_'};
for i = 1:length(MetNames_prefix)
    pfix = MetNames_prefix{i,1};
    for j=1:length(MetNames)
        name = MetNames{j,1};
        if contains(name, pfix)
            MetNames{j} = strrep(name, pfix, '');
        end
    end
end
reducedModel.metNames = MetNames;

save('../models/reduced_coreModel_0325.mat', 'reducedModel')

% Please refer to GECKO3.0 to generate ecModel.
% https://github.com/SysBioChalmers/GECKO/tree/main/tutorials/light_ecModel
load('../models/ecreduced_coreModel.mat')
model = ecModel;
model = changeObjective(model, 'WBM-ATPM');
res_WBM = optimizeCbModel(model, 'max', 'one');



%% Set constraints for fasting simulation

% Exchange rxns
[selExc, ~] = findExcRxns(model);
ExcRxns = model.rxns(selExc);
ExcRxns(end) = []; % Not include "prot_pool_exchange"

% close uptake rxns bounds
model = changeRxnBounds(model, ExcRxns, 0, 'l');

% ExcLb = model.lb(selExc);
% model.lb(selExc) = 0;

% open H2O, O2, Glycogen and Palmitate
openRxns = {'EX_MAM02630e'; 'EX_MAM01596e'; 'EX_MAM02040e'};
model = changeRxnBounds(model, openRxns, -1000, 'l');

% open substrates uptake rxns during fasting
substrates_rxns = {'EX_Liver_MAM03161e'; 'EX_MAM02674e'};
model = changeRxnBounds(model, substrates_rxns, -1000, 'l');

% Lactate transport
Lac_e = {'Heart_Exchange_Lac'; 'Lung_Exchange_Lac';...
         'Brain_Exchange_Lac'; 'Liver_Exchange_Lac';...
         'Kidney_Exchange_Lac'; 'Skin_Exchange_Lac';...
         'Muscle_Exchange_Lac'};
model = changeRxnBounds(model, Lac_e, -1000, 'l');
model = changeRxnBounds(model, Lac_e, 1000, 'u');
Lac_index = findIndex(model.rxns, Lac_e);

model = changeRxnBounds(model, 'Muscle_Exchange_Lac', 0, 'l');
model = changeRxnBounds(model, 'Muscle_Exchange_Lac', 1000, 'u');

% model = changeRxnBounds(model, 'Brain_Exchange_Lac', 0, 'l');
% model = changeRxnBounds(model, 'Brain_Exchange_Lac', 1000, 'u');

% To facilitate visualization, reversed uptake rxns are turned off
Lac_e_rev = {'Heart_Exchange_Lac_REV'; 'Lung_Exchange_Lac_REV';...
             'Brain_Exchange_Lac_REV'; 'Liver_Exchange_Lac_REV';...
             'Kidney_Exchange_Lac_REV'; 'Skin_Exchange_Lac_REV'};
model = changeRxnBounds(model, Lac_e_rev, 0, 'b');
Lac__rev_index = findIndex(model.rxns, Lac_e_rev);

% Glucose transport
Glu_e = {'Liver_Trans_MAM01965e'; 'LungGlucose_Trans';...
         'Brain_Trans_MAM01965e'; 'HeartGlucose_Trans';...
         'Muscle_Trans_MAM01965e'; 'SkinGlucose_Trans';...
         'Muscle_Trans_MAM01965e'; 'Kidney_Trans_MAM01965e'};
Glu_index = findIndex(model.rxns, Glu_e);
model = changeRxnBounds(model, Glu_e, -1000, 'l');

model = changeRxnBounds(model, 'Muscle_Trans_MAM01965e', -1000, 'l');
model = changeRxnBounds(model, 'Muscle_Trans_MAM01965e', 0, 'u');

%model = changeRxnBounds(model, 'Brain_Trans_MAM01965e', -1000, 'l');
%model = changeRxnBounds(model, 'Brain_Trans_MAM01965e', 1000, 'u');

% To facilitate visualization, reversed uptake rxns are turned off
Glu_e_rev = {'Liver_MAR05027_REV'; 'LungGlucose_Trans_REV';...
             'HeartGlucose_Trans_REV'; 'SkinGlucose_Trans_REV';...
             'MAR09034'; 'Diet_EX_glc_D[d]'; 'DM_MAM01253e'}; % The normal body does not allow glucose to be excreted
Glu__rev_index = findIndex(model.rxns, Glu_e_rev);
model = changeRxnBounds(model, Glu_e_rev, 0, 'b');


% ketone
ketone_e = {'LiverBHB_Trans'; 'LiverACAC_Trans';...
         'LiverACET_Trans'; 'BrainBHB_Trans';...
         'BrainBHBACAC_Trans'};
ketone_index = findIndex(model.rxns, ketone_e);
model = changeRxnBounds(model, ketone_e, -1000, 'l');

ketone_e_rev = {'LiverBHB_Trans_REV'; 'LiverACAC_Trans_REV';...
             'LiverACET_Trans_REV'; 'BrainBHB_Trans_REV';...
             'BrainBHBACAC_Trans_REV'};
ketone__rev_index = findIndex(model.rxns, ketone_e_rev);
model = changeRxnBounds(model, ketone_e_rev, 0, 'b');


% close some duplicated transport rxns
Glu_duplicated_rxns = {'Kidney_Trans_MAM01965e'; 'Brain_Trans_MAM01965e';...
                        'Liver_Trans_MAM01965e'; 'Muscle_Trans_MAM01965e'};


model = changeRxnBounds(model, 'Liver_MAR06033_REV', 0, 'u');
model = changeRxnBounds(model, 'Liver_MAR06001_REV', 0, 'u');
model = changeRxnBounds(model, 'DM_MAM01253e', 0, 'u');
model = changeRxnBounds(model, 'DM_MAM00157e', 0, 'u');

model = changeRxnBounds(model, 'LiverACAC_Trans', 0, 'u');
model = changeRxnBounds(model, 'MAR09132', 1000, 'u');
model = changeRxnBounds(model, 'DM_MAM00157e', 0, 'u');

model = changeRxnBounds(model, 'BrainBHBACAC_Trans', 0, 'u');

model = changeRxnBounds(model, 'MAR09086', 0, 'b');
model = changeRxnBounds(model, 'Diet_EX_ac[d]', 0, 'l');
model = changeRxnBounds(model, 'Diet_EX_ac[d]', 1000, 'u');

model = changeRxnBounds(model, 'LungHe_trans', 0, 'b');
model = changeRxnBounds(model, 'HeartHe_trans', 0, 'b');


% 
% model = changeRxnBounds(model, 'BrainHe_trans_REV', 0, 'b');
% 
% model = changeRxnBounds(model, 'LiverHe_trans_REV', 0, 'b');
% 
% model = changeRxnBounds(model, 'LungHe_trans_REV', 0, 'b');
% model = changeRxnBounds(model, 'HeartHe_trans_REV', 0, 'b');
% model = changeRxnBounds(model, 'MuscleHe_trans_REV', 0, 'b');
% model = changeRxnBounds(model, 'KidneyHe_trans_REV', 0, 'b');
% model = changeRxnBounds(model, 'SkinHe_trans_REV', 0, 'b');

% model = changeRxnBounds(model, 'Brain_Trans_MAM01252e', 0, 'b');



% ¸Ä±äAKG->succoa Hm
% change_rxns = {'Liver_MAR06409', 'Lung_MAR06409', 'Heart_MAR06409', 'Muscle_MAR06409', 'Brain_MAR06409', 'Skin_MAR06409', 'Kidney_MAR06409'};
% met_id = {'Liver_MAM02039m', 'Lung_MAM02039m', 'Heart_MAM02039m', 'Muscle_MAM02039m', 'Brain_MAM02039m', 'Skin_MAM02039m', 'Kidney_MAM02039m'};
% 
% for i = 1:length(change_rxns)
%     r = change_rxns{1,i};
%     m = met_id{1,i};
%     r_index = findIndex(model.rxns, r);
%     m_index = findIndex(model.mets, m);
%     model.S(m_index,r_index) = -2;
% end

%%
% [model, ~] = addReaction(model, 'Liver_MAR05018_REV', 'Liver_MAM02039c + Liver_MAM00157c + 6.00594 prot_pool ->	Liver_MAM02039e + Liver_MAM00157e');
% [model, ~] = addReaction(model, 'Brain_MAR04976_REV', 'Brain_MAM02039c + Brain_MAM01253c + 1.10737 prot_pool ->	Brain_MAM02039e + Brain_MAM01253e');
% [model, ~] = addReaction(model, 'Brain_MAR04977_REV', 'Brain_MAM02039m + Brain_MAM01253m + 4.34244 prot_pool ->	Brain_MAM02039c + Brain_MAM01253c');
% % [model, ~] = addReaction(model, 'Liver_MAR00174', 'Liver_MAM01371c + Liver_MAM01597c + Liver_MAM02642c + 0.294871 prot_pool ->	Liver_MAM01334c + Liver_MAM02759c + Liver_MAM02644c');
% [model, ~] = addReaction(model, 'Liver_MAR00174_REV', 'Liver_MAM01334c + Liver_MAM02759c + Liver_MAM02644c + 0.294871 prot_pool ->	Liver_MAM01371c + Liver_MAM01597c + Liver_MAM02642c');
% 
% % [model, ~] = addReaction(model, 'Kidney_MAR00174', 'Kidney_MAM01371c + Kidney_MAM01597c + Kidney_MAM02642c + 0.294871 prot_pool ->	Kidney_MAM01334c + Kidney_MAM02759c + Kidney_MAM02644c');
% [model, ~] = addReaction(model, 'Kidney_MAR00174_REV', 'Kidney_MAM01334c + Kidney_MAM02759c + Kidney_MAM02644c + 0.294871 prot_pool  -> Kidney_MAM01371c + Kidney_MAM01597c + Kidney_MAM02642c');


%% change abnormal prot_mass

for i=1:length(model.rxns)
    r = model.rxns{i,1};
    m = 'prot_pool';
    r_index = findIndex(model.rxns, r);
    m_index = findIndex(model.mets, m);
    if abs(model.S(m_index,r_index)) > 50
        model.S(m_index,r_index) = model.S(m_index,r_index)*0.001;
    end
end

%% Exchange Index
[selExc, selUpt] = findExcRxns(model);
index1 = find(selExc);
% index2 = find(selUpt);
% all_ex_index = [index1; index2];
all_ex_index = index1;

% find transport rxns
[metList, ~] = findMetsFromRxns(model, model.rxns); 
trans_rxn_check = zeros(length(model.rxns), 1);
for i=1:length(metList)
    RxnMets = metList{i,1};
    if numel(RxnMets) == 2
        MetComps = model.metComps(findIndex(model.mets, RxnMets));
        if all(MetComps == 1, 'all')
            trans_rxn_check(i,1) = 1;
        end
    end
end

trans_rxn_index = find(trans_rxn_check == 1);
trans_rxn_id = model.rxns(trans_rxn_index);


%% fasting ATPM simulation

% Determination of glycogen parameter
Sg=zeros(100,1); % Glycogen concentration
Sg(1)=100; % initial Glycogen concentration for the first step (mmol/gDCW)

v1=zeros(102,1); % Flux of glycogen degradation
% This Km value refers to the glycogen phosphorylase (EC: 2.4.1.1) data
% in the Brenda database. The Km value of glycogen phosphorylase when 
% glycogen is the substrate is between 0.0138-40 (mM), and we choose 
% 20 mM as the reference.
% vmaxg = 20; 
vmaxg = 25; 
Kmg = 20;


% Determination of glycogen parameter
St=zeros(100,1); % Glycogen concentration
St(1)=100; % initial FA concentration for the first step (mmol/gDCW)

% Multiple tissues can uptake fat. In order to simplify, enzyme parameters
% related to fatty acid activation (6.2.1.3) were used to constraints the 
% corresponding exchange reation of fat.
v2=zeros(102,1); % Flux of Fat uptake
% This Km value refers to the long-chain-fatty-acid-CoA ligase (EC: 6.2.1.3)
% data in the Brenda database. The Km value of the ligase is 0.00021 - 6(mM)
% with palmitate as the substrate, and we choose 3 mM as the reference.
vmaxt = 25; % same as glycogen degradation
Kmt=3;


% Lactate
Slg=zeros(102,1);
Slg(1)=0; % initial lactate concentration for the first step (mmol/L)
vlg=zeros(102,1);

% O2
SOg=zeros(102,1);
SOg(1)=0; % initial O2 concentration for the first step (mmol/L)
vOg=zeros(102,1);

% CO2
SCg=zeros(102,1);
SCg(1)=0; % initial CO2 concentration for the first step (mmol/L)
vCg=zeros(102,1);

% LiverACAC_Trans
SLiver_acacg=zeros(102,1);
SLiver_acacg(1)=0; 
vACACg=zeros(102,1);

% LiverBHB_Trans
SLiver_BHBg=zeros(102,1);
SLiver_BHBg(1)=0;
vBHBg=zeros(102,1);

v_all_exchange=zeros(170,102);
v_all_trans=zeros(320,102);

miu=zeros(102,1);
All_ATPM=zeros(102,1);

SLacg=zeros(102,6);
vLacg=zeros(7,102);

SGlu=zeros(102,6);
vGlus=zeros(7,102);

% FBAg=zeros(9526,102);
FBAg=zeros(length(model.rxns),102);

model = changeRxnBounds(model, 'prot_pool_exchange', -1000, 'l');
model = changeRxnBounds(model, 'EX_MAM02674e', -1000, 'l');
model = changeRxnBounds(model, 'EX_MAM02674e', 1000, 'u');
model = changeRxnBounds(model, 'EX_Liver_MAM03161e', -1000, 'l');
model = changeRxnBounds(model, 'EX_Liver_MAM03161e', 1000, 'u');

model = changeObjective(model, 'WBM-ATPM');
res_WBM = optimizeCbModel(model, 'max', 'one');
lowerBound = floor(res_WBM.f);
upperBound = ceil(res_WBM.f);

for j=1:102
    model_ALE(j) = model;
    
    if j==1
        % model_ALE(j) = changeObjective(model_ALE(j), 'WBM_ATPM');
        % res_WBM = optimizeCbModel(model_ALE(j), 'max', 'one');
        % v_ATP = res_WBM.x(findIndex(model_ALE(j).rxns, 'WBM_ATPM'));
        % v1(j) = res_WBM.x(findIndex(model.rxns, 'EX_MAM03161e'));
        % v2(j) = res_WBM.x(findIndex(model.rxns, 'EX_MAM02674e'));
    
        model_ALE(j) = changeRxnBounds(model_ALE(j),'WBM-ATPM', lowerBound, 'l');
        model_ALE(j) = changeRxnBounds(model_ALE(j),'WBM-ATPM', upperBound, 'u');
        model_ALE(j) = changeRxnBounds(model_ALE(j),'EX_Liver_MAM03161e', -1000,'l');
        model_ALE(j) = changeRxnBounds(model_ALE(j),'EX_Liver_MAM03161e', 1000,'u');
        model_ALE(j) = changeRxnBounds(model_ALE(j),'EX_MAM02674e',-1000,'l');
    	model_ALE(j) = changeRxnBounds(model_ALE(j),'EX_MAM02674e', 1000,'u');
        
        
        model_ALE(j) = changeObjective(model_ALE(j),'prot_pool_exchange');
        FBAsolution = optimizeCbModel(model_ALE(j),'max','one');
        
        % vlg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM02403e')); % Lactate
        vCg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM01596e')); % CO2
        vOg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM02630e')); % O2
        vACACg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'BrainBHBACAC_Trans')); 
        vBHBg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'LiverBHB_Trans'));
        
        v_all_exchange(:,j)=FBAsolution.x(all_ex_index);
        v_all_trans(:,j)=FBAsolution.x(trans_rxn_index);
        
        vLacg(:,j) = FBAsolution.x(Lac_index);
        vLacg(:,j) = FBAsolution.x(Lac_index);
        
        vGlus(:,j) = FBAsolution.x(Glu_index);
        
        miu(j)=FBAsolution.x(findIndex(model.rxns, 'prot_pool_exchange'));
        All_ATPM(j)=FBAsolution.x(findIndex(model.rxns, 'WBM-ATPM'));
        
        v1(j) = FBAsolution.x(findIndex(model.rxns, 'EX_Liver_MAM03161e'));
        v2(j) = FBAsolution.x(findIndex(model.rxns, 'EX_MAM02674e'));
        
        Slg(j+1)=Slg(j)+vlg(j)*0.1;% set Lactate concentration for the next step (mmol/gDCW)
        SCg(j+1)=SCg(j)+vCg(j)*0.1;% set Lactate concentration for the next step (mmol/gDCW)
        SOg(j+1)=SOg(j)-vOg(j)*0.1;% set O2 concentration for the next step (mmol/gDCW)
        SLiver_acacg(j+1)=SLiver_acacg(j)+vACACg(j)*0.1;
        SLiver_BHBg(j+1)=SLiver_BHBg(j)+vBHBg(j)*0.1;
        
    
        Sg(j+1)=Sg(j)+v1(j)*0.1;% set Glycogen concentration for the next step (mmol/gDCW)
        St(j+1)=St(j)+v2(j)*0.1;% set FA concentration for the next step (mmol/gDCW)
        
        v1(j+1) = vmaxg*Sg(j+1)/(Kmg+Sg(j+1));
        v2(j+1) = 1000;
        
        FBAg(:,j)=FBAsolution.x; 
        
    else
        if Sg(j) >= 0.0001
            % if abs(St(j)-St(j-1)) <= 0.001
                model_ALE(j) = changeRxnBounds(model_ALE(j), 'EX_Liver_MAM03161e', -v1(j), 'b');
                model_ALE(j) = changeRxnBounds(model_ALE(j), 'EX_MAM02674e', -v2(j), 'l');
                
                model_ALE(j) = changeObjective(model_ALE(j), 'WBM-ATPM');
                res_WBM = optimizeCbModel(model_ALE(j), 'max', 'one');
                atp = res_WBM.x(findIndex(model_ALE(j).rxns, 'WBM-ATPM'));
            
                model_ALE(j) = changeRxnBounds(model_ALE(j), 'WBM-ATPM', atp, 'l');
            
                model_ALE(j) = changeObjective(model_ALE(j), 'prot_pool_exchange');
                FBAsolution = optimizeCbModel(model_ALE(j), 'max', 'one');
                
                % vlg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM02403e')); % Lactate
                vCg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM01596e')); % CO2
                vOg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM02630e')); % O2
                vACACg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'BrainBHBACAC_Trans'));
                vBHBg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'LiverBHB_Trans'));
                
                v_all_exchange(:,j)=FBAsolution.x(all_ex_index);
                v_all_trans(:,j)=FBAsolution.x(trans_rxn_index);
                
                vLacg(:,j) = FBAsolution.x(Lac_index);
                
                vGlus(:,j) = FBAsolution.x(Glu_index);
                
                miu(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'prot_pool_exchange'));
                All_ATPM(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'WBM-ATPM'));
        
                v1(j) = FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_Liver_MAM03161e'));
                v2(j) = FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM02674e'));
        
                Slg(j+1)=Slg(j)+vlg(j)*0.1;% set Lactate concentration for the next step (mmol/gDCW)
                SCg(j+1)=SCg(j)+vCg(j)*0.1;% set Lactate concentration for the next step (mmol/gDCW)
    
                Sg(j+1)=Sg(j)+v1(j)*0.1;% set Glycogen concentration for the next step (mmol/gDCW)
                St(j+1)=St(j)+v2(j)*0.1;% set FA concentration for the next step (mmol/gDCW)
                SOg(j+1)=SOg(j)-vOg(j)*0.1;% set O2 concentration for the next step (mmol/gDCW)
                SLiver_acacg(j+1)=SLiver_acacg(j)+vACACg(j)*0.1;
                SLiver_BHBg(j+1)=SLiver_BHBg(j)+vBHBg(j)*0.1;
        
                v1(j+1) = vmaxg*Sg(j+1)/(Kmg+Sg(j+1));
                v2(j+1) = vmaxt*St(j+1)/(Kmt+St(j+1));
                FBAg(:,j)=FBAsolution.x;
                
        else
            model_ALE(j) = changeRxnBounds(model_ALE(j), 'EX_Liver_MAM03161e', 0, 'l');
            model_ALE(j) = changeRxnBounds(model_ALE(j), 'EX_Liver_MAM03161e', 1000, 'u');
%             model_ALE(j) = changeRxnBounds(model_ALE(j), 'EX_MAM02674e', -v2(j), 'l');
            model_ALE(j) = changeRxnBounds(model_ALE(j), 'EX_MAM02674e', -v2(j), 'b');
            
            model_ALE(j) = changeObjective(model_ALE(j), 'WBM-ATPM');
            res_WBM = optimizeCbModel(model_ALE(j), 'max', 'one');
            atp = res_WBM.x(findIndex(model_ALE(j).rxns, 'WBM-ATPM'));
            
            model_ALE(j) = changeRxnBounds(model_ALE(j), 'WBM-ATPM', atp, 'l');
            
            model_ALE(j) = changeObjective(model_ALE(j), 'prot_pool_exchange');
            FBAsolution = optimizeCbModel(model_ALE(j), 'max', 'one');
            
            % vlg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM02403e')); % Lactate
            vCg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM01596e')); % CO2
            vOg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM02630e')); % O2
            vACACg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'BrainBHBACAC_Trans'));
            vBHBg(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'LiverBHB_Trans'));
            
            v_all_exchange(:,j)=FBAsolution.x(all_ex_index);
            v_all_trans(:,j)=FBAsolution.x(trans_rxn_index);
            
            vLacg(:,j) = FBAsolution.x(Lac_index);
            
            vGlus(:,j) = FBAsolution.x(Glu_index);
            
            miu(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'prot_pool_exchange'));
            All_ATPM(j)=FBAsolution.x(findIndex(model_ALE(j).rxns, 'WBM-ATPM'));
        
            v1(j) = FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_Liver_MAM03161e'));
            v2(j) = FBAsolution.x(findIndex(model_ALE(j).rxns, 'EX_MAM02674e'));
        
            Slg(j+1)=Slg(j)+vlg(j)*0.1;% set Lactate concentration for the next step (mmol/gDCW)
            SCg(j+1)=SCg(j)+vCg(j)*0.1;% set Lactate concentration for the next step (mmol/gDCW)
            SOg(j+1)=SOg(j)-vOg(j)*0.1;% set O2 concentration for the next step (mmol/gDCW)
            SLiver_acacg(j+1)=SLiver_acacg(j)+vACACg(j)*0.1;
            SLiver_BHBg(j+1)=SLiver_BHBg(j)+vBHBg(j)*0.1;
    
            Sg(j+1)=0;% set Glycogen concentration for the next step (mmol/gDCW)
            St(j+1)=St(j)+v2(j)*0.1;% set FA concentration for the next step (mmol/gDCW)
        
            v1(j+1) = 0;
            v2(j+1) = vmaxt*St(j+1)/(Kmt+St(j+1));
        
            FBAg(:,j)=FBAsolution.x; 
            
            
        end
    end
end

%% Plot
%% figure plot1 (Changes of glycogen and fatty acid concentration)
figure;
hold on;  % this is to hold the plotted stuff and keep them showing

time = 0:0.1:10.2;
% plot(time,SCg,'linewidth', 2, 'color', "#0072BD");
% plot(time,SCg,'linewidth', 2, 'color', [17 115 187]/256);
% plot(time,SOg,'linewidth', 2, 'color', [0 0.8 0]);
plot(time,Sg,'linewidth', 2, 'color', [17 115 187]/256);
plot(time,St,'linewidth', 2, 'color', [0.85 0.325 0.098]);
axis([0,12,0,150]);
legend({'Glycogen','Fatty acid'},'Location','northwest'); 
xlabel('Time(h-1)') 
ylabel('Concentration(mmol/gDCW)') 

%% figure plot 2 (The flux changes of lactate exchange reaction in each tissue)
b = vLacg;
b(:,[1]) = [];

figure;
hold on; 

time = 0:0.1:10.0;
plot(time, b(1,:), 'linewidth', 2, 'color', [17 115 187]/256);
plot(time, b(2,:),'linewidth', 2, 'color', [0 0.8 0]);
plot(time, b(3,:),'linewidth', 2, 'color', [0.6 0.5 0.3]);
plot(time, b(4,:),'linewidth', 2, 'color', [0.85 0.325 0.098]);
plot(time, b(5,:),'linewidth', 2, 'color', [0.635 0.078 0.184]);
plot(time, b(6,:),'linewidth', 2, 'color', [0.929 0.694 0.125]);
plot(time, b(7,:),'linewidth', 2, 'color', [0.229 0.594 0.325]);

legend(model.rxns(Lac_index));
xlabel('Time(h-1)') 
ylabel('Flux(mmol/gDCW/h)') 

%% figure plot 3 (The flux changes of glucose exchange reaction in each tissue)
c = vGlus;
c(:,[1]) = [];

figure;
hold on;

time = 0:0.1:10.0;
plot(time, c(1,:), 'linewidth', 2, 'color', [17 115 187]/256);
plot(time, c(2,:),'linewidth', 2, 'color', [0 0.8 0]);
plot(time, c(3,:),'linewidth', 2, 'color', [0.6 0.5 0.3]);
plot(time, c(4,:),'linewidth', 2, 'color', [0.85 0.325 0.098]);
plot(time, c(5,:),'linewidth', 2, 'color', [0.635 0.078 0.184]);
plot(time, c(6,:),'linewidth', 2, 'color', [0.929 0.694 0.125]);
plot(time, c(7,:),'linewidth', 2, 'color', [0.229 0.594 0.325]);

legend(model.rxns(Glu_index));
xlabel('Time(h-1)')
ylabel('Flux(mmol/gDCW/h)') 

%% figure plot 4
% pyr = {'LiverBHB_Trans'; 'LiverACAC_Trans';...
%        'LiverACET_Trans'; 'BrainBHB_Trans';...
%        'BrainBHBACAC_Trans'};
ketone = {'LiverBHB_Trans'; 'BrainBHB_Trans'};
   
% pyr = {'LiverBHB_Trans'; 'LiverACAC_Trans';...
%        'LiverACET_Trans'; 'BrainBHB_Trans';...
%        'BrainACAC_Trans'};
   
% pyr = {'LiverACAC_Trans'; 'LiverBHB_Trans'; 'BrainBHBACAC_Trans_REV'; 'BrainBHB_Trans'; 'LiverBHB_Trans_REV'; 'LiverACAC_Trans_REV'; 'LiverACET_Trans_REV'; 'BrainBHB_Trans_REV'};
ketone_index = findIndex(model.rxns, ketone);
p = FBAg(ketone_index,:);

p(:,[1]) = [];

figure;
hold on;

time = 0:0.1:10.0;
plot(time, p(1,:), 'linewidth', 2, 'color', [17 115 187]/256);
plot(time, p(2,:),'linewidth', 2, 'color', [0.929 0.694 0.125]);
% plot(time, p(3,:),'linewidth', 2, 'color', [0.6 0.5 0.3]);
% plot(time, p(4,:),'linewidth', 2, 'color', [0.85 0.325 0.098]);
% plot(time, p(5,:),'linewidth', 2, 'color', [0.635 0.078 0.184]);
% plot(time, p(6,:),'linewidth', 2, 'color', [0 0.8 0]);
% plot(time, p(7,:),'linewidth', 2, 'color', [0.229 0.594 0.325]);

% legend({'Pyr exchange in Brain', 'Pyr exchange in Heart', 'Pyr exchange in Kidney', 'Pyr exchange in Liver', 'Pyr exchange in Lung', 'Pyr exchange in Skin', 'Pyr exchange in Muscle'},'Location','northeast'); 
% legend({'Ketone(BHB) exchange in Liver'; 'Ketone(ACAC) exchange in Liver'; 'Ketone(ACET) exchange in Liver'; 'Ketone(BHB) exchange in Brain'; 'Ketone(ACAC) exchange in Brain'})
legend({'Ketone(BHB) exchange in Liver'; 'Ketone(BHB) exchange in Brain'})
xlabel('Time(h-1)')
ylabel('Flux(mmol/gDCW/h)') 

