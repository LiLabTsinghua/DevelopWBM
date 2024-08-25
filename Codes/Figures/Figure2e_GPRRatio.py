import math
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.io import load_yaml_model, save_yaml_model, load_matlab_model, save_matlab_model
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import stats
import seaborn as sns
import pandas as pd
from scipy.stats import ranksums

from proplot import rc
from pylab import *


def SubsystemRatio(model):
    all_sub = []
    for r in model.reactions:
        # for s in r.subsystem:
        if str(r.subsystem) not in all_sub:
            all_sub.append(str(r.subsystem))

    dict_sub_rxn_all = {key: 0 for key in all_sub}
    dict_sub_rxn_gpr = {key: 0 for key in all_sub}
    for i in all_sub:
        for j in model.reactions:
            if j.gene_reaction_rule != '':
                if i in j.subsystem:
                    dict_sub_rxn_all[i] += 1
                    dict_sub_rxn_gpr[i] += 1
            else:
                if i in j.subsystem:
                    dict_sub_rxn_all[i] += 1

    dict_sub_gpr_ratio = {}
    for k, v in dict_sub_rxn_gpr.items():
        dict_sub_gpr_ratio[k] = v/dict_sub_rxn_all[k]

    return dict_sub_rxn_all, dict_sub_gpr_ratio


def SubRatio_plot(df):
    # plt.style.use("seaborn-darkgrid")
    # sns.set_style("darkgrid")
    bubble_scale = 60 
    plt.rcParams.update({'font.size': 6})
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(1.8, 3))


    plt.scatter(
        df['Human1_rxn_num'],
        df['Subsystem'],
        s=df['Human1_sub_ratio'] * bubble_scale,
        color='#d7191c',
        # label = 'Human1',
        # alpha=0.7 
    )

    plt.scatter(
        df['Human2_rxn_num'],
        df['Subsystem'],
        s=df['Human2_sub_ratio'] * bubble_scale,
        color='#2c7bb6',
        # label = 'Human2',
        # alpha=0.7
    )

    # plt.rcParams['figure.facecolor'] = '#f0f0f0'

    size_values = [1, 1, 0.2, 0.6, 1]
    size_values_label = [element * bubble_scale for element in size_values]
    size_labels = ['Human2', 'Human1', '0.2', '0.6', '1']
    # size_labels = [str(element) for element in size_values]

    # set labels
    for label, size in zip(size_labels, size_values_label):
        if label == 'Human2':
            plt.scatter([], [], c='#2c7bb6', s=size, label=str(label))
        elif label == 'Human1':
            plt.scatter([], [], c='#d7191c', s=size, label=str(label))
        else:
            # plt.scatter([], [], c='k', alpha=0.3, s=size, label=str(label))
            plt.scatter([], [], c='#d6d2d6', s=size, edgecolors='black', linewidths=0.5, label=str(label))



    my_x_ticks = np.arange(0, 140, 20)
    plt.xticks(my_x_ticks, fontsize=6)


    plt.xlabel('Reactions')
    plt.grid(True) 

    # plt.colorbar(label='number of papers')

    plt.tick_params(direction='in') 
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    plt.legend(frameon=False, bbox_to_anchor=(1, 1))

    minorticks_off()

    ax.grid(False)
    ax = plt.gca()

    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    # plt.savefig("D:/All_Human_GTEx/Human2_Paper_and_Suplementary/Figures/Fig_2e_temp_2.pdf", dpi=400, bbox_inches = 'tight')

    plt.show()


Human2 = load_matlab_model('../../Models/'+'Human2.mat')
Human1 = load_matlab_model('../../Models/'+'Human1.mat')

Human2_SubRxns_all, Human2_SubGPR_rxns = SubsystemRatio(Human2)
Human1_SubRxns_all, Human1_SubGPR_rxns = SubsystemRatio(Human1)

Subs = ['Purine metabolism', 
        'Arginine and proline metabolism', 
        'Tyrosine metabolism', 
        'Pyrimidine metabolism', 
        'Acyl-CoA hydrolysis', 
        'Glycine, serine and threonine metabolism', 
        'Valine, leucine, and isoleucine metabolism', 
        'Cysteine and methionine metabolism', 
        'Glycolysis / Gluconeogenesis', 
        'Lysine metabolism', 
        'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism', 
        'Pyruvate metabolism', 
        'Pentose phosphate pathway', 
        'Protein assembly', 
        'Histidine metabolism', 
        'Tryptophan metabolism']

Human2_SubRxns = [Human2_SubRxns_all[i] for i in Subs]
Human2_SubRxns_ratio = [Human2_SubGPR_rxns[i] for i in Subs]
Human1_SubRxns = [Human1_SubRxns_all[i] for i in Subs]
Human1_SubRxns_ratio = [Human1_SubGPR_rxns[i] for i in Subs]

df = pd.DataFrame({'Subsystem':Subs,
                  'Human2_rxn_num':Human2_SubRxns,
                  'Human2_sub_ratio':Human2_SubRxns_ratio,
                  'Human1_rxn_num':Human1_SubRxns, 
                  'Human1_sub_ratio':Human1_SubRxns_ratio})

SubRatio_plot(df)