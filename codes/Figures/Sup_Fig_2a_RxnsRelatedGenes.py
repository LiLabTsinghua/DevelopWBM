import pandas as pd
import numpy as np
import re
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.io import load_yaml_model, save_yaml_model, load_matlab_model, save_matlab_model
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline


def Rxns_genes(model):
    dict_Nums = {'Rxn_1':0,
                 'Rxn_2':0,
                 'Rxn_3':0,
                 'Rxn_4':0,
                 'Rxn_5':0,
                 'Rxn_6':0,
                 'Rxn_7':0,
                 'Rxn_8':0,
                 'Rxn_9':0,
                 'Rxn_10':0}
    
    for g in model.genes:
            if len(g.reactions) > 0:
                if len(g.reactions) == 1:
                    dict_Nums['Rxn_1'] += 1
                elif len(g.reactions) == 2:
                    dict_Nums['Rxn_2'] += 1
                elif len(g.reactions) == 3:
                    dict_Nums['Rxn_3'] += 1
                elif len(g.reactions) == 4:
                    dict_Nums['Rxn_4'] += 1
                elif len(g.reactions) == 5:
                    dict_Nums['Rxn_5'] += 1
                elif len(g.reactions) == 6:
                    dict_Nums['Rxn_6'] += 1
                elif len(g.reactions) == 7:
                    dict_Nums['Rxn_7'] += 1
                elif len(g.reactions) == 8:
                    dict_Nums['Rxn_8'] += 1
                elif len(g.reactions) == 9:
                    dict_Nums['Rxn_9'] += 1
                elif len(g.reactions) >= 10:
                    dict_Nums['Rxn_10'] += 1
            
    return dict_Nums


def plot_bar(dict_Nums_human1, dict_Nums_human2):
    # Generate sample data
    # np.random.seed(42)
    data1 = [v for k, v in dict_Nums_human1.items()]
    data2 = [v for k, v in dict_Nums_human2.items()]

    # Create the figure and axis objects
    fig, ax1 = plt.subplots(figsize=(8, 3.5))

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype']= 42

    # Plot histograms
    bins = list(range(10))
    bins = [x+1 for x in bins]
    # bins[-1] = '>=10'
    ax1.bar(bins, data1, alpha=1, width = 1, label='Human1', color='#d7191c', edgecolor='black', linewidth=1.2)
    ax1.bar(bins, data2, alpha=1, width = 1, label='Human2', color='#2c7bb6', edgecolor='black', linewidth=1.2)
    # ax1.hist(data2, bins, alpha=0.7, label='Human2', density=True, color='#2c7bb6', edgecolor='black', linewidth=1.2)

    # Plot cumulative distribution functions
    x = np.linspace(0, 30, 200)
    x = [t+2 for t in x]

    indices = np.arange(len(bins))

    # Set labels and title
    ax1.set_xlabel('Rxns associated with a gene')
    ax1.set_ylabel('Gene counts')

    ax1.set_ylim(0, 1500)

    # plt.xticks([x+2 for x in bins])

    # Add legend
    ax1.legend(loc='upper right',frameon=False)

    for i in ['bottom', 'top', 'right', 'left']:
        ax1.spines[i].set_color('black')
        ax1.spines[i].set_visible('black')
        ax1.spines[i].set_linewidth(0.5)

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    # Adjust layout and display the plot
    plt.savefig("../../results/Figures/Fig_Gene_Rxns_nums.pdf", transparent = True, dpi=400, bbox_inches = 'tight')
    plt.tight_layout()
    plt.show()

def main():
    filepath = '../../models/'
    model_2 = cobra.io.load_matlab_model(filepath+'Human2.mat')
    model_1 = cobra.io.load_matlab_model(filepath+'Human1.mat')

    dict_Nums_human1 = Rxns_genes(model_1)
    dict_Nums_human2 = Rxns_genes(model_2)

    plot_bar(dict_Nums_human1, dict_Nums_human2)

if __name__ == '__main__':
    main()