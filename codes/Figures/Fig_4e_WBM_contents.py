import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from proplot import rc
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import stats
import seaborn as sns
import pandas as pd
from scipy.stats import ranksums

from proplot import rc
from pylab import *


def Plot_bar(data_list, content, Y_slim):
    plt.rcParams.update({'font.size': 6})
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42

    models = ['Adult male', 'Elderly male', 'Adult female','Elderly female', 'Fetus']

    colors = ['#2c7bb6', '#abd9e9', '#d7191c', '#fdae61', '#bababa']

    fig, ax = plt.subplots(figsize=(1.8, 1.8))

    fig.patch.set_alpha(0)

    ax.grid(False)
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    plt.tick_params(direction='in') 
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    for g in range(len(models)):
        ax.bar(models[g], data_list[g], color=colors[g], edgecolor = 'black')

    ax.set_ylim(Y_slim[0], Y_slim[1])
    # ax.set_xticklabels(fontsize=6)

    ax.yaxis.grid(True)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(models, fontsize=6, rotation=30, ha='right')
    # plt.xticks(fontsize=8, rotation=30)
    # plt.yticks(fontsize=8)

    ax.grid(False)
    minorticks_off()


    plt.ticklabel_format(axis='y', style='sci', scilimits=(-1,2))
    # ax.set_yticklabels(np.array([0, 1, 2, 3]), fontsize=6)

    ax.set_ylabel(content, fontdict={'size': 6})

    # plt.savefig("../../results/Figures/Fig_WBM_genes.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()


def main():
    content = 'Genes' 

    # Because WBM is so large that it takes a long time to load, here is the WBM content sorted out in advance.
    dict_WBM = {'Genes':[2418, 2419, 2417, 2417, 2417],
                'Metabolites':[125960, 125823, 137419, 137339, 115679],
                'Reactions':[157479, 157273, 171813, 171616, 145005]}
    dict_Y_slim = {'Genes':[0,2600],
                   'Metabolites':[0,160000],
                   'Reactions':[0,190000]}
    
    Plot_bar(dict_WBM[content], content, dict_Y_slim[content])

if __name__ == "__main__":
    main()
    