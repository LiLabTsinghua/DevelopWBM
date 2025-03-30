import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import stats
import seaborn as sns
import pandas as pd
from scipy.stats import ranksums

from proplot import rc
from pylab import *


def main():
    mets = [0.7788, 0.885]
    models = ['Human1', 'Human2']
    colors = ['#d7191c', '#2c7bb6']

    plt.rcParams.update({'font.size': 8})
    plt.rcParams['font.family'] = 'Arial'

    fig, ax = plt.subplots(figsize=(2.1, 3.5))

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
        ax.bar(models[g], mets[g], color=colors[g], edgecolor = 'black', width=0.5)

    ax.set_ylim(0.6, 0.9)
    ax.yaxis.grid(True)
    ax.set_xticks([0, 1])
    plt.xticks(fontsize=8)

    ax.grid(False)
    minorticks_off()

    ax.set_ylabel('GPT check accuracy', fontdict={'size': 8})

    # plt.savefig("../../results/Figures/Fig_3b_GPT_Acc_0309.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

if __name__ == '__main__' :
    main()