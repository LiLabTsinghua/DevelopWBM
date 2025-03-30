import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import stats
import seaborn as sns
import pandas as pd
from scipy.stats import ranksums
import statistics

from proplot import rc
from pylab import *


def Depmap_vailin(data):

    plt.rcParams.update({'font.size': 8})
    plt.rcParams['font.family'] = 'Arial'
    # set colors
    fig, ax = plt.subplots(figsize=(3, 3.5))
    violin_color = ['#d7191c', '#2c7bb6']
    outer_line_color = '#000000'
    scatter_color = '#000000'
    scatter_size = 2

    # plot violin
    vplot = sns.violinplot(data=data, ax=ax, inner="box", cut=0, 
                        linewidth=1, edgecolor="black", palette=violin_color)


    for violin, color in zip(vplot.collections, violin_color):
        # violin.set_facecolor(color)
        violin.set_edgecolor('black')
        violin.set_linewidth(1)

    for i, d in enumerate(data):
        x = np.random.normal(i + 0, 0.1, len(d)) 
        ax.scatter(x, d, color='black', s=scatter_size)

    for i in ['bottom', 'top', 'right', 'left']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_visible('black')
        ax.spines[i].set_linewidth(0.5) 

    for i, d in enumerate(data):
        median = np.median(d)
        ax.hlines(median, i-0.2, i+0.2, color='#bababa', linewidth=2)

    ax.set_xticks([0,1])
    ax.set_xticklabels(['Human1', 'Human2'], fontsize=8)
    # ax.set_yticks([1])
    ax.set_ylabel('Metthews correlation coefficient', fontsize=8)
    ax.grid(False)

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    minorticks_off()

    # plt.savefig("../../results/Figures/Fig_3a_Depmap_Human1_Human2.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():

    filepath = '../../results/eGenes/'
    # filepath2 = 'D:\\Human_issue_check\\Human1.3_ftINIT\\'
    df_human1 = pd.read_csv(filepath+'results_Human1_DepMap.txt', sep = '\t')

    Human1_MCC = df_human1['MCC'].values.tolist()
    filtered_list2 = [x for x in Human1_MCC if not np.isnan(x)]
    median = statistics.median(filtered_list2)
    print(median)
    ndarray2 = np.array(filtered_list2)

    df_human2 = pd.read_csv(filepath+'results_Human2_Depmap.txt', sep = '\t')

    Human2_MCC = df_human2['MCC'].values.tolist()
    filtered_list3 = [x for x in Human2_MCC if not np.isnan(x)]
    median = statistics.median(filtered_list3)
    print(median)
    ndarray3 = np.array(filtered_list3)

    data = [ndarray2, ndarray3]

    Depmap_vailin(data)

if __name__== "__main__": 
    main()