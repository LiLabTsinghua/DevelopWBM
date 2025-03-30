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


def figure_bar(df_1, df_2):

    # check dataFrames
    assert df_1.index.equals(df_2.index), "The indices of the two DataFrames are not equal."
    assert df_1.columns.equals(df_2.columns), "The columns of the two DataFrames are not equal."

    metrics = ['accuracy', 'sensitivity', 'specificity', 'F1', 'MCC']
    cellLines = ['DLD1', 'GBM', 'HCT116', 'HELA', 'RPE1']

    y_min = [0.6, 0, 0.7, 0, 0]
    y_max = [0.9, 0.5, 1, 0.5, 0.5]

    # Create graphs and axes
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(7, 8.5))
    plt.rcParams.update({'font.size': 8})
    plt.rcParams['font.family'] = 'Arial'

    for idx, ax in enumerate(axs.flat):
        if idx < len(metrics):
            metric = metrics[idx]
            human1_data = df_1[metric].values.tolist()[:-1]
            human2_data = df_2[metric].values.tolist()[:-1]

            # set width of bar
            bar_width = 0.35
            ind = np.arange(len(cellLines))
            r1 = [x - (bar_width/2) for x in ind]
            r2 = [x + (bar_width/2) for x in ind]

            # Draw a bar chart
            ax.bar(r1, human1_data, bar_width, color='#d73027', edgecolor = 'black', label = 'Human1')
            ax.bar(r2, human2_data, bar_width, color='#4575b4', edgecolor = 'black', label = 'Human2')

            # Add some labels
            if metric != 'MCC':
                ax.set_ylabel(metric.capitalize())
            else:
                ax.set_ylabel(metric)

            ax.set_xticks(ind)
            ax.set_xticklabels(tuple(cellLines))

            ax.spines['bottom'].set_linewidth(0.5)
            ax.spines['left'].set_linewidth(0.5)
            ax.spines['top'].set_linewidth(0.5)
            ax.spines['right'].set_linewidth(0.5)

            my_y_ticks = np.arange(y_min[idx], y_max[idx], 0.1)
            ax.set_ylim(y_min[idx], y_max[idx])
            ax.set_yticks(my_y_ticks)

            ax.tick_params(direction='in')
            ax.tick_params(which='major',length=1.5)
            ax.tick_params(which='major',width=0.4)

            ax.legend(ncol=2, fontsize=6, frameon=False)

            ax.grid(False)
            ax.minorticks_off()

            filepath = 'D:\\All_Human_GTEx\\Human2_Paper_and_Suplementary\\New_Human_GEM\\results\\Figures\\'
            # plt.savefig(filepath+"Supple_Fig3_f_temp_2.pdf", dpi=400, bbox_inches = 'tight')
        else:
            # fig.delaxes(axs[2, 1])
            ax.axis('off')

    # show the plot
    plt.show()


def main():
    datasetName = 'Hart2015'
    df_1 = pd.read_csv('../../results/eGenes/Human1_Hart2015_essentialGenes_ftINIT.tsv', index_col = 0, sep = '\t')
    df_2 = pd.read_csv('../../results/eGenes/Human2_Hart2015_essentialGenes_ftINIT.tsv', index_col = 0, sep = '\t')

    figure_bar(df_1, df_2)

if __name__ == '__main__' :
    main()