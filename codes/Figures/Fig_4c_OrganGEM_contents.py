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


def plot_figure(alldata, Ticks, content):

    colors = ['#2c7bb6', '#abd9e9', '#d7191c', '#fdae61', '#bababa',  '#f4a582', '#8dd3c7']

    # Create canvas
    fig, ax = plt.subplots(figsize=(1.9, 1.9))
    plt.grid(False) # No grid

    plt.rcParams.update({'font.size': 6})
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype']= 42

    # Draw a box diagram
    sns.boxplot(data=alldata, showmeans=False, width=0.6, linewidth=0.8, boxprops=dict(facecolor="None"))

    # Mark point
    for column in alldata:
        column_r = alldata.columns.get_loc(column)
        ax.plot(np.random.normal(column_r, 0.05, size=len(alldata[column])), alldata[column], 'o', color = colors[column_r], alpha=0.8, markersize='5')
        # ax.plot(np.random.normal(column_r, 0.02, size=10), alldata[column], 'o', alpha=0.3)

    # for i, column in enumerate(alldata.columns):
    #     plt.setp(ax.artists[i], color=colors[i], alpha=0.6)

    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5) 
    # ax.spines['bottom'].set_color('black') 
    ax.spines['left'].set_linewidth(0.5)
    # ax.spines['left'].set_color('black') 
    ax.spines['top'].set_linewidth(0.5)
    # ax.spines['top'].set_color('black') 
    ax.spines['right'].set_linewidth(0.5)
    # ax.spines['right'].set_color('black') 


    # set ticks 
    ax.set_yticks(Ticks)

    plt.xticks(fontsize=6, rotation=30, ha='right')
    plt.yticks(fontsize=6)

    ax.tick_params(axis='x', direction='in', length=3)
    ax.tick_params(axis='y', direction='in', length=3)

    # set labels
    ax.set_ylabel(content+' counts', fontdict={'size': 6})


    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    minorticks_off()

    # plt.savefig("../../results/Figures/Fig_4_Organ_"+content+".pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    content = 'Reactions'

    filepath = '../../data/'
    df = pd.read_excel(filepath + 'OrganModel_contents.xlsx', sheet_name=content, index_col=0)

    # Considering that the fetus has only 10 representative organ metabolite models, only the contents of these 10 organ-specific models are compared here.
    Common_organs = ['Brain', 'Heart', 'Kidney', 'Liver', 'Lung', 'Muscle', 'Pancreas', 'sIEC', 'Spleen', 'Stomach']

    dict_WBM_rxns = {}

    for i in df.columns.to_list():
        organ_rxns = []
        for j in Common_organs:
            organ_rxns.append(df.loc[j, i])

        dict_WBM_rxns[i] = organ_rxns

    alldata = pd.DataFrame(dict_WBM_rxns)

    if content == 'Genes':
        Ticks = [800,1200,1600,2000,2400]
    elif content == 'Metabolites':
        Ticks = [2000,3000,4000,5000,6000]
    elif content == 'Reactions':
        Ticks = [2000,4000,6000,8000,10000]

    plot_figure(alldata, Ticks, content)


if __name__== "__main__": 
    main()