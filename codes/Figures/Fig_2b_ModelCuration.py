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

def figure_dot_plot(df, Labels, Colors):
    # set fig size
    fig, ax = plt.subplots(figsize=(5, 3))


    for i in range(len(Labels)):
        label_name = Labels[i]
        ax.plot(df['version'], df[label_name], marker='o', color=Colors[i], label=label_name)

    for i in ['bottom', 'top', 'right', 'left']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_visible('black')
        ax.spines[i].set_linewidth(0.5) 

    # add legend
    # ax.legend(fontsize=7, frameon = False, loc='upper left', bbox_to_anchor=(1.005, 1))
    ax.legend(fontsize=6, frameon = False)

    # add label
    # plt.title('Memote Score Over Versions')
    plt.xlabel('Version', fontsize=6)
    plt.ylabel('Number', fontsize=6)

    plt.xticks(rotation=30, fontsize=6)
    plt.xticks(fontsize=6)

    plt.tick_params(direction='in') 
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    # displays the X-axis grid
    ax.set_axisbelow(True) 
    ax.grid(axis='y')

    minorticks_off()

    # show the fig
    # plt.savefig("../../results/Figures/Fig_2b.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    # load data
    filepath = '../../data/All_model_changes.tsv'
    df_origin_data = pd.read_csv(filepath, sep='\t')
    data = {
        'version': df_origin_data['version'].values.tolist(),
        'Modify rxns': df_origin_data['fixed rxns'].values.tolist(),
        'Revise gpr': df_origin_data['revised gpr'].values.tolist(),
        'Add rxns': df_origin_data['rxn add'].values.tolist(),
        'Add mets': df_origin_data['met add'].values.tolist(),
        'Add gene': df_origin_data['gene add'].values.tolist(),
        'Remove rxns': df_origin_data['remove rxns(deplicate + deadend)'].values.tolist(),
        'Remove mets': df_origin_data['remove mets(duplicated + deadend)'].values.tolist(),
        'Remove genes': df_origin_data['gene remove'].values.tolist(),
    }

    Labels = ['Modify rxns', 'Revise gpr', 'Add rxns', 'Add mets', 'Add gene', 'Remove rxns', 'Remove mets', 'Remove genes']
    Colors = ['#d7191c', '#fdae61', '#66bd63', '#abd9e9', '#2c7bb6', '#80cdc1', '#de77ae', '#8073ac']
    df = pd.DataFrame(data)

    figure_dot_plot(df, Labels, Colors)

if __name__ == '__main__' :
    main()