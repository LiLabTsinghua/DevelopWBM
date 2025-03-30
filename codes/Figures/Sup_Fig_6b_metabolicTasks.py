import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial import distance
from matplotlib.colors import LinearSegmentedColormap

from proplot import rc
from pylab import *


def Tasks_plot(df):
    rows = df.index.values.tolist()
    columns = df.columns.values.tolist()

    data = df.values

    plt.rcParams.update({'font.size': 8})
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype']= 42

    fig, ax = plt.subplots(figsize=(6, 5))


    # plot dots
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i, j]:
                ax.scatter(j, data.shape[0]-1-i, s=6, color='#2c7bb6')

    ax.set_xlim(-2, len(columns)+1)
    ax.set_ylim(-1, len(rows)+1)

    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([''] * len(columns))
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([''] * len(rows))

    ax.tick_params(axis='both',
                   direction='in',
                   length=1.5,
                   width=0.4)

    ax = plt.gca()
    for side in ['left', 'right', 'top', 'bottom']:
        ax.spines[side].set_color('black')
        ax.spines[side].set_linewidth(0.5)

    plt.grid(False)

    minorticks_off()

    plt.tight_layout()
    # plt.savefig('../../results/Figures/Support_Fig_Organ_models_MetabolicTasks.pdf', dpi=400, bbox_inches='tight')
    plt.show()

def main():
    filepath = '../../results/'
    df = pd.read_csv(filepath + 'FullTasks_Results.tsv', sep='\t', index_col=0) # load data
    Tasks_plot(df)

if __name__ == "__main__":
    main()