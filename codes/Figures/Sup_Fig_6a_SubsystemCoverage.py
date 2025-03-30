import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial import distance
from matplotlib.colors import LinearSegmentedColormap

from proplot import rc
from pylab import *


def sub_plot(df):
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 7 

    # clustering of rows
    row_linkage = hierarchy.linkage(distance.pdist(df), method='average')

    # set color
    colors = ["#023858", "#3690c0", "#ffffff", "#fc9272", "#800026"]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    plt.figure(figsize=(8.75, 5))

    # plot
    g = sns.clustermap(
        df,
        figsize=(8.75, 5),
        row_linkage=row_linkage,
        col_cluster=False,
        cmap=cmap,
        center=df.values.mean(),
        linewidths=0,
        yticklabels=True,
        xticklabels=True,
        dendrogram_ratio=(0.05, 0),
        # cbar_pos=(1.02, 0.2, 0.03, 0.6), 
        cbar_pos=(0.2, 0.06, 0.3, 0.02),
        cbar_kws={"orientation": "horizontal"},
        tree_kws={'linewidths': 0.5}
    )

    g.ax_heatmap.set_xticklabels([])  
    g.ax_heatmap.set_yticklabels([])

    # set edgecolor
    g.ax_heatmap.add_patch(plt.Rectangle(
        (0, 0),
        len(df.columns),
        len(df.index),
        fill=False,
        edgecolor='black',
        lw=1
    ))

    # add ticks
    g.ax_heatmap.tick_params(axis='both', which='major', direction='in', length=1.5, width=0.4, color='black', labelcolor='black')
    g.ax_heatmap.tick_params(axis='both', which='minor', length=0)
    g.ax_heatmap.set_ylabel("")


    g.cax.tick_params(axis='x', which='major', direction='in', length=1.5, width=0.4, color='black', labelcolor='black')
    g.cax.tick_params(axis='x', which='minor', direction='in', length=0)



    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    minorticks_off()

    plt.tight_layout()
    
    # get cluster result
    row_idx = g.dendrogram_row.reordered_ind
    clustered_row_order = df.index[row_idx]
    print(clustered_row_order)

    # plt.savefig("../../results/Figures/Supple_Fig6_temp_1_2.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

    # save file
    # clustered_row_order.to_csv('clustered_row_order.csv', index=False)

def main():
    filepath = '../../results/'
    df = pd.read_csv(filepath + 'SubCoverage.tsv', sep='\t', index_col=0) # load data
    sub_plot(df)

if __name__ == "__main__":
    main()