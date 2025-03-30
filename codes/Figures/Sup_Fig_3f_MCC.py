import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from proplot import rc
from pylab import *

def plot_MCC(df):
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["#d7191c", "#2c7bb6"])
    norm = plt.Normalize(df["MCC"].min(), df["MCC"].max())

    fig, ax = plt.subplots(figsize=(3.2, 2.5))
    plt.rcParams.update({'font.size': 8})
    plt.rcParams['font.family'] = 'Arial'

    bars = ax.bar(df["Version"], df["MCC"], color=cmap(norm(df["MCC"])))

    plt.xticks(rotation=45, ha="right", fontsize=6)
    plt.yticks(fontsize=8)
    plt.xlabel("Version", fontsize=8)
    plt.ylabel("MCC", fontsize=8)

    ax.tick_params(direction='in')
    ax.tick_params(which='major',length=1.5)
    ax.tick_params(which='major',width=0.4)

    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    ax.set_ylim(0.1, 0.325)

    ax.grid(False)

    ax.minorticks_off()

    # plt.savefig("../../results/Figures/Supple_Fig_Allversions_MCC.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    filepath = '../../results/eGenes/'
    df_1 = pd.read_csv(filepath+'HumanGEM_allVersions_MCC.tsv', index_col = 0, sep = '\t')
    data = {}
    data['Version'] = [str(i) for i in df_1.index.tolist()]
    data['MCC'] = df_1['MCC'].values.tolist()

    df = pd.DataFrame(data)

    plot_MCC(df)

if __name__ == '__main__' :
    main()