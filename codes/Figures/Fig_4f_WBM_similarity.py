import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt

from proplot import rc
from pylab import *


def Similarity_heatMap(harvest,modelNames_1,modelNames_2):
    fig, ax = plt.subplots(figsize=(1.7, 1.6))
    im = ax.imshow(harvest, cmap='OrRd')
    cbar = fig.colorbar(im, fraction=0.044, pad=0.05)

    model1 = modelNames_1
    model2 = modelNames_2
    plt.xticks(np.arange(len(model1)), labels=model1, fontsize=7, 
                         rotation=30, rotation_mode="anchor", ha="right")
    plt.yticks(np.arange(len(model2)), labels=model2, fontsize=7)     
    # plt.title("Harvest of local farmers (in tons/year)")

    for i in range(len(model1)):
        for j in range(len(model1)):
            if i == j:
                text = plt.text(j, i, "%.2f" % harvest[i, j], ha="center", va="center", color="w", fontsize=7, fontweight = 'bold')
            else:
                text = plt.text(j, i, "%.2f" % harvest[i, j], ha="center", va="center", color="black", fontsize=7, fontweight = 'bold')

    plt.grid(False)

    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(0.5)

    cbar.minorticks_off()

    ax = plt.gca()
    for side in ['left', 'right', 'top', 'bottom']:
        ax.spines[side].set_color('black')
        ax.spines[side].set_linewidth(0.5)

    plt.tick_params(direction='in') 
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    minorticks_off()

    # plt.savefig("D:/All_Human_GTEx/Papers_code_and_data/Figures/Fig4e_hotmap.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()
    
def main():
    df = pd.read_csv('../../results/'+'WBM_Similarity.tsv',sep = '\t', index_col = 0)
    modelNames_1 = df.index.tolist()
    modelNames_2 = df.columns.tolist()
    harvest = df.to_numpy()
    Similarity_heatMap(harvest,modelNames_1,modelNames_2)

if __name__ == '__main__' :
    main()
