import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from proplot import rc
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt

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
    mets = [0.2834, 0.4142, 0.4469, 0.605, 0.654]
    models = ['Recon3D', 'Human1', 'Human2', 'ecHuman1', 'ecHuman2']
    colors = ['#ffffbf', '#2c7bb6', '#d7191c', '#2c7bb6', '#d7191c']

    plt.rcParams.update({'font.size': 8})
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42

    fig, ax = plt.subplots(figsize=(2.8, 3.5))

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
        ax.bar(models[g], mets[g], color=colors[g], edgecolor = 'black')

    ax.set_ylim(0, 0.7)

    ax.yaxis.grid(True)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(models, fontsize=8, rotation=30, ha='right')

    ax.grid(False)
    minorticks_off()

    ax.set_ylabel('IEM accuracy', fontdict={'size': 8})

    # plt.savefig("../../results/Figures/Fig_3c_IEM_acc_0309.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

if __name__ == '__main__' :
    main()