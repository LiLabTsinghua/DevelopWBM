import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import stats
import seaborn as sns
import pandas as pd
from scipy.stats import ranksums

from proplot import rc
from pylab import *


def figure_barh2(gem_info, number1, number2):
        
    colors = ['#fb746c', '#629afa', '#629afa', '#629afa', '#04b031', '#04b031']

    # rc["font.family"] = "Helvetica"
    plt.rcParams.update({'font.size': 6})
    # plt.rcParams.update({'font.size': 7})
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype']= 42

    fig, ax = plt.subplots(figsize=(1.8, 3))

    fig.patch.set_alpha(0)

    ax.grid(False)
    ax = plt.gca()
    # plt.legend(loc='upper right', fontsize=6, frameon=True, fancybox=True, framealpha=0.2, borderpad=0.3,
    #            ncol=1, markerfirst=True, markerscale=1, numpoints=1, handlelength=3.5)

    ax.barh(gem_info[0], number1[0], 0.5, color=colors[0], edgecolor = 'black', label = 'genes')

    ax.barh(gem_info[1], number1[1], 0.5, color=colors[1], edgecolor = 'black', label = 'rxns')

    ax.barh(gem_info[2], number1[2], 0.5, color=colors[2], edgecolor = 'black')

    ax.barh(gem_info[3], number1[3], 0.5, color=colors[3], edgecolor = 'black')

    ax.barh(gem_info[4], number1[4], 0.5, color=colors[4], edgecolor = 'black', label = 'mets')

    ax.barh(gem_info[5], number1[5], 0.5, color=colors[5], edgecolor = 'black')

    plt.legend(loc='best', fontsize=6, bbox_to_anchor=(0.535,1), frameon = False)

    # autolabel(ax)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    my_x_ticks = np.arange(0, 25000, 5000)
    plt.xticks(my_x_ticks)
    # plt.yticks(my_y_ticks)

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)


    plt.xticks(fontsize=6, rotation=30, ha='right')
    plt.yticks(fontsize=6)

    minorticks_off()

    # plt.savefig("../../Figures/Fig_2c_temp.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    number1 = [2890, 12995, 8066, 1658, 8456, 4156]
    number2 = [3627, 13520, 8321, 0, 10103, 4166]
    gem_info = ["Total genes", "Total rxns", "Gene-associated rxns", "Exchange rxns", "Total mets", "Unique mets"]
    figure_barh2(gem_info, number1, number2)

if __name__ == '__main__' :
    main()