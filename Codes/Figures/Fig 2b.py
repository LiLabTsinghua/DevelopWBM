
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

from collections import Counter


def Compartments_Statistics(df_location):
    Genes = df_location['genes'].values.tolist()
    Compartments = df_location['compartments'].values.tolist()

    All_Genes = []
    All_Comps = []
    for i in range(len(Compartments)):
        gene = Genes[i]
        if str(Compartments[i]) != 'nan':
            comps = Compartments[i].split(';')
        
            for j in comps:
                All_Genes.append(gene)
                All_Comps.append(j)
        else:
            All_Genes.append(gene)
            All_Comps.append('Others')

    count = Counter(All_Comps)
    # sorted_dict = sorted(count.items(), key=lambda x: x[1])
    All_comps = sorted(count, key=count.get)
    All_nums = sorted(count.values())

    return All_comps, All_nums

def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(pct*total/100.0)
        # pct_distance = 0.6 + pct * 0.2
        return '{v:d}'.format(v=val)
    return my_autopct


def figure_pie(All_comps, All_nums):
    # Create a pie chart
    plt.rcParams.update({'font.size': 6})
    # plt.rcParams.update({'font.size': 7})
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype']= 42
    fig, ax = plt.subplots(figsize=(2, 2))

    colors = ['#f8746b', '#d19100', '#91a700', '#00b834', '#00bf9e', '#00b8e3', '#609cff', '#d870f9', '#ff60c1']
    wedgeprops = {'edgecolor': 'black', 'linewidth': 1}
    label_distances = [1.5 if num < 166 else 1.1 for num in All_nums]
    # ax.pie(All_nums, labels=None, autopct=make_autopct(All_nums), colors = colors, labeldistance=label_distances, startangle=90, wedgeprops=wedgeprops)
    ax.pie(All_nums, labels=None, autopct=make_autopct(All_nums), colors = colors, pctdistance = 0.6, labeldistance=10, startangle=90, wedgeprops=wedgeprops)
    ax.axis('equal')  # Equal aspect ratio ensures the pie chart is circular
    # ax.set_title('Protein location')

    # Add a legend
    plt.legend(All_comps, loc='center left',
            labelspacing=0,
            columnspacing = 0.4,
            fontsize=6,
            prop={'size': 6, 'family': 'Arial'},
            bbox_to_anchor=(0, -0.15),
            ncol=2,
            frameon=False)

    # plt.savefig("../../Figures/Fig_2b.pdf", transparent = True, dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    # data_preparation
    df_location = pd.read_csv('../../Data/'+'genes.tsv', sep = '\t')
    All_comps, All_nums = Compartments_Statistics(df_location)

    # plot
    figure_pie(All_comps, All_nums)

if __name__ == '__main__' :
    main()