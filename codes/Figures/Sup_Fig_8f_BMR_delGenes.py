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


def data_categories(df_1):
    dict_1 = {'<95%':0,
            '95-97%':0,
            '97-98%':0,
            '98-99%':0,
            '99-100%':0,
            '100-101%':0,
            '101-102%':0,
            '102-103%':0,
            '103-105%':0,
            '>105%':0}

    refer_value = 1634.5


    decrease_genes = []
    increse_genes = []
    for i in df_1.index:
        if df_1.loc[i, 'BMR'] != 'NA':
            if df_1.loc[i, 'BMR'] < 0.95 * refer_value:
                dict_1['<95%'] += 1
                decrease_genes.append(i)
            elif 0.95 * refer_value <= df_1.loc[i, 'BMR'] < 0.97 * refer_value:
                dict_1['95-97%'] += 1
            elif 0.97 * refer_value <= df_1.loc[i, 'BMR'] < 0.98 * refer_value:
                dict_1['97-98%'] += 1
            elif 0.98 * refer_value <= df_1.loc[i, 'BMR'] < 0.99 * refer_value:
                dict_1['98-99%'] += 1
            elif 0.99 * refer_value <= df_1.loc[i, 'BMR'] < 1.0 * refer_value:
                dict_1['99-100%'] += 1
            elif 1.0 * refer_value <= df_1.loc[i, 'BMR'] < 1.01 * refer_value:
                dict_1['100-101%'] += 1
            elif 1.01 * refer_value <= df_1.loc[i, 'BMR'] < 1.02 * refer_value:
                dict_1['101-102%'] += 1
            elif 1.02 * refer_value <= df_1.loc[i, 'BMR'] < 1.03 * refer_value:
                dict_1['102-103%'] += 1
            elif 1.03 * refer_value <= df_1.loc[i, 'BMR'] < 1.05 * refer_value:
                dict_1['103-105%'] += 1
            elif 1.05 * refer_value <= df_1.loc[i, 'BMR']:
                dict_1['>105%'] += 1
                increse_genes.append(i)

    return dict_1

def BMR_plot(dict_1):
    categories = list(dict_1.keys())
    values = list(dict_1.values())

    fig, ax = plt.subplots(figsize=(3, 2))
    plt.rcParams.update({'font.size': 8})
    plt.rcParams['font.family'] = 'Arial'

    plt.bar(categories, values, color='#2171b5', edgecolor='black', linewidth=1)

    plt.xlabel('BMR Categories',fontsize=8)
    plt.ylabel('Gene counts',fontsize=8)

    ax.grid(False)
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    plt.tick_params(direction='in') 
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    plt.xticks(rotation=30, ha='right')

    ax.grid(False)
    minorticks_off()

    # plt.tight_layout()
    # plt.savefig(filepath+"../../results/Figures/Support_Fig_Del_BMR.pdf", transparent = True, dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    df_1 = pd.read_csv('../../results/AdultMale_BMR_del_genes.tsv', sep = '\t', index_col = 0)
    dict_1 = data_categories(df_1)
    BMR_plot(dict_1)

if __name__ == '__main__' :
    main()