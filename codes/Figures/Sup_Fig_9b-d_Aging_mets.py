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

def change_mets(df_adult, df_old, df_class):
    adult_mets = df_adult.index.tolist()
    old_mets = df_old.index.tolist()

    both_mets = [value for value in adult_mets if value in set(old_mets)]

    dict_met_class = dict(zip(df_class['WBM_met'].values.tolist(), df_class['super_class'].values.tolist()))

    uni_class = list(set(df_class['super_class'].values.tolist()))
    change_type = ['increase', 'stability', 'decrease']
    df_1 = pd.DataFrame(0, index=uni_class, columns=change_type)

    dict_change_mets = {}
    for i in both_mets:
        if i in dict_met_class:
            adult_value = df_adult.loc[i, 'DemandValue']
            old_value = df_old.loc[i, 'DemandValue']
            # print(adult_value)
            try:
                if not np.isnan(adult_value) and not np.isnan(old_value):
                    if abs(old_value) - abs(adult_value) > 1e-6:
                        dict_change_mets[i] = 'increase'
                    elif abs(old_value) - abs(adult_value) < -1e-6:
                        dict_change_mets[i] = 'decrease'
                    elif -1e-6 < abs(old_value) - abs(adult_value) < 1e-6:
                        dict_change_mets[i] = 'stability'
            except:
                print(adult_value)
                print(old_value)


    for j in uni_class:
        for k, v in dict_change_mets.items():
            if dict_met_class[k] == j:
                df_1.loc[j, v] = df_1.loc[j, v] + 1

    return df_1

def Aging_mets(df):

    colors = {"decrease": "#2c7bb6", "stability": "#bababa", "increase": "#d7191c"}
    labels = ["increase", "stability", "decrease"]

    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    left = [0] * len(df)

    for label in labels:
        ax.barh(
            df.index, 
            df[label], 
            left=left, 
            color=colors[label], 
            label=label,
            height=0.6 
        )
        left = [sum(x) for x in zip(left, df[label])]

    ax.set_xlim(0, 100)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    ax.set_xlabel("Metabolite Change Ratio", fontname='Arial', fontsize=6)

    ax.yaxis.set_tick_params(labelsize=6)

    ax.tick_params(
        axis='both',
        which='major',
        direction='in',
        length=1.5,
        width=0.4,
        bottom=True,
        left=False
    )

    legend = ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.25),
        ncol=3,
        frameon=False
    )

    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    ax.grid(False)

    minorticks_off()

    plt.tight_layout()
    # plt.savefig(filepath+"Support_Fig_8_Male_u_2.pdf", transparent = True, dpi=400, bbox_inches = 'tight')
    plt.show()


def main():
    df_adult = pd.read_csv('../../results/Aging_mets/AdultMale_bc_MetsDemand.tsv', sep = '\t', index_col=0)
    df_old = pd.read_csv('../../results/Aging_mets/ElderlyMale_bc_MetsDemand.tsv', sep = '\t', index_col=0)

    df_class = pd.read_csv('../../data/AgingMets_classify.tsv', sep='\t')

    df_1 = change_mets(df_adult, df_old, df_class)
    df_1 = df_1[df_1.sum(axis=1) >= 20] # Metabolite types with a number of less than 20 are not considered
    df_2 = df_1.apply(lambda row: (row / row.sum()) * 100, axis=1)
    df_3 = df_2.sort_values(by='increase', ascending=True)

    Aging_mets(df_3)

if __name__ == "__main__":
    main()
