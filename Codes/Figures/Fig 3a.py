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



def figure_violin(data):
    fig, ax = plt.subplots(figsize=(3, 3.5))

    plt.rcParams['font.family'] = 'Helvetica'

    # set colors
    violin_color = ['#1a9641', '#d7191c', '#2c7bb6']
    outer_line_color = '#000000'
    scatter_color = '#000000'
    scatter_size = 2.5

    # plot violin
    vplot = ax.violinplot(data, showmeans=False, showmedians=True,
                        showextrema=False, vert=True,
                        widths=0.7, bw_method=0.5)

    for patch, color in zip(vplot['bodies'], violin_color):
        patch.set_facecolor(color)


    for body in vplot['bodies']:
        body.set_edgecolor(outer_line_color)

    for i, d in enumerate(data):
        x = np.random.normal(i + 1, 0.1, len(d)) 
        ax.scatter(x, d, color=scatter_color, s=scatter_size, alpha=0.5)

    for i in ['bottom', 'top', 'right', 'left']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_visible('black')
        ax.spines[i].set_linewidth(0.5) 


    # ax.set_title('Violin Plot with Scattered Points')
    ax.yaxis.grid(True)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['Recon3D', 'Human 1', 'Human 2'], fontsize=6)
    # ax.set_yticks([1])
    ax.set_ylabel('Metthews correlation coefficient', fontsize=6)
    ax.set_yticklabels(np.array([0, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]), fontsize=6)

    ax.grid(False)

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    minorticks_off()

    plt.savefig("D:/All_Human_GTEx/Papers_code_and_data/Figures/Fig_3a_temp.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    # 生成数据
    # filepath = 'D:\\浏览器下载\\Human1_Publication_Data_Scripts\\Human1_Publication_Data_Scripts\\tINIT_GEMs\\gene_essentiality_outputs\\comparison_results\\'
    df_Recon3D = pd.read_csv(filepath+'results_Recon3D_DepMap_06.txt', sep = '\t')
    Recon3D_MCC = df_Recon3D['MCC'].values.tolist()
    filtered_list = [x for x in Recon3D_MCC if not np.isnan(x)]
    ndarray1 = np.array(filtered_list)

    # filepath = 'D:\\All_Human_GTEx\\Cell_line\\Depmap_results\\'
    df_human1 = pd.read_csv(filepath+'results_humanGEM_DepMap_06.txt', sep = '\t')
    Human1_MCC = df_human1['MCC'].values.tolist()
    filtered_list2 = [x for x in Human1_MCC if not np.isnan(x)]
    ndarray2 = np.array(filtered_list2)

    # filepath1 = 'D:\\All_Human_GTEx\\Cell_line\\Depmap_results\\'
    df_human2 = pd.read_csv(filepath+'results_Depmap_humanGEM_1.18.txt', sep = '\t')
    Human2_MCC = df_human2['MCC'].values.tolist()
    filtered_list3 = [x for x in Human2_MCC if not np.isnan(x)]
    ndarray3 = np.array(filtered_list3)

    data = [ndarray1, ndarray2, ndarray3]

    figure_violin(data)


if __name__== "__main__": 
    main()
