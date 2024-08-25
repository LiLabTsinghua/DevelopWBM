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

def figure_dot_plot(df, Labels, Colors):
    # 设置图形大小
    fig, ax = plt.subplots(figsize=(6, 3))


    # 绘制每个系列的数据
    for i in range(len(Labels)):
        label_name = Labels[i]
        ax.plot(df['version'], df[label_name], marker='o', color=Colors[i], label=label_name)

    for i in ['bottom', 'top', 'right', 'left']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_visible('black')
        ax.spines[i].set_linewidth(0.5) 

    # 添加图例
    # ax.legend(fontsize=7, frameon = False, loc='upper left', bbox_to_anchor=(1.005, 1))
    ax.legend(fontsize=6, frameon = False)

    # 添加标题和标签
    # plt.title('Memote Score Over Versions')
    plt.xlabel('Version', fontsize=6)
    plt.ylabel('Number', fontsize=6)

    # 调整x轴标签角度
    plt.xticks(rotation=30, fontsize=6)
    plt.xticks(fontsize=6)

    plt.tick_params(direction='in') 
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    # 显示网格
    ax.set_axisbelow(True) 
    ax.grid(axis='y')

    minorticks_off()

    # 显示图形
    plt.savefig("D:/All_Human_GTEx/Papers_code_and_data/Figures/Fig_2d_temp.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    # load data
    filepath = 'C:\\Users\\DELL\\Pictures\\Fig\\Figure3\\'
    df_origin_data = pd.read_excel(filepath+'All_change_class.xlsx', sheet_name = 'Sheet3')
    data = {
        'version': df_origin_data['version'].values.tolist(),
        'Modify rxns': df_origin_data['fixed rxns'].values.tolist(),
        'Revise gpr': df_origin_data['revised gpr'].values.tolist(),
        'Add rxns': df_origin_data['rxn add'].values.tolist(),
        'Add mets': df_origin_data['met add'].values.tolist(),
        'Add gene': df_origin_data['gene add'].values.tolist(),
        'Remove rxns': df_origin_data['remove rxns(deplicate + deadend)'].values.tolist(),
        'Remove mets': df_origin_data['remove mets(duplicated + deadend)'].values.tolist(),
        'Remove genes': df_origin_data['gene remove'].values.tolist(),
    }

    Labels = ['Modify rxns', 'Revise gpr', 'Add rxns', 'Add mets', 'Add gene', 'Remove rxns', 'Remove mets', 'Remove genes']
    Colors = ['#d7191c', '#fdae61', '#66bd63', '#abd9e9', '#2c7bb6', '#80cdc1', '#de77ae', '#8073ac']
    # 创建数据框
    df = pd.DataFrame(data)

    figure_dot_plot(df, Labels, Colors)

if __name__ == '__main__' :
    main()