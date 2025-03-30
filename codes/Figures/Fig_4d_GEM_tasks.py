import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from proplot import rc
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt

from proplot import rc
from pylab import *

def Task_plot(harvest, list1, list2):
    model1 = list1
    model2 = list2

    plt.rcParams.update({'font.size': 6})
    plt.rcParams['font.family'] = 'Arial'

    # plt.figure(figsize=(0.5, 0.4))
    fig, ax = plt.subplots(figsize=(1.7, 1.6))
    im = ax.imshow(harvest, cmap='Blues')
    cbar = fig.colorbar(im, fraction=0.044, pad=0.05)

    plt.xticks(np.arange(len(model1)), labels=model1, fontsize=6, 
                        rotation=30, rotation_mode="anchor", ha="right")
    plt.yticks(np.arange(len(model2)), labels=model2, fontsize=6)

    for i in range(len(model1)):
        for j in range(len(model1)):
            if harvest[i, j] > 0.6:
                text = plt.text(j, i, "%.2f" % harvest[i, j], ha="center", va="center", color="w", fontsize=6, fontweight = 'bold')
            else:
                text = plt.text(j, i, "%.2f" % harvest[i, j], ha="center", va="center", color="black", fontsize=6, fontweight = 'bold')

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

    # plt.savefig("D:/All_Human_GTEx/Human2_Paper_and_Suplementary/codes/Figures/Fig4e_hotmap_2.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()


def main():
    filepath = '../../results/metabolicTasks_comparation/'
    Sex = ['Male', 'Female']
    Tissue = ['Liver', 'Heart', 'Kidney', 'Adipocytes']

    Models = ['Harvey', 'Harvetta', 'Adult Male', 'Adult Female']
    df = pd.DataFrame(index=Tissue, columns=Models)
    for i in Sex:
        df_temp = pd.read_csv(filepath+'Adult'+i+'_MetabolicTasks.tsv',sep = '\t', index_col = 0)
        for j in Tissue:
            if i == 'Male':
                df.loc[j, 'Adult '+i] = float(df_temp.loc['Human2', j])
                df.loc[j, 'Harvey'] = float(df_temp.loc['Recon3D', j])
            elif i == 'Female':
                df.loc[j, 'Adult '+i] = float(df_temp.loc['Human2', j])
                df.loc[j, 'Harvetta'] = float(df_temp.loc['Recon3D', j])

    harvest = np.array([df.iloc[0].tolist(), 
                        df.iloc[1].tolist(), 
                        df.iloc[2].tolist(),
                        df.iloc[3].tolist()])
    Task_plot(harvest, Models, Tissue)

if __name__ == '__main__' :
    main()
