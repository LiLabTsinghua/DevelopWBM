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


def Fig_flux(df_1):
    mets = df_1['Metabolite'].values.tolist()
    Measured = df_1['Measured'].values.tolist()
    Predicted  = df_1['Predicted'].values.tolist()
    Cell_lines = df_1['Cell.Line'].values.tolist()
    colors  = df_1['colors'].values.tolist()


    Measured_pre = []
    for i in Measured:
        if i != 0:
            Measured_pre.append(math.log10(abs(i)))
        else:
            Measured_pre.append(-8)

    Predicted_pre = []
    for j in Predicted:
        if j != 0:
            Predicted_pre.append(math.log10(abs(j)))
        else:
            Predicted_pre.append(-8)

    x = np.arange(-8, 0, 0.1)
    y = x

    fig, ax = plt.subplots(figsize=(3, 3))

    plt.rcParams['font.family'] = 'Helvetica'

    fig.patch.set_alpha(0)

    ax.grid(False)
    ax = plt.gca()
    for i in range(len(Measured)):
        if i == 0 or i == 27 or i == 54 or i == 81 or i == 108 or i == 135 or i == 162 or i == 189 or i == 216 or i == 243 or i == 270:
            plt.scatter(Measured_pre[i], Predicted_pre[i], s = 8, alpha = 0.8, color = colors[i], label = Cell_lines[i])
        else:
            plt.scatter(Measured_pre[i], Predicted_pre[i], s = 8, alpha = 0.8, color = colors[i])

        ax.legend(bbox_to_anchor=(-0.05, -0.3), 
                loc='lower left',
                ncol=3,
                labelspacing=0,
                columnspacing = 0.4,
                fontsize=6,
                frameon=False)

    plt.plot(x, y, color='#000000',linewidth=1)



    for i in ['bottom', 'top', 'right', 'left']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_visible('black')
        ax.spines[i].set_linewidth(0.5) 

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    plt.xlabel('Measured flux (log10 abs. val.)', fontsize=6)
    plt.ylabel('Predicted flux (log10 abs. val.)', fontsize=6)

    my_x_ticks = np.arange(-8, 2, 2)
    my_y_ticks = np.arange(-8, 2, 2)
    plt.xticks(my_x_ticks, fontsize=6)
    plt.yticks(my_y_ticks, fontsize=6)

    minorticks_off()

    plt.savefig("D:/All_Human_GTEx/Papers_code_and_data/Figures/Fig_3c_temp.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    filepath = 'D:\\All_Human_GTEx\\EC_models\\cell_line\\simulation_dir\\Results\\'
    df_1 = pd.read_excel(filepath+'human2_ecModel_full_changedMass.xlsx')

    Fig_flux(df_1)

if __name__ == '__main__' :
    main()