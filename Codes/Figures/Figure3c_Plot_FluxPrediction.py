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


def process_flux_data(file_path):
    # Load experimental and predicted exchange fluxes
    fluxes = pd.read_csv(file_path, sep='\t')
    fluxes = fluxes.drop('exchangeIDs', axis=1)  # Remove unneeded column
    
    samples = fluxes.columns
    sample_types = pd.Categorical(['Predicted' if s.startswith('exp_') else 'Measured' for s in samples])
    
    # Rearrange data
    fluxdata = []
    for i in range(1, len(fluxes.columns), 2):
        add_data = fluxes.iloc[:, [0, i, i+1]]
        add_data['Cell.Line'] = fluxes.columns[i+1]
        add_data.columns = ['Metabolite', 'Measured', 'Predicted', 'Cell.Line']
        fluxdata.append(add_data)
    
    fluxdata = pd.concat(fluxdata, ignore_index=True)
    fluxdata['Cell.Line'] = fluxdata['Cell.Line'].astype('category')
    fluxdata['Metabolite'] = fluxdata['Metabolite'].astype('category')
    
    # Log-transform fluxes
    log_fluxdata = fluxdata.copy()
    log_fluxdata.iloc[:, 1:3] = np.log10(np.abs(fluxdata.iloc[:, 1:3]))
    log_fluxdata.iloc[:, 1:3] = log_fluxdata.iloc[:, 1:3].replace([np.inf, -np.inf], -8)
    
    return fluxdata, log_fluxdata

def Fig_flux(df_1):
    Measured = df_1['Measured'].values.tolist()
    Predicted  = df_1['Predicted'].values.tolist()
    Cell_lines = df_1['Cell.Line'].values.tolist()
    # define cell colors
    dict_colors = {'HS_578T':'#f44336',
                   'RPMI_8226':'#c2185b',
                   'HT29':'#ab47bc',
                   'MALME_3M':'#4527a0',
                   'SR':'#3f51b5',
                   'UO_31':'#1e88e5',
                   'MDMAMB_231':'#4db6ac',
                   'HOP62':'#2e7d32',
                   'NCI_H226':'#c0ca33',
                   'HOP92':'#f9a825',
                   'O_786':'#795548'}

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

    plt.rcParams.update({'font.size': 6})
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype']= 42

    fig.patch.set_alpha(0)

    ax.grid(False)
    ax = plt.gca()
    for i in range(len(Measured)):
        if i == 0 or i == 27 or i == 54 or i == 81 or i == 108 or i == 135 or i == 162 or i == 189 or i == 216 or i == 243 or i == 270:
            plt.scatter(Measured_pre[i], Predicted_pre[i], s = 8, alpha = 1, color = dict_colors[Cell_lines[i]], label = Cell_lines[i])
        else:
            plt.scatter(Measured_pre[i], Predicted_pre[i], s = 8, alpha = 1, color = dict_colors[Cell_lines[i]])

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

    plt.savefig("../../Figures/Fig_3c.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():

    # load the predicted results from "predict_cellLines_gRates.m"
    fluxdata, log_fluxdata = process_flux_data('../../Results/11_cellLines_NCI60/'+'ecModels_const_0_exchangeFluxesComp.txt')
    Fig_flux(fluxdata)

if __name__ == '__main__' :
    main()

