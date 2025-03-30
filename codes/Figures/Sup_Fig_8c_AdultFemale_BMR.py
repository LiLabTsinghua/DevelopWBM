from scipy import stats
from scipy.stats import ranksums
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error,r2_score

from proplot import rc
from pylab import *


def Data_plot(x,y):
    True_energy = x
    predicted_energy = y

    correlation, p_value1 = stats.pearsonr(True_energy, predicted_energy)
    # p_value1 = stats.ttest_ind(True_energy, predicted_energy)
    r2 = r2_score(True_energy,predicted_energy)
    rmse = np.sqrt(mean_squared_error(True_energy,predicted_energy))

    print('The overall r is %.4f' % correlation)
    print('The overall P value is', p_value1)
    print('The overall R2 is %.4f' % r2)
    print('The overall RMSE is %.4f' % rmse)

    True_energy1 = True_energy
    predicted_energy1 = predicted_energy

    True_energy_ = np.array(True_energy1)
    parameter = np.polyfit(True_energy, predicted_energy, 1)
    r_line = parameter[0] * True_energy_  +  parameter[1]

    predicted_energy_ = np.array(predicted_energy1)

    correlation = np.corrcoef(predicted_energy_, r_line)[0,1] 
    R2 = correlation**2

    print('The line R2 is %.4f' % R2)
    
    fig, ax = plt.subplots(figsize=(3, 3))

    plt.rcParams.update({'font.size': 8})
    plt.rcParams['font.family'] = 'Arial'

    # Plot
    # ax.plot(True_energy, predicted_energy, 'o', markersize = 5, alpha=0.7)
    ax.plot(True_energy, predicted_energy, 'o', color = '#2171b5', markersize = 5, alpha=1)
    plt.plot(True_energy_, r_line, color='#000000',linewidth=1, linestyle = 'dotted')

    ax.grid(False)

    for i in ['bottom', 'top', 'right', 'left']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_visible('black')
        ax.spines[i].set_linewidth(0.5) 

    ax.set_xlim(1000, 1700)
    ax.set_ylim(1400, 1600)

    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)

    plt.tick_params(direction='in') 
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    plt.xlabel('Measured BMR (kcal)', fontsize=8)
    plt.ylabel('Predicted BMR by WBM (kcal)', fontsize=8)

    minorticks_off()

    # plt.savefig("../../results/Figures/Support_Fig_7c_female_BMR.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    df_corr = pd.read_csv('../../results/BMR/Result_AdultFemale_BMR.tsv', sep = '\t')
    x = df_corr['BMR_kcal_'].values.tolist()
    y = df_corr['PredictedBMR'].values.tolist()

    Data_plot(x,y)

if __name__ == '__main__':
    main()