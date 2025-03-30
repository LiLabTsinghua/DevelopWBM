# 散点图
from scipy import stats
from scipy.stats import ranksums
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error,r2_score

from proplot import rc
from pylab import *


def Food_energy(Food_id, True_energy, predicted_energy):
    
    correlation, p_value = stats.pearsonr(True_energy, predicted_energy)
    p_value1 = stats.ttest_ind(True_energy, predicted_energy)
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

    plt.figure(figsize=(2,2))

    plt.grid(False)

    rc('font',**{'family':'serif','serif':['Arial']})
    plt.rcParams['pdf.fonttype'] = 42

    # plt.axes([0.12,0.12,0.83,0.83])

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    True_energy1 = []
    predicted_energy1 = []
    
    # potato_id = '31036'
    # potato_energy_x = []
    # potato_predicted_energy = []

    Cornstarch_id = '20027'
    Cornstarch_energy_x = []
    Cornstarch_predicted_energy = []

    Fish = '15018'
    Fish_energy_x = []
    Fish_predicted_energy = []

    Butter_id = '1145'
    Butter_energy_x = []
    Butter_predicted_energy = []

    for f in range(len(True_energy)):
        id = Food_id[f]
        Te = True_energy[f]
        Pe = predicted_energy[f]
        if id == Cornstarch_id:
            Cornstarch_energy_x.append(Te)
            Cornstarch_predicted_energy.append(Pe)
        elif id == Fish:
            Fish_energy_x.append(Te)
            Fish_predicted_energy.append(Pe)
        elif id == Butter_id:
            Butter_energy_x.append(Te)
            Butter_predicted_energy.append(Pe)
        else:
            True_energy1.append(Te)
            predicted_energy1.append(Pe)


    # plt.scatter(True_energy1, predicted_energy1, s = 2, alpha = 0.4, color = '#2c7bb6')
    plt.scatter(True_energy1, predicted_energy1, s = 2, alpha = 1, color = '#8e9083')
    
    plt.scatter(Cornstarch_energy_x, Cornstarch_predicted_energy, s = 10, edgecolors='black', linewidths=0.5, alpha = 1, color = '#fdae61')
    plt.scatter(Fish_energy_x, Fish_predicted_energy, s = 10, edgecolors='black', linewidths=0.5, alpha = 1, color = '#b2182b')
    plt.scatter(Butter_energy_x, Butter_predicted_energy, s = 10, edgecolors='black', linewidths=0.5, alpha = 1, color = '#762a83')
    plt.plot(True_energy_, r_line, color='#000000',linewidth=1, linestyle = 'dotted')

    plt.xlabel('Energy(Kcal)', fontsize=6)
    plt.ylabel('Predict ATP(mmol/gDCW/h)', fontsize=6)

    my_x_ticks = np.arange(0, 1000, 200)
    my_y_ticks = np.arange(0, 500, 100)
    plt.xticks(my_x_ticks, fontsize=6)
    plt.yticks(my_y_ticks, fontsize=6)

    minorticks_off()

    plt.savefig("../../results/Figures/Fig_5_Food_corr_1.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def main():
    filepath = 'D:\\All_Human_GTEx\\Human2_Paper_and_Suplementary\\New_Human_GEM\\results\\'
    df_food = pd.read_csv(filepath+'Result_FoodEnergy_Prediction.tsv', sep = '\t')
    temp_Food_id = df_food['food_id'].values.tolist()
    Food_id = [str(i) for i in temp_Food_id]
    True_energy = df_food['Energy'].values.tolist()
    predicted_energy = df_food['PrdictEnergy'].values.tolist()
    Food_energy(Food_id, True_energy, predicted_energy)

if __name__ == '__main__' :
    main()
