import shap
import numpy as np
import pandas as pd 
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from matplotlib import rc
from matplotlib import pyplot as plt

from proplot import rc
from pylab import *


def ABS_SHAP(df_shap,df):
    #import matplotlib as plt
    # Make a copy of the input data
    shap_v = pd.DataFrame(df_shap)
    feature_list = df.columns
    shap_v.columns = feature_list
    df_v = df.copy().reset_index().drop('index',axis=1)
    
    # Determine the correlation in order to plot with different colors
    corr_list = list()
    for i in feature_list:
        b = np.corrcoef(shap_v[i],df_v[i])[1][0]
        corr_list.append(b)
    corr_df = pd.concat([pd.Series(feature_list),pd.Series(corr_list)],axis=1).fillna(0)
    # Make a data frame. Column 1 is the feature, and Column 2 is the correlation coefficient
    corr_df.columns  = ['Variable','Corr']
    corr_df['Sign'] = np.where(corr_df['Corr']>0,'red','blue')

    shap_abs = np.abs(shap_v)
    k=pd.DataFrame(shap_abs.mean()).reset_index()
    k.columns = ['Variable','SHAP_abs']
    k2 = k.merge(corr_df,left_on = 'Variable',right_on='Variable',how='inner')
    k2 = k2.sort_values(by='SHAP_abs',ascending = True)

    return k2


def polt_SHAP(df):

    number = df['SHAP_abs'].values.tolist()
    gem_info = df['Variable'].values.tolist()
    Corr = df['Corr'].values.tolist()

    fig, ax = plt.subplots(figsize=(1.5, 3))
    
    plt.rcParams.update({'font.size': 6})
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype']= 42

    ax.grid(False)
    ax = plt.gca()
    # plt.legend(loc='upper right', fontsize=6, frameon=True, fancybox=True, framealpha=0.2, borderpad=0.3,
    #            ncol=1, markerfirst=True, markerscale=1, numpoints=1, handlelength=3.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)

    # for g in range(len(number)):
    #     ax.barh(gem_info[g], number[g], color=colors[g], label)
    bottom = 2

    for i in range(len(number)):
        cor = Corr[i]
        if cor > 0:
            ax.barh(gem_info[i], number[i], 0.7, color='#d73027', edgecolor = 'black')
        else:
            ax.barh(gem_info[i], number[i], 0.7, color='#2c7bb6', edgecolor = 'black')

    for i in ['bottom', 'top', 'right', 'left']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_visible('black')
        ax.spines[i].set_linewidth(0.5) 

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    # plt.xticks(fontsize=7, rotation=30, ha='right')
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)

    plt.xlabel('SHAP Value')  # X轴标题
    plt.ylabel('Variable')     # Y轴标题

    minorticks_off()

    # plt.savefig("../../results/Figures/Fig_SHAP.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()


def main():
    df_SF = pd.read_csv('../../results/Result_BMR_Paramter_male.tsv', sep = '\t')

    all_data = df_SF[df_SF['BMR'] != 0]

    feature_names = all_data.columns.values.tolist()[1:]
    feature_data = pd.DataFrame(all_data,columns=feature_names)

    labels = all_data.BMR

    X_train, X_test, y_train, y_test = train_test_split(feature_data, labels, test_size=0.20, random_state=42)
    # print(X_test)

    rf = RandomForestRegressor(n_estimators=10)
    rf.fit(X_train, y_train)

    explainer = shap.TreeExplainer(rf)
    shap_values = explainer.shap_values(X_train)
    # print(shap_values)
    corr_df = ABS_SHAP(shap_values,X_train)
    polt_SHAP(corr_df)

if __name__== "__main__":
    main()