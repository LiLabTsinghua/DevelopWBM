import matplotlib.pyplot as plt
from proplot import rc
import pandas as pd
from collections import Counter
from pylab import *

# rc["font.family"] = "Helvetica"

def figure_barh(dict_GPT):
    Sub_info = [i for i, v in dict_GPT.items()]
    number = [v for i, v in dict_GPT.items()]
    rc["font.family"] = "Arial"
    plt.rcParams['pdf.fonttype']= 42

    fig, ax = plt.subplots(figsize=(1.5, 2.5))

    fig.patch.set_alpha(0)

    ax.grid(False)
    ax = plt.gca()
    # plt.legend(loc='upper right', fontsize=6, frameon=True, fancybox=True, framealpha=0.2, borderpad=0.3,
    #            ncol=1, markerfirst=True, markerscale=1, numpoints=1, handlelength=3.5)

    # for g in range(len(number)):
    #     ax.barh(gem_info[g], number[g], color=colors[g], label)

    for i in range(len(number)):
        # ax.barh(Sub_info[i], number[i], 0.7, color=colors[i], edgecolor = 'black')
        ax.barh(Sub_info[i], number[i], 0.7, color='#2c7bb6', edgecolor = 'black')

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    for i in ['bottom', 'top', 'right', 'left']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_visible('black')
        ax.spines[i].set_linewidth(0.5) 

    # plt.legend(loc='best', fontsize=6, bbox_to_anchor=(0.7,1))

    my_x_ticks = np.arange(0, 250, 50)
    plt.xticks(my_x_ticks)
    # plt.yticks(my_y_ticks)

    plt.xticks(fontsize=6, rotation=30, ha='right')
    plt.yticks(fontsize=6)

    minorticks_off()

    # plt.savefig("../../Figures/Fig_2a.pdf", dpi=400, bbox_inches = 'tight')
    plt.show()

def data_sort(df):
    pending_Yes = 0
    pending_No = 0
    manual_Right = 0
    manual_TBD = 0
    dict_GPT = {}
    for i in df.index:
        if df.loc[i, 'Correlation'] == 'Pending(Yes)':
            pending_Yes += 1
        else:
            pending_No += 1
            if df.loc[i, 'Manual check'] == 1:
                manual_Right += 1
                sub = df.loc[i, 'Subsystem']
                if sub not in dict_GPT:
                    dict_GPT[sub] = 1
                else:
                    dict_GPT[sub] += 1
            else:
                manual_TBD += 1

    print('There checked {} rxn-gene pairs in Human-GEM.'.format(len(df.index)))
    print('{} were checked the rxn can be associated linked to the gene by GPT.'.format(pending_Yes))
    print('{} were checked the rxn can nor be associated linked to the gene by GPT.'.format(pending_No))
    print('Of the items with no correlation, {} were correct by manual check.'.format(manual_Right))
    print('Of the items with no correlation, {} were not true by manual check.'.format(manual_TBD))
    return dict_GPT

def main():

    df = pd.read_csv('../results/'+'GPT_check.tsv', sep = '\t')
    dict_GPT = data_sort(df)
    top_15 = dict(sorted(dict_GPT.items(), key=lambda item: item[1], reverse=True)[:15])
    figure_barh(top_15)

if __name__ == '__main__' :
    main()