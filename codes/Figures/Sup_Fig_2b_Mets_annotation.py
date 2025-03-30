import numpy as np
import matplotlib.pyplot as plt

def main():
    data = {
        "MetaNetX": [0.81, 0.72],
        "PubChem": [0.61, 0.44],
        "ChEBI": [0.53, 0.4],
        "KEGG": [0.5, 0.45],
        "HMDB": [0.48, 0.27],
        "LipidMaps": [0.27, 0.17],
        "ModelSEED": [0.52, 0],
    }
    labels = ["Human2", "Human1"]

    categories = list(data.keys())
    values = np.array(list(data.values()))

    num_vars = len(categories)

    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

    values = np.append(values, [values[0]], axis=0)
    angles += angles[:1]


    colors = ['#2c7bb6', '#d7191c']
    plt.figure(figsize=(2.8, 2.8))
    ax = plt.subplot(111, polar=True)

    plt.rcParams.update({'font.size': 8})
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype']= 42


    for i in range(values.shape[1]):
        # ax.fill(angles, values[:, i], color=colors[i], alpha=0.3, label = labels[i])
        ax.fill(angles, values[:, i], color=colors[i], alpha=0.6, label = labels[i])
        ax.plot(angles, values[:, i], color=colors[i], linewidth=0.5)
        # ax.scatter(angles, values[:, i], color=colors[i], edgecolors='black', s=10, alpha=0.6)
        ax.scatter(angles, values[:, i], color=colors[i], edgecolors='black', s=10, alpha=0.6)

    ax.spines['polar'].set_color('black')
    ax.spines['polar'].set_linewidth(1)

    plt.tick_params(axis='x', which='both', pad=10)

    ax.yaxis.grid(True, color='#d9d9d9', linestyle='dashed', linewidth=0.8, alpha=0.7)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.set_yticklabels(["0.1", "0.3", "0.5", "0.7", "0.9"])

    # ax.legend(loc='upper right')

    plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1.35), fontsize=8, frameon=False)
    plt.tight_layout()

    # plt.savefig("../../results/Figures/Fig_MetsAnnotation.pdf", transparent = True, dpi=400, bbox_inches = 'tight')
    plt.show()

if __name__ == '__main__':
    main()
