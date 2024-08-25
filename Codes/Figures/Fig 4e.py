import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from proplot import rc
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt

model1 = ['Fetus', 'Adult male', 'Adult female', 'Elderly male', 'Elderly female']
model2 = ['Fetus', 'Adult male', 'Adult female', 'Elderly male', 'Elderly female']

harvest = np.array([[1.0, 0.737, 0.6809, 0.7108, 0.6714],
                    [0.737, 1.0, 0.7126, 0.9218, 0.7054],
                    [0.6809, 0.7126, 1.0, 0.6967, 0.9017],
                    [0.7108, 0.9218, 0.6967, 1.0, 0.7009],
                    [0.6714, 0.7054, 0.9017, 0.7009, 1.0]])

# plt.figure(figsize=(0.5, 0.4))
fig, ax = plt.subplots(figsize=(1.7, 1.6))
im = ax.imshow(harvest, cmap='OrRd')
cbar = fig.colorbar(im, fraction=0.044, pad=0.05)
# 'pad': 调整colorbar和热图之间的距离
# 'fraction': 调整colorbar的大小

# plt.xticks(np.arange(len(model1)), labels=model1, fontsize=7, 
#                      rotation=30, rotation_mode="anchor", ha="right", fontweight = 'bold')
# plt.yticks(np.arange(len(model2)), labels=model2, fontsize=7, fontweight = 'bold')  
plt.xticks(np.arange(len(model1)), labels=model1, fontsize=7, 
                     rotation=30, rotation_mode="anchor", ha="right")
plt.yticks(np.arange(len(model2)), labels=model2, fontsize=7)     
# plt.title("Harvest of local farmers (in tons/year)")

for i in range(len(model1)):
    for j in range(len(model1)):
        if i == j:
            text = plt.text(j, i, "%.2f" % harvest[i, j], ha="center", va="center", color="w", fontsize=7, fontweight = 'bold')
        else:
            text = plt.text(j, i, "%.2f" % harvest[i, j], ha="center", va="center", color="black", fontsize=7, fontweight = 'bold')

plt.grid(False) # 无网格


# 设置边框颜色为黑色
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

plt.savefig("D:/All_Human_GTEx/Papers_code_and_data/Figures/Fig4e_hotmap.pdf", dpi=400, bbox_inches = 'tight')
plt.show()
