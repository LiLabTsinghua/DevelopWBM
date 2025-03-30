# Introduction
The `DevelopWBM` repo contains all data, codes and models analyzed in "Human2: Improved reconstruction of human metabolism through combining a community effort with large language models"

# Codes
`./codes/` contains all the codes for calculations and analysis in the article. These codes mainly include:
- Codes for evaluation of Human-GEM;
- Codes for cell- or organ-specific models generation;
- Codes for comparation and analysis of context-specific models;
- Codes for whole-body models reconstruction, analysis and simulation;

For details, please refer to ["./codes/readme.md"](https://github.com/LiLabTsinghua/DevelpWBM/tree/main/codes)

# Models
- The generic genome-scale metabolic model of _Homo sapiens_ is [Human2](https://github.com/SysBioChalmers/Human-GEM). If you have any doubts about the Human2, don't hesitate to ask your questions directly in [here](https://github.com/SysBioChalmers/Human-GEM/issues).

- The organ-specific GEMs are stored in ['./models/'](https://github.com/LiLabTsinghua/DevelpWBM/tree/main/models). These organ-specific GEMs included a total of 24 organs and 4 groups (`Adult male`, `Adult Female`, `Elderly male`, `Elderly female`). Overview of these GEMs(genes/reactions/metabolites):
| **Organ**      | **AdultMale**  | **ElderlyMale**| **AdultFemale**| **ElderlyFemale**| **Fetus**    |
|----------------|----------------|----------------|----------------|----------------|----------------|
| Adipocytes     | 2259/7992/5718 | 2270/8015/5726 | 2262/7978/5712 | 2265/7997/5728 | 2258/7972/5713 |
| Agland         | 2243/7966/5709 | 2239/7962/5706 | 2243/7975/5716 | 2250/8000/5728 | 2243/7924/5694 |
| Brain          | 2270/7830/5618 | 2267/7799/5591 | 2272/7844/5616 | 2271/7832/5616 | 2208/7986/5657 |
| Breast         | NA             | NA             | 2272/8125/5771 | 2270/8118/5777 | NA             |
| Cervix         | NA             | NA             | 2281/8159/5766 | 2268/7919/5670 | NA             |
| Colon          | 2261/8061/5755 | 2268/8021/5741 | 2266/8104/5778 | 2268/8080/5756 | 2256/8062/5760 |
| Heart          | 2214/7519/5494 | 2220/7528/5492 | 2216/7532/5501 | 2220/7540/5508 | 2206/8012/5723 |
| Kidney         | 2293/8234/5790 | 2295/8229/5796 | 2283/8240/5788 | 2299/8292/5821 | 2217/8136/5762 |
| Liver          | 2274/8074/5755 | 2264/8061/5747 | 2278/8138/5759 | 2273/8060/5750 | 2241/8615/5849 |
| Lung           | 2270/8011/5719 | 2275/8047/5740 | 2276/8025/5730 | 2273/8044/5745 | 2218/8068/5763 |
| Muscle         | 2193/7266/5334 | 2204/7355/5372 | 2194/7299/5354 | 2208/7372/5383 | 2209/7919/5704 |
| Ovary          | NA             | NA             | 2274/7934/5689 | 2258/7894/5662 | NA             |
| Pancreas       | 2234/7784/5648 | 2227/7789/5644 | 2228/7741/5629 | 2224/7795/5639 | 2239/8024/5736 |
| Prostate       | 2275/8060/5772 | 2275/8058/5766 | NA             | NA             | NA             |
| Scord          | 2265/7862/5633 | 2260/7849/5634 | 2260/7773/5617 | 2261/7848/5638 | 2256/7786/5600 |
| sIEC           | 2275/8273/5783 | 2294/8269/5787 | 2279/8271/5782 | 2272/8266/5788 | 2241/8302/5819 |
| Skin           | 2272/8025/5728 | 2276/8044/5740 | 2272/8003/5722 | 2274/8049/5733 | 2268/7999/5707 |
| Spleen         | 2254/7907/5684 | 2262/7912/5688 | 2255/7911/5684 | 2254/7930/5701 | 2223/8096/5773 |
| Stomach        | 2259/8047/5761 | 2258/8006/5745 | 2251/8026/5747 | 2252/8017/5754 | 2229/8177/5809 |
| Testis         | 2283/8021/5718 | 2285/8034/5722 | NA             | NA             | NA             |
| Thyroidgland   | 2263/7996/5729 | 2267/8043/5753 | 2263/7988/5732 | 2265/8023/5737 | 2263/7992/5726 |
| Urinarybladder | 2274/8082/5771 | 2236/7756/5588 | 2270/8023/5736 | 2249/7853/5651 | 2267/7990/5726 |
| Uterus         | NA             | NA             | 2262/7951/5711 | 2268/7923/5681 | NA             |

- The whole-body models of are stored in ['./models/WBMs/'](https://github.com/LiLabTsinghua/DevelpWBM/tree/main/models/WBMs/). Overview of these WBMs(genes/reactions/metabolites):
|             |  Fetus    | AdultMale   | AdultFemale | ElderlyMale   | ElderlyFemale |
|-------------|-----------|-------------|-------------|---------------|--------|
| Genes       | 2506      | 2383        | 2384        | 2386          | 2382   |
| Metabolites | 105829    | 104823      | 113880      | 103615        | 113272 |
| Reactions   | 133609    | 125360      | 135322      | 123730        | 135139 |

# Citation
- Please cite the paper ["Human2: Improved reconstruction of human metabolism through combining a community effort with large language models"](links)

# Contact
- Jiahao Luo([@Jiahao](https://github.com/JHL-452b)), Institute of Biopharmaceutical and Health Engineering, Tsinghua Shenzhen International Graduate School, Tsinghua University, Shenzhen, China
- Feiran Li([@feiran](https://github.com/feiranl)), Institute of Biopharmaceutical and Health Engineering, Tsinghua Shenzhen International Graduate School, Tsinghua University, Shenzhen, China
