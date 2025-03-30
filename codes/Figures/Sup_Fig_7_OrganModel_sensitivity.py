import pandas as pd
import numpy as np
import xlrd
import re
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.io import load_yaml_model, save_yaml_model, load_matlab_model, save_matlab_model
from tqdm import tqdm
import os
import csv

def Sensitivity_Process(Gender, model):
    # filepath_1 = 'D:\\All_Human_GTEx\\Human2_Paper_and_Suplementary\\New_Human_GEM\\results\\Sensitivity\\'
    # filepath = filepath_1+Gender+'_sensitivity\\'
    filepath = '../../results/Sensitivity/'+Gender+'_sensitivity/'
    fileNames = os.listdir(filepath)
    fileNames = [f.replace('_sensitivity.tsv', '') for f in fileNames]

    subsystems = [r.subsystem[0] for r in model.reactions]
    subsystems = list(set(subsystems))
    df_1 = pd.DataFrame(0, index=subsystems, columns=fileNames)
    df_2 = pd.DataFrame(0, index=subsystems, columns=fileNames)
  
    sub_num = []
    for i in subsystems:
        for j in fileNames:
            df = pd.read_csv(filepath+j+'_sensitivity.tsv', sep = '\t')
            dict_rxn_corr = dict(zip(df['Rxn_id'].values.tolist(), df['Values'].values.tolist()))

            temp_cor = []
            for k, v in dict_rxn_corr.items():
                if pd.isna(v):
                    # print("The variable is NaN")
                    pass
                else:
                    if k in model.reactions:
                        if i == model.reactions.get_by_id(k).subsystem[0]:
                            temp_cor.append(v)

            # sub_num.append(len(temp_cor))
            df_1.loc[i,j] = np.median(temp_cor)
            df_2.loc[i,j] = len(temp_cor)


    max_values = df_2.max(axis=1)
    filtered_subs = max_values[max_values <= 10].index # Subsystems with number of reaction less than 10 are not considered
    filtered_subs = filtered_subs.tolist()
    df_1.drop(filtered_subs, inplace=True)

    df_1 = df_1.transpose()

    cols_to_drop = df_1.columns[(df_1 > 0.95).all()] # Subsystems with all values ​​greater than 0.95 are not considered
    df_1 = df_1.drop(columns=cols_to_drop)

    return df_1


def main():
    filepath = '../../results/Sensitivity/'
    Gender = 'Female'
    
    model = load_matlab_model('../../models/Human2.mat')
    df_1 = Sensitivity_Process(Gender, model)
    df_1.to_csv(filepath+Gender+'_sensitivity_results_0330.tsv', sep='\t')
    
    # You can run the following code in R to visualize the results:
    # library(pheatmap)
    # Sex = 'Male' # or Female
    # # Input_file = paste(Sex, '_sensitivity_results_10_remove_0.95.tsv', sep="")
    # Input_file = paste(Sex, '_sensitivity_results_0330.tsv', sep="")
    # Output_file = paste(Sex, '_sensitivity_results_10_remove_0.95', '_clustered_data.tsv',  sep="")

    # # load data
    # data <- read.delim(Input_file, row.names = 1)
    # data[is.na(data)] <- 0

    # # Generate heatmap and return values
    # # pheatmap_result <- pheatmap(data, fontsize = 7, width = 8, height = 5)
    # pheatmap_result <- pheatmap(data, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", fontsize = 7, width = 7, height = 4, angle_col = 45)


    # # get row names
    # clustered_row_names <- rownames(data)[pheatmap_result$tree_row$order]

    # # get column names
    # clustered_col_names <- colnames(data)[pheatmap_result$tree_col$order]

    # # get ranked data
    # clustered_data <- data[clustered_row_names, clustered_col_names]

    # # store data
    # write.table(clustered_data, file = Output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
    
if __name__== "__main__":
    main()