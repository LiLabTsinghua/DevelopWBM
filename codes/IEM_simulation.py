import pandas as pd
import numpy as np
import cobra
import re
from cobra import Model, Reaction, Metabolite
from cobra.io import load_yaml_model, save_yaml_model
from cobra.io import load_matlab_model, save_matlab_model
from cobra.io import read_sbml_model
from cobra.io import load_json_model
from tqdm import tqdm
import os
import csv

def load_metabolic_model(modelfile, model_id):
    # Converts the file name to lowercase
    model_name = model_id.split('.')[0]
    if model_name.startswith('ec'):
        model_name = model_name.replace('ec', '')
    else:
        model_name = model_name

    file_extension = model_id.lower().split('.')[-1]
    
    if file_extension == 'xml' or file_extension == 'sbml':
        try:
            model = read_sbml_model(modelfile+model_id)
            return model, model_name
        except:
            print("Something error with the model!")
    elif file_extension == 'json':
        try:
            model = load_json_model(modelfile+model_id)
            return model, model_name
        except:
            print("Something error with the model!")
    elif file_extension == 'yml':
        try:
            model = load_yaml_model(modelfile+model_id)
            return model, model_name
        except:
            print("Something error with the model!")
    elif file_extension == 'mat':
        try:
            model = load_matlab_model(modelfile+model_id)
            return model, model_name
        except:
            print("Something error with the model!")
    else:
        return "Unsupported file format!"
    

def effected_rxns(gene_list, model):
    """
    The affected reactions were defined as the reaction in which 
    the reaction boundary changed after the corresponding gene 
    was knocked out
    """

    dict_original_rxn_lb = {}
    dict_original_rxn_ub = {}
    for r in model.reactions:
        dict_original_rxn_lb[r.id] = r.lower_bound
        dict_original_rxn_ub[r.id] = r.upper_bound

    effected_rxns_list = []
    for g in gene_list:
        if g in model.genes:
            for r in model.genes.get_by_id(g).reactions:
                effected_rxns_list.append(r.id)
    
    return effected_rxns_list


def Met_in_model(m, model, model_name):
    """
    Check the biomarker metabolites whether in the model. Identifies 
    of metabolites should be got ready and compartment 'c' was given 
    priority, followed by 'e'. 
    """

    if model_name == 'Human2':
        model_m_c = str(m)+'c'
        model_m_e = str(m)+'e'

    elif model_name == 'Human1':
        model_m_c = str(m)+'c'
        model_m_e = str(m)+'s'
    
    elif model_name == 'Recon3D_301':
        model_m_c = str(m)+'[c]'
        model_m_e = str(m)+'[e]'
    else:
        model_m_c = str(m)+'c'
        model_m_e = str(m)+'e'

    if model_m_c in model.metabolites:
        return model_m_c
    elif model_m_e in model.metabolites and model_m_c not in model.metabolites:
        return model_m_e
    else:
        return 'Not_In_Model'
        
    # if model_m_e in model.metabolites:
    #     return model_m_e
    # elif model_m_c in model.metabolites and model_m_e not in model.metabolites:
    #     return model_m_c
    # else:
    #     return 'Not_In_Model'
    

def Gene_del_in_model(g, model):
    """
    Check the error genes whether in the model. Identifies of genes 
    should be got ready. 'g' belong to str, if there are multiple 
    genes in 'g' then need to use ';' to split it out.
    """

    del_gene_list = []
    if ';' in g:
        temp = g.split(';')
        for del_g in temp:
            if del_g in model.genes:
                del_gene_list.append(del_g)
            else:
                pass
    elif g in model.genes and ';' not in g:
        del_gene_list.append(g)

    else:
        pass

    return del_gene_list


def Hams_media(df_media, model, model_name, model_type):
    """
    set Ham's media constraints to models.
    """

    # Colse all exchange rxns firsttly.
    # for ex in model.exchanges:
    #     ex.bounds = (0,1000)
    for r in model.reactions:
        if len(r.metabolites) == 1:
            r.bounds = (0,1000)
    
    media_mets = df_media[model_name].values.tolist()
    # Open media mets
    for m in media_mets:
        # for r in model.exchanges:
        for r in model.reactions:
            if len(r.metabolites) == 1:
                if m in r.reaction:
                    if re.search(r"C\d", str(model.metabolites.get_by_id(m).formula)) and model_type == 'Model':
                    # if re.search('C' in str(model.metabolites.get_by_id(m).formula)) and model_type == 'Model':
                        r.bounds = (-10,1000)
                    elif re.search(r"C\d", str(model.metabolites.get_by_id(m).formula)) and model_type == 'ecModel':
                    # elif re.search('C' in str(model.metabolites.get_by_id(m).formula)) and model_type == 'ecModel':
                        r.bounds = (-1000,1000)
                        # r.bounds = (-10,1000)
                    else:
                        r.bounds = (-1000,1000)
                else:
                    pass

    return model


def Calculation(m_id, model):
    """
    add demand reaction for biomarker metabolite and set it as objective to simulate it.
    """

    if m_id != 'Not_In_Model':
        # model.add_boundary(model.metabolites.get_by_id(m_id), type='demand')
        rxn_id = 'DM_IEM_' + m_id
        reaction = Reaction(rxn_id)
        model.add_reactions([reaction])
        reaction.build_reaction_from_string(m_id+' --> ')
        model.objective = model.reactions.get_by_id(rxn_id)
        sol = model.optimize()
        sol_value = sol.objective_value
        sol_status = sol.status
    else:
        sol_value = 'Met_not_in_Model'
        sol_status = 'Met_not_in_Model'

    return sol_value, sol_status


def health_model_define(gene_list, model):
    """
    In the healthy state, a maximum flux of the affected reactions are required.
    """

    # identify effected rxns
    with model as cmodel:
        uni_effected_rxns = effected_rxns(gene_list, cmodel)
        # print(str(uni_effected_rxns))

    # add pesudo met and rxn into model
    new_met = Metabolite('pseudo_met',
                      formula = 'X',
                      name = 'pseudo_met_health',
                      compartment = 'e', 
                      charge = 0)
    model.add_metabolites(new_met)
    # model.add_boundary(new_met, type='exchange')
    pseudo_rxn = 'EX_IEM_' + new_met.id
    # print(model.reactions.get_by_id(pseudo_rxn).bounds)
    reaction = Reaction(pseudo_rxn)
    model.add_reactions([reaction])
    reaction.build_reaction_from_string('pseudo_met <=> ')

    # add pesudo_met into effected_rxns
    for rxn in uni_effected_rxns:
        model.reactions.get_by_id(rxn).add_metabolites({'pseudo_met': 1})

    try:
        model.objective = model.reactions.get_by_id(pseudo_rxn)
        sol = model.optimize()
        pseudo_rxn_flux = int(sol.objective_value)
        # print(pseudo_rxn_flux)
    except:
        print(str(gene_list)+': '+str(uni_effected_rxns))
        pseudo_rxn_flux = -1000

    model.reactions.get_by_id(pseudo_rxn).lower_bound = pseudo_rxn_flux

    return model


def disease_model_define(gene_list, model):
    """
    In the disease state, the flux of the affected reactions are set to zero.
    """

    # identified effected rxns by gene_list
    with model as cmodel:
        uni_effected_rxns = effected_rxns(gene_list, cmodel)

    # add pesudo met and rxn into model
    new_met = Metabolite('pseudo_met',
                      formula = 'X',
                      name = 'pseudo_met_health',
                      compartment = 'e', 
                      charge = 0)
    model.add_metabolites(new_met)
    # model.add_boundary(new_met, type='exchange')
    pseudo_rxn = 'EX_IEM_' + new_met.id
    # print(model.reactions.get_by_id(pseudo_rxn).bounds)
    reaction = Reaction(pseudo_rxn)
    model.add_reactions([reaction])
    reaction.build_reaction_from_string('pseudo_met <=> ')

    # add pesudo_met into effected_rxns
    for rxn in uni_effected_rxns:
        model.reactions.get_by_id(rxn).add_metabolites({'pseudo_met': 1})

    model.reactions.get_by_id(pseudo_rxn).lower_bound = 0
    model.reactions.get_by_id(pseudo_rxn).upper_bound = 0

    return model


def main():

    folder = '../models/'
    print('Please put your model in '+ folder)

    model_id = input('Please input model name: ')
    print(model_id)
    file_path = os.path.join(folder, model_id)

    if os.path.isfile(folder+model_id):
        # load model
        ihuman, model_name = load_metabolic_model(folder, model_id)
        print(model_name)

        # check model type: "ecModel" or "Model"
        if 'prot_pool_exchange' in ihuman.reactions:
            model_type = 'ecModel'
            model_prefix = 'ec_'
        else:
            model_type = 'Model'
            model_prefix = ''

        # set Ham's media
        df_media = pd.read_csv('../data/'+'Hams_media.tsv', sep = '\t')
        ihuman = Hams_media(df_media, ihuman, model_name, model_type)

        if model_type == 'ecModel':
            ihuman.reactions.get_by_id('prot_pool_exchange').lower_bound = -500
            # ihuman.reactions.get_by_id('prot_pool_exchange').lower_bound = -200
        elif model_type == 'Model':
            if model_name == 'Human2' or model_name == 'Human-GEM':
                ihuman.reactions.get_by_id('MAR13082').bounds = (1,1)
            elif model_name == 'Human1':
                ihuman.reactions.get_by_id('biomass_human').bounds = (1,1)
            elif model_name == 'Recon3D_301':
                ihuman.reactions.get_by_id('biomass_reaction').bounds = (1,1)
            else:
                pass

        # load IEM data
        df = pd.read_csv('../data/'+'IEM_data.tsv', sep = '\t')

        if 'Human' in model_name:
            met_col_name = 'met_'+model_name
            df_result = df[['Gene', met_col_name, 'In_Vivo']]
            df_result['col1'] = np.nan
            df_result['col2'] = np.nan
            new_column_names = {'Gene': 'Gene_id', met_col_name: 'Biomarker', 'In_Vivo': 'Exp_(abnormal/normal)', 'col1': 'Simulation', 'col2': 'Comparation'}
            df_result = df_result.rename(columns=new_column_names)

        elif 'Recon' in model_name:
            met_col_name = 'met_'+model_name
            df_result = df[['Gene_Recon', met_col_name, 'In_Vivo']]
            df_result['col1'] = np.nan
            df_result['col2'] = np.nan
            new_column_names = {'Gene_Recon': 'Gene_id', met_col_name: 'Biomarker', 'In_Vivo': 'Exp_(abnormal/normal)', 'col1': 'Simulation', 'col2': 'Comparation'}
            df_result = df_result.rename(columns=new_column_names)
        
        All_gene_list = df_result['Gene_id'].values.tolist()
        All_biomarker_mets = df_result['Biomarker'].values.tolist()
    
        for i in range(len(All_gene_list)):
            
            del_gene_list = Gene_del_in_model(All_gene_list[i], ihuman)
            model_m_id = Met_in_model(All_biomarker_mets[i], ihuman, model_name)
            
            
            if len(del_gene_list) > 0 and model_m_id != 'Not_In_Model':
                # origin_model simulation
                with ihuman as model_health:
                    model_health = health_model_define(del_gene_list, model_health)
                    sol_normal_value, sol_normal_status = Calculation(model_m_id, model_health)


                with ihuman as model_del:
                    model_del = disease_model_define(del_gene_list, model_del)
                    sol_abnormal_value, sol_abnormal_status = Calculation(model_m_id, model_del)

                if isinstance(sol_normal_value, float) and isinstance(sol_abnormal_value, float) and sol_normal_value > 0:
                    fc = round((sol_abnormal_value-sol_normal_value)/sol_normal_value, 4)
                    # print(fc)
                    df_result.loc[i, 'Simulation'] = fc
                    if fc > 0 and df_result.loc[i, 'Exp_(abnormal/normal)'] == 'increase':
                        df_result.loc[i, 'Comparation'] = 1
                    elif fc < 0 and df_result.loc[i, 'Exp_(abnormal/normal)'] == 'decrease':
                        df_result.loc[i, 'Comparation'] = 1
                    else:
                        df_result.loc[i, 'Comparation'] = 0
                else:
                    fc = 'No_result'
                    df_result.loc[i, 'Simulation'] = fc
                    df_result.loc[i, 'Comparation'] = fc
            else:
                df_result.loc[i, 'Simulation'] = 'Gene or met not in model!'
                df_result.loc[i, 'Comparation'] = 'Gene or met not in model!'
        
        df_result.to_csv('../results/IEM/'+model_prefix+model_name+'_IEM.tsv', sep = '\t', index = False)

        count_ones = (df_result['Comparation'] == 1).sum()
        proportion = round(count_ones / len(df_result.index), 4)
        print("The accuracy of "+model_prefix+model_name+" is: ", proportion)
    
    else:
        print(model_name + " is not in " + folder)

if __name__ == '__main__':
    main()