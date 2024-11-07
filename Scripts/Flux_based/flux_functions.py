# Extracting GPR's from a model


import cobra
import pandas as pd
import os

def extract_gprs(model):
    """
    Extracts GPRs from a model

    :param model: cobra.Model
    :return: pd.DataFrame
    """
    gpr_list = []
    for reaction in model.reactions:
        gpr = reaction.gene_reaction_rule
        if gpr:
            gpr_list.append([reaction.id, gpr])
    df = pd.DataFrame(gpr_list, columns=['reaction', 'gpr'])
    # save to file here: /Users/nfarooqi/AntiGEM/Data/intermediate/Model_GPRs
    df.to_csv(('/Users/nfarooqi/AntiGEM/Data/intermediate/Model_GPRs/Model_GPRs.csv'), index=False)
    return df

model = cobra.io.read_sbml_model('/Users/nfarooqi/AntiGEM/Models/subtilis.xml')

extract_gprs(model)



def normalize_columns(df):
    """
    Normalize the columns of gene abundance data to be between 0 and 1
    """
    for column in df.columns:
        # if column is numeric, normalize it
        if df[column].dtype in ['int64', 'float64']:
            df[column] = (df[column] - df[column].min()) / (df[column].max() - df[column].min())
    return df

