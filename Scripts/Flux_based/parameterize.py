import cobra
import pandas as pd
import riptide


def parameterize_model(model, transcriptomic):

    # Set model medium

    medium = model.medium
    medium['EX_glu__L_e'] = 1
    medium['EX_trp__L_e'] = 1
    medium['EX_cit_e'] = 1
    model.medium = medium
    
    # Prepare abundance data

    transcriptomic = transcriptomic.iloc[:, [1, 2]]  
    transcriptomic.columns = ['gene', 'value'] 
    transcriptomic['gene'] = transcriptomic['gene'].astype(str)
    transcriptomic['value'] = pd.to_numeric(transcriptomic['value'], errors='coerce')

    gene_df = pd.read_csv('/Users/yahyafarooqi/Documents/AntiGEM/Data/intermediate/gene_df2.tsv', sep='\t')
    
    filtered_gene_df = gene_df[gene_df['gene'].isin(transcriptomic['gene'])]

    merged_data = pd.merge(filtered_gene_df, transcriptomic, on='gene', how='inner')
    merged_data = merged_data.drop(columns=['gene'])
    merged_data.to_csv('/Users/yahyafarooqi/Documents/AntiGEM/Data/intermediate/riptide_transcripts/{transcriptomic}_updated_abundance.tsv', sep='\t', index=False)

    # Parameterize model

    path = '/Users/yahyafarooqi/Documents/AntiGEM/Data/intermediate/riptide_transcripts/{transcriptomic}_updated_abundance.tsv'

    tr = riptide.read_transcription_file(path, header = True)   
    riptide = riptide.maxfit(model, tr)

    return riptide

model = cobra.io.read_sbml_model('/Users/yahyafarooqi/Documents/AntiGEM/models/subtilis_ iYO844.xml')
transcriptomic = pd.read_csv('/Users/yahyafarooqi/Documents/AntiGEM/Data/intermediate/transcriptomics_data.tsv', sep='\t')


object = parameterize_model(model, transcriptomic)



    











