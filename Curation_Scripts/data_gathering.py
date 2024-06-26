import os
import pandas as pd

def download_and_save(AIDs, path, col_list, smi_col, id_col, activity_col): 
    """
    Retrieve desired PubChem datasets from PubChem and saves them into a destined path. 
    Input: 
        AIDs: list of desired AIDs.
        path 
        col_list: list of desired columns to extract
        smi_col: name of column containing SMILES strings
        id_col: name of column containing identifiers (e.g, CIDs)
        activity_col: name of column containing activity outcomes
    """
    count = 0
    for AID in AIDs:
        url = f'https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?query=download&record_type=datatable&actvty=all&response_type=save&aid={AID}'
        assay = pd.read_csv(url, usecols=col_list)

        #convert SMILES to string
        assay[smi_col] = assay[smi_col].astype(str)

        #delete rows with nan values
        assay = assay.dropna(subset=[smi_col, activity_col, id_col])

        #convert cids to int:
        assay.loc[:, id_col] = assay[id_col].astype(int)

        #reindex assay dataframe from 0:
        assay.reset_index(drop=True, inplace=True)

        #Create a new folder to save the data:
        if not os.path.exists(f'{path}/before_finished/step_1'):
            os.makedirs(f'{path}/before_finished/step_1')

        #save assay data
        assay.to_csv(f'{path}/before_finished/step_1/AID{AID}.csv', index=False)
        count += 1
        print(f'Completed {count} out of {len(AIDs)} datasets.')