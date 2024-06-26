import os
import pandas as pd
from rdkit import Chem
from thermo import functional_groups
from rdkit.Chem import rdMolDescriptors, Descriptors, Lipinski, Crippen, inchi
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

### STEP 6 - RDKIT PARSE FILTER ##############################################

def rdkit_parse(df, smi_col, id_col):
    """
    Check if SMILES strings in a Dataframe parsable by RDKit and chemistry problems detected by RDkit
    Input: 
        df: Pandas dataframe
        smi_col (str): name of column containing smi
        id_col (str): name of column containing IDs (e.g, CIDs)
    Output: 
        problem_dict (dict): problems detected by RDKit (id -> problem)
        cannot_parse (list): list of IDs that non-parsable by RDKit
    """
    problem_dict = {}
    cannot_parse = []

    for i in df[id_col]:

        # Convert each smi to molecule:
        smi = df.loc[df[id_col] == i, smi_col].values[0]
        m = Chem.MolFromsmi(smi, sanitize=True)

        if m is None:
            cannot_parse.append(i)  # Save if molecule is non-parsable
        else:
            problem = Chem.DetectChemistryProblems(m)  # Identify and capture error messages when creating mol objects.
            if problem:
                problem_dict[i] = problem
    
    if problem_dict:
        print(problem_dict)
    else:
        print("No problems detected")

    return problem_dict, cannot_parse

### STEP 7 - INORGANIC FILTER ##############################################

def inorganic_filter(df, chem_rep_col, id_col, type):
    """
    Check if compounds in a Dataframe are organic or inorganic.
    Input:
        df: Pandas dataframe
        chem_rep_col (str): name of the smi or InChI col that contains the chemical representations
        id_col (str): name of the column that contain IDs (e.g, CIDs)
        type (str): either 'smi' or 'inchi'
    Output: 
        inorganic_ids (list): IDs of compounds that are inorganic
        organic_ids (list): IDs of compounds that are organic
    """
    inorganic_ids = []
    organic_ids = []

    for index, row in df.iterrows():
        id = row[id_col]
        chem_rep = row[chem_rep_col]

        if type == 'smi':
            mol = Chem.MolFromsmi(chem_rep, sanitize=True)
        elif type == 'inchi':
            mol = Chem.MolFromInchi(chem_rep, sanitize=True)
        else:
            print('Filter for this type of chemical representation not available.')
            continue

        if functional_groups.is_inorganic(mol):
            inorganic_ids.append(id)
        else:
            organic_ids.append(id)
    
    return inorganic_ids, organic_ids

### STEP 8 - MIXTURE FILTER ##############################################

def quick_check_mixtures(name, smi_list):
    """
    Fast counter for mixtures present in a list of SMILES string. 
    """
    count=0
    for smi in smi_list:
        if '.' in smi:
            count+=1
    print(f"Total number of mixtures in {name} is {count}")

def process_smi_mixtures(df, smi_col, id_col):
    """
    From a given Dataframe, detect and handle mixtures based on SMILES representation.
    Input: 
        df: Pandas dataframe
        smi_col (str): name of column containing smi
        id_col (str): name of column containing IDs (e.g, CIDs)
    Output: 
        (dictionary: ID -> smi)
            processed: non-mixture forms of every smi in the given dataset
            removed: mixture components or mixtures removed from the original mixture smi
            small_organic: small organic molecules removed from a size-imbalanced mixtures of organic molecules
            small_inorganic: small inorganic molecules removed from a mixture of both organic and inorganic molecules
            big_organic_not_lipinski: large organic molecules not kept due to not passing the lipinski criteria
            cleaned_from_mixtures: non-mixture forms after handling of the orignal mixtures
    """
    processed = {}
    removed = {}
    small_inorganic = {}
    small_organic = {}
    big_organic_not_lipinski = {}
    cleaned_from_mixtures = {}
    
    for index, row in df.iterrows():
        id = row[id_col]
        smi = row[smi_col]
        
        if '.' in smi: # Check for mixtures
            molecules = smi.split('.')
            mols = [Chem.MolFromSmiles(mol) for mol in molecules]
            num_atoms = [mol.GetNumAtoms() for mol in mols if mol is not None]
            
            if all(x == num_atoms[0] for x in num_atoms): # Check if all molecules have the same number of atoms. If yes, keep one of them
                processed[id] = molecules[0]
                cleaned_from_mixtures[id] = molecules[0]
                removed[id] = '.'.join(molecules[1:])
            else:
                if max(num_atoms) - min(num_atoms) <= 5:
                    removed[id] = smi # Remove the mixture if the difference in number of atoms is less than 5
                    print(f"Cannot decide between {molecules} for CID {id}")
                else:
                    max_index = num_atoms.index(max(num_atoms))
                    min_index = num_atoms.index(min(num_atoms))
                    if functional_groups.is_inorganic(mols[min_index]) == True:
                        processed[id] = molecules[max_index]
                        cleaned_from_mixtures[id] = molecules[max_index]
                        removed[id] = molecules[min_index] # Keep the organic molecule and remove the inorganic one
                        small_inorganic[id] = molecules[min_index]
                    else:
                        big_molecule = mols[max_index]
                        
                        # Calculate properties for Lipinski's rule of five
                        mw = Descriptors.MolWt(big_molecule)
                        hbd = rdMolDescriptors.CalcNumHBD(big_molecule)
                        hba = rdMolDescriptors.CalcNumHBA(big_molecule)
                        logp = Crippen.MolLogP(big_molecule)

                        # Check Lipinski's criteria
                        if mw <= 500 and hbd <= 5 and hba <= 10 and logp <= 5:
                            processed[id] = molecules[max_index] # Keep the big organic molecule if it passes Lipinski's rule of five
                            cleaned_from_mixtures[id] = molecules[max_index]
                            removed[id] = '.'.join([molecules[i] for i in range(len(molecules)) if i != max_index])
                            small_organic[id] = molecules[min_index]
                        else:
                            removed[id] = smi # Remove the mixture if the big organic molecule does not pass Lipinski's rule of five
                            big_organic_not_lipinski[id] = molecules[max_index]
                            print(f"Big organic molecule for CID {id} does not pass Lipinski's rule of five")

        else:
            processed[id] = smi
            
    return processed, removed, small_organic, small_inorganic, big_organic_not_lipinski, cleaned_from_mixtures

def process_mixture_df(name, df, processed, removed, small_organic, small_inorganic, smi_col, id_col):
    """
    Update a given Dataframe with information on mixture handling. 
    Input: 
        Name: name for the dataset
        df: Pandas Dataframe with at least SMILES and ID columns 
        processed, removed, small_organic, small_inorganic: dictionaries returned by process_smi_mixtures()
        smi_col (str): name of column containing smi
        id_col (str): name of column containing IDs (e.g, CIDs)
    Output:
        new_df: Dataframe with SMILES mixtures processed and additional columns containing mixture handling information
    """
    indices_to_drop = []  # List to keep track of row indices that should be dropped
    
    for index, row in df.iterrows():
        id = row[id_col]
        if id in processed:
            df.loc[index, smi_col] = processed[id]
            if id in removed:
                df.loc[index, 'Mol removed from mixture'] = removed[id]
            if id in small_organic:
                df.loc[index, 'Small organic molecule'] = small_organic[id]
            if id in small_inorganic:
                df.loc[index, 'Small inorganic molecule'] = small_inorganic[id]
        else:
            indices_to_drop.append(index)
            print(f"{id} has been removed from {name} because it is a mixture with less than 5 atoms difference or the big organic molecule does not pass Lipinski's rule of five.")
    
    # Drop rows outside the loop and reset index if needed
    new_df = df.drop(indices_to_drop).reset_index(drop=True)
    
    return new_df

### STEP 11.2 - AUTOFLUORESCENCE FILTER ##############################################

def load_autofluorescence_cids(data_folder):
    """
    Return a list of CIDs that were confirmed to be autofluorescence interference compounds, from a series of HTS experiments. 
    """
    autofluorescence_aids = ['587', '588', '590', '591', '592', '593', '594']
    autofluorescence_cids = []
    col_list = ['PUBCHEM_CID', 'PUBCHEM_ACTIVITY_OUTCOME']

    count = 0
    autofluorescence_path = os.path.join(data_folder, 'autofluorescence', 'autofluorescence_cids.txt')

    # Check if autofluorescence_cids.txt is available, if not download it
    if not os.path.exists(autofluorescence_path):
        os.makedirs(os.path.dirname(autofluorescence_path), exist_ok=True)
        for AID in tqdm(autofluorescence_aids, desc="Downloading autofluorescence data"):
            url = f'https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?query=download&record_type=datatable&actvty=all&response_type=save&aid={AID}'
            autofluorescence_df = pd.read_csv(url, usecols=col_list)
            # Delete rows with NaN values
            autofluorescence_df = autofluorescence_df.dropna(subset=['PUBCHEM_CID', 'PUBCHEM_ACTIVITY_OUTCOME'])

            # Convert CIDs to int
            autofluorescence_df['PUBCHEM_CID'] = autofluorescence_df['PUBCHEM_CID'].astype(int)

            # Keep only rows that are "Active"
            autofluorescence_df = autofluorescence_df[autofluorescence_df['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active']
            autofluorescence_cids.extend(autofluorescence_df['PUBCHEM_CID'].tolist())

        # Drop duplicates in the list
        autofluorescence_cids = list(set(autofluorescence_cids))
        
        # Export as txt file
        with open(autofluorescence_path, 'w') as file:
            for cid in autofluorescence_cids:
                file.write(f"{cid}\n")
    else:
        # If file exists, read the CIDs from the file
        with open(autofluorescence_path, 'r') as file:
            autofluorescence_cids = [int(line.strip()) for line in file]
    print('Autofluorescence compounds imported successfully!')

    return autofluorescence_cids

### STEP 11.3 - RDKIT PAIN FILTER ##############################################

def detect_pains(df, smi_col, id_col):
    """
    Use the FilterCatalogs module from RDkit to detect PAIN molecules based on substructure matching. 
    Input: 
        Pandas Dataframe 
        smi_col: name of column containing SMILES strings
        id_col: name of column containing identifiers (e.g, CIDs)
    Output: 
        list of ids corresponding to PAINs
    """
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    pains = []
    count_pains = 0
    count_not_pains = 0
    for i in df.index:
        smile = str(df.loc[i, smi_col])
        id = df.loc[i, id_col]
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            if catalog.HasMatch(mol):
                pains.append(id)
                count_pains += 1
                print(f'{count_pains} pains detected')
            else: 
                count_not_pains += 1
    return pains

### STEP 12 - DRUG-LIKENESS FILTER ##############################################

def drug_likeness_filter(smiles):
    """
    Check if a given SMILES satisfies the common standard conditions for drug-likeness. 
    Input:
        SMILES (str)
    Output: 
        Result (bool)
    """

    # Convert SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False  # Return False if the molecule cannot be parsed
    
    # Check molecular weight
    mw = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if not (150 < mw < 800):
        return False
    
    # Check AlogP
    logp = Chem.Crippen.MolLogP(mol)
    if not (-0.3 < logp < 5):
        return False
    
    # Check number of rotatable bonds
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    if rotatable_bonds >= 15:
        return False
    
    # Check H-bond acceptor count and H-bond donor count
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)
    if hba >= 15 or hbd >= 15:
        return False
    
    # Check total formal charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if not (-2 < total_charge < 2):
        return False
    
    # If all filters passed, return True
    return True

def drug_likeness_filter_multiprocessing(df, smi_col, id_col):
    """
    Update a given dataframe by dropping molecules that don't pass the drug-likeness filter.
    Output: 
        List of IDs corresponding to molecules that didn't pass. 
    """
    to_drop = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        # Create future tasks for each SMILES string in the dataframe
        futures = {executor.submit(drug_likeness_filter, row[smi_col]): row[id_col] for index, row in df.iterrows()}
        
        # Use tqdm to display progress bar
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing SMILES"):
            id = futures[future]
            if not future.result():
                to_drop.append(id)
    return to_drop
