from rdkit import Chem
from chembl_structure_pipeline import checker
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

### STEP 9 - NEUTRALIZATION ####################################################

def neutralize_atoms(mol):
    """
    Code adapted from https://www.rdkit.org/docs/Cookbook.html. 
    Source: https://baoilleach.blogspot.com/2019/12/no-charge-simple-approach-to.html
    (Noel O’Boyle, 2019)

    Return a neutralized molecule for a given input Mol object. 
    Additional handling was added for molecules with tetracoordinated boron. 
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            
            #Skip adjustment for tetracoordinated boron
            if atom.GetAtomicNum() == 5 and atom.GetDegree() == 4: #ADD COMMENT
                continue  # Just bypass the problematic atom

            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

### STEP 10 - AROMATIZATION ####################################################

def aromatize_smile(mol):
    """
    Dekekulize an input Mol object and return the aromatic form of isomeric SMILES. 
    """
    Chem.Kekulize(mol)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL)
    aromatic_smiles = Chem.MolToSmiles(mol, isomericSmiles = True)
    return aromatic_smiles

### STEP 13 - CHEMBL CURATION PIPELINE CHECKER ####################################################

def checker_score(smiles, id):
    """
    Return the penalty score associated with a given compound id, computed by the ChemBL Curation Pipeline Checker.
    Input: 
        smiles (str)
        id: identifier associated with the compound (e.g, CID)
    """
    result = checker.check_molblock(Chem.MolToMolBlock(Chem.MolFromSmiles(smiles)))
    if result == ():
        penalty_score = 0
    else:
        penalty_score = result[0][0]
    return id, penalty_score

def checker_multiprocessing(df, smi_col, id_col):
    """
    Return a dictionary of penalty score associated with unique ids in a dataframe. 
    """
    chembl_score_dict = {}
    with ThreadPoolExecutor(max_workers=10) as executor:
        # Map futures to CIDs directly for easier reference
        futures = {executor.submit(checker_score, row[smi_col], row[id_col]): row[id_col] for _, row in df.iterrows()}
        # Properly use tqdm to create a progress bar
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing SMILES"):
            cid, penalty_score = future.result()
            chembl_score_dict[cid] = penalty_score
    return chembl_score_dict