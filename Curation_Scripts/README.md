
## Curation Scripts

### Package Requirements

1. RDKit
2. thermo
3. tqdm
4. Biopython
5. chembl_structure_pipeline

### File Description

**Jupiter Notebook files (.ipynb)**: These are the notebook with scripts for curation process. The general pipeline includes 15 main steps:

- Step 0 - Initial Set-up: Import all necessary packages and set up basic variables such as dataset name, .csv file column format, and directory.

- Step 1 - Data Gathering: Retrieve datasets from PubChem based on a predefined list of assay identifiers (AIDs).

- Step 2 - Isomeric SMILES: Update all SMILES strings in all datasets using the isomeric SMILES output option from the PubChem Identifier Exchange Service. This step requires manual interaction with the PubChem server. 

- Step 3 - Import InChIs: Import standard InChIs corresponding to CIDs from the PubChem Identifier Exchange Service. This step also requires manual interaction with the PubChem server.

- Step 4 - Check Duplicates: Identify duplicates in the datasets based on CIDs, SMILES, and InChIs formats. Examine redundancy in each type of chemical identifier (e.g., compounds having the same CIDs but different SMILES).

- Step 5 - Hierarchical Curation: Follow the hierarchical curation process described in the associated paper's supplement.

- Step 6 - RDKit Parse Check: Ensure all compounds are parsable by RDKit.

  

- Step 7 - Inorganics Filter: Filter out all inorganic molecules from the datasets.

- Step 8 - Mixture Handling: Screen through all SMILES strings in the datasets and remove or update mixture compounds according to predefined rules.

- Steps 9 & 10 - Neutralization & Aromatization of Molecules: Neutralize all charged molecules and aromatize the SMILES strings. Perform these steps based on user needs.

- Step 11 - PAIN Filters: Apply a series of three filters for Pan-Assay Interference Compounds: Promiscuity Filter (or Frequency of Hits Filter), Autofluorescence Filter, and a comprehensive PAIN pattern matching filter.

- Step 12 - Drug-likeness Filter: Keep only molecules that exhibit known molecular property ranges for drugs.

- Step 13 - ChemBL Curation Pipeline Checker: Use a checker to identify problematic molecules for further manual inspection before exporting the datasets.

- Step 14 - Final Handling of Chemical Representations: Perform final handling of chemical identifiers due to missing data or redundancies resulting from previous molecular processing. Refer to the description for this step inside each notebook for details.

- Step 15- Additional Modifications: Adjust datasets to be more user-friendly. Export control datasets (without curation) for experiments.

  

**data_gathering.py**: supporting file for Step 2.

  

**filters.py**: supporting file for Filters used in the pipeline.

  

**utils.py**: supporting file for Step 9, 10, and 13.

  

**pubchem_sum**: this folder will be generated when running any of the Jupiter Notebook, associated with Step 11 (11.1).

  

**data**: this folder will be generated when running any of the Jupiter Notebook, containing all intermediate data and logger/statistics files for all curation steps.