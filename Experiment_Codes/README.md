# NeurIPS_Benchmark

## Guideline for Running Scheduler 


### Requirement 
1. Pytorch 
2. Torch-geometric 
3. tqdm
4. matplotlib 


### Data Preparation
1. Download the `dataset` folder and paste in the parent directory of NeurIPS_Benchmark
2. Check processed 2D/3D dataset e.g. `dataset/processed/2dmol/processed_2dmol_AID1798.pt` exists in the dataset
3. Check the split files e.g. `dataset/split/random/AID1798_2d_random_cv1.pt`

### Update and Run Scheduler 
1. Go over comments in `scheduler-ndslab4-gin-finetune.py` in the local_scripts
2. Update `back_to_experiment_folder()` function
3. In the NeurIPS_Benchmark, run by `local_scripts/scheduler-ndslab4-gin-finetune.py`

### Different Experiments
1. Fine Tuning: Add parameters options in the tunable lists in the `generate_configs` function
2. Scaffold vs. Random: Add different split scheme in the split_scheme
3. One_hot vs. Rich Features: change one_hot list 