#!/usr/bin/python3
import scanpy as sc
import numpy as np
import pandas as pd
import scanpy.external as sce
import muon as mu
from muon import prot as pt
import mudata as md
import gc
import random
import os
import concurrent.futures
from pathlib import Path

def find_outs_directories(root_path, max_depth=None, parallel=False):
    """Find outs directories, optionally limiting depth or searching in parallel."""
    outs_directories = []
    
    if not os.path.exists(root_path):
        print(f"Error: path {root_path} does not exist")
        return outs_directories
    
    base_depth = len(Path(root_path).parts)
    
    def _search_dir(current_path, current_depth):
        if max_depth is not None and current_depth - base_depth > max_depth:
            return []
        
        results = []
        try:
            with os.scandir(current_path) as entries:
                for entry in entries:
                    if entry.is_dir():
                        if entry.name == "outs":
                            results.append(entry.path)
                        results.extend(_search_dir(entry.path, current_depth + 1))
        except PermissionError:
            print(f"Warning: cannot access directory {current_path}")
        return results
    
    def _search_dir_parallel(dir_path):
        return _search_dir(dir_path, base_depth + 1)
    
    if parallel:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            top_level_dirs = []
            try:
                with os.scandir(root_path) as entries:
                    for entry in entries:
                        if entry.is_dir():
                            top_level_dirs.append(entry.path)
            except PermissionError:
                print(f"Warning: cannot access root directory {root_path}")
            
            results = list(executor.map(_search_dir_parallel, top_level_dirs))
            outs_directories = [item for sublist in results for item in sublist]
    else:
        outs_directories = _search_dir(root_path, base_depth)
    
    return outs_directories
# Run from the project root; adjust the relative Cell Ranger input below.
random.seed(123)
##############################################################################################################################
# Cell Ranger input root.
root_directory = "data/cellranger/last_cellranger"

outs_dirs = find_outs_directories(root_directory, parallel=True)

if outs_dirs:
    print(f"Found {len(outs_dirs)} outs directories:")
    for idx, directory in enumerate(outs_dirs, 1):
        print(f"{idx}. {directory}")
else:
    print("No outs directories found")    
exclude_dirs = [
    "data/cellranger/last_cellranger/12_lanes_0606_L002/0606-A1_GS_L002/outs",
    "data/cellranger/last_cellranger/12_lanes_0606_L002/0606-A4_GS_L002/outs",
    "data/cellranger/last_cellranger/12_lanes_0606_L002/0606-B1_GS_L002/outs",
    "data/cellranger/last_cellranger/12_lanes_0606_L002/0606-B3_GS_L002/outs",
    "data/cellranger/last_cellranger/12_lanes_0606_L002/0606-B2_GS_L002/outs",
    "data/cellranger/last_cellranger/12_lanes_0606_L002/0606-A3_GS_L002/outs",
    "data/cellranger/last_cellranger/12_lanes_0606_L002/0606-A2_GS_L002/outs",
]

filtered_dirs = [d for d in outs_dirs if d not in exclude_dirs]
filtered_dirs = sorted(filtered_dirs, key=lambda x: (0, x) if 'PBMC' in x else (1, x))
lane =   ["PBMC_01","PBMC_02",
          "0407_01","0407_02",
          "0410_A3","0410_A4","0410_B1","0410_B2","0410_B3","0410_B4",
          "0604_A1","0604_A2","0604_A3","0604_A4","0604_B1","0604_B2","0604_B3","0604_B4",
          "0606_B4","0606_A1","0606_A2","0606_A3","0606_A4","0606_B1","0606_B2","0606_B3",
          "0609_A1","0609_A2","0609_A3","0609_A4","0609_B1","0609_B2","0609_B3","0609_B4",
          "0611_A1","0611_A2","0611_A3","0611_A4","0611_B1","0611_B2","0611_B3","0611_B4",
          "0628_A1","0628_A2","0628_A3","0628_A4","0628_B1","0628_B2","0628_B3","0628_B4",
          "0630_A1","0630_A2","0630_A3","0630_A4","0630_B1","0630_B2","0630_B3","0630_B4",
          "0704_A1","0704_A2","0704_A3","0704_A4","0704_B1","0704_B2","0704_B3","0704_B4",
          "0407_01_24","0407_02_24","0407_03_24","0407_04_24",
          "0606_A1_24","0606_A2_24","0606_A3_24","0606_A4_24","0606_A5_24","0606_B1_24","0606_B2_24","0606_B3_24","0606_B4_24",
         ]

date =   ["PBMC","PBMC",
          "0407","0407",
          "0410","0410","0410","0410","0410","0410",
          "0604","0604","0604","0604","0604","0604","0604","0604",
          "0606","0606","0606","0606","0606","0606","0606","0606",
          "0609","0609","0609","0609","0609","0609","0609","0609",
          "0611","0611","0611","0611","0611","0611","0611","0611",
          "0628","0628","0628","0628","0628","0628","0628","0628",
          "0630","0630","0630","0630","0630","0630","0630","0630",
          "0704","0704","0704","0704","0704","0704","0704","0704",
          "0407_24","0407_24","0407_24","0407_24",
          "0606_24","0606_24","0606_24","0606_24","0606_24","0606_24","0606_24","0606_24","0606_24",
         ]

Classification=  ["PBMC","PBMC",
                  "0407","0407",
                  "0410A","0410A","0410B","0410B","0410B","0410B",
                  "0604A","0604A","0604A","0604A","0604B","0604B","0604B","0604B",
                  "0606B","0606A","0606A","0606A","0606A","0606B","0606B","0606B",
                  "0609A","0609A","0609A","0609A","0609B","0609B","0609B","0609B",
                  "0611A","0611A","0611A","0611A","0611B","0611B","0611B","0611B",
                  "0628A","0628A","0628A","0628A","0628B","0628B","0628B","0628B",
                  "0630A","0630A","0630A","0630A","0630B","0630B","0630B","0630B",
                  "0704A","0704A","0704A","0704A","0704B","0704B","0704B","0704B",
                  "0407","0407","0407","0407",
                  "0606A","0606A","0606A","0606A","0606A","0606B","0606B","0606B","0606B",
                 ]

for i in range(len(filtered_dirs)):
    data_dir = filtered_dirs[i]
    mdata = mu.read_10x_mtx(os.path.join(data_dir,"filtered_feature_bc_matrix"))
    mdata_raw = mu.read_10x_mtx(os.path.join(data_dir,"raw_feature_bc_matrix"))
    mdata['prot'].layers['counts'] = mdata['prot'].X
    mdata['rna'].layers['counts'] = mdata['rna'].X
    
    # Demultiplex cells using hashtag counts.
    cell_hashing_columns=["Hashtag1","Hashtag2","Hashtag3","Hashtag4","Hashtag5","Hashtag6","Hashtag7","Hashtag8","Hashtag9","Hashtag10","Hashtag11","Hashtag12"]
    Hashtag = mdata['prot'].layers['counts'][:,0:12]
    Hashtag = Hashtag.toarray()
    mdata.obs[cell_hashing_columns] = pd.DataFrame(Hashtag[:,0:12]).values
    sce.pp.hashsolo(mdata,cell_hashing_columns=cell_hashing_columns)
    
    mdata.obs["Classification"] = Classification[i]+"."+mdata.obs["Classification"]
    mdata.obs["barcode"] = mdata.obs.index.values
    mdata.obs=mdata.obs[["barcode","most_likely_hypothesis","Classification"]]
    mdata.obs["lane"] = lane[i]
    mdata.obs["date"] = date[i]
    # Normalize ADT counts with DSB.
    save_dir = os.path.join("data/hashsolo_dsb_single", lane[i])
    os.makedirs(save_dir, exist_ok=True)
    mdata_raw['rna'].obs["log10umi"] = np.array(np.log10(mdata_raw['rna'].X.sum(axis=1) + 1)).reshape(-1)

    isotypes=[gene for gene in mdata.var_names if 'Isotype' in gene]
    pt.pp.dsb(mdata, mdata_raw, empty_counts_range=(1,2), isotype_controls=isotypes, random_state=1)
    mdata['prot'].layers['dsb'] = mdata['prot'].X

    # Prefix barcodes with lane IDs and make feature names unique.
    mdata.obs.index = lane[i]+"_"+mdata.obs.index
    mdata.mod["rna"].obs.index=lane[i]+"_"+mdata.mod["rna"].obs.index
    mdata.mod["prot"].obs.index=lane[i]+"_"+mdata.mod["prot"].obs.index
    mdata.var_names_make_unique()
    mdata.mod["prot"].var_names_make_unique()
    mdata.mod["rna"].var_names_make_unique()
    
    # Keep singlets and remove hashtag/isotype features.
    mdata_singlet = mdata[mdata.obs["most_likely_hypothesis"]==1.0,]
    columns_to_delete = cell_hashing_columns+isotypes
    var_names = mdata_singlet.var_names.tolist()
    indexes_to_delete = [var_names.index(col) for col in columns_to_delete if col in var_names]
    
    new_var_names = [i for i in range(len(var_names)) if i not in indexes_to_delete]
    mdata_singlet = mdata_singlet[:,new_var_names]
    print(mdata_singlet)

    # Save processed lane objects.
    
    md.set_options(pull_on_update=False)
    mdata.write(os.path.join(save_dir, lane[i]+"_mdata_raw.h5mu"))
    mdata_singlet.write(os.path.join(save_dir, lane[i]+"_mdata_singlet.h5mu"))
    print(lane[i]+" save Done!")
    gc.collect()
print("All Done!")
