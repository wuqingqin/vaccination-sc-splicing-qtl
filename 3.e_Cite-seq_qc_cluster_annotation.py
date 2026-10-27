#!/usr/bin/python3
import scanpy as sc
import numpy as np
import pandas as pd
import scanpy.external as sce
import muon as mu
import mudata as md
import gc
import random
import os
import re
from collections import Counter

# Run from the project root.
random.seed(123)
read_dir = "data/concat_all_singlets/mdata_All.h5mu"
mdata_All=mu.read(read_dir)

# Add timepoints and phenotype groups.
input_list=mdata_All.obs.Classification
output_list = []
for item in input_list:
    if re.search(r'PBMC.Hashtag1|PBMC.Hashtag2|PBMC.Hashtag3', item):
        output_list.append(re.sub(r'.*PBMC.Hashtag[123]', 'day0', item))
    elif re.search(r'PBMC.Hashtag6|PBMC.Hashtag4|PBMC.Hashtag5', item):
        output_list.append(re.sub(r'.*PBMC.Hashtag[456]', 'day7', item))
    else:
        output_list.append(item)
new_output_list = []
for item in output_list:
    if re.search(r'Hashtag7|Hashtag8|Hashtag9|Hashtag10|Hashtag12|Hashtag13', item):
        new_output_list.append(re.sub(r'.*.Hashtag([7-9]|10|12|13)', 'day7', item))
    elif re.search(r'Hashtag1|Hashtag2|Hashtag3|Hashtag4|Hashtag5|Hashtag6',item):
        new_output_list.append(re.sub(r'.*.Hashtag[123456]', 'day0', item))
    else:
        new_output_list.append(item)

Counter(new_output_list)
mdata_All.obs["timepoint"] = new_output_list
mdata_All.obs["barcode"] = mdata_All.obs.index.values

Phenotype_All=pd.read_csv("data/phenotype/Phenotype_All.csv")
Hashtag_e1=Phenotype_All["Hashtag_e1"]
Hashtag_e2=Phenotype_All["Hashtag_e2"]
Hashtag_e1 = [element.replace('_', '.') for element in Hashtag_e1]
Hashtag_e2 = [element.replace('_', '.') for element in Hashtag_e2]
four_Groups=Phenotype_All["four_Groups"]
Group_TA=Phenotype_All["Group_TA"]
Group_T=Phenotype_All["Group_T"]
people_id=mdata_All.obs["Classification"]
replacement_dict = dict(zip(Hashtag_e1, four_Groups))

people_id = [replacement_dict.get(val, val) for val in people_id]

replacement_dict = dict(zip(Hashtag_e2, four_Groups))

people_id = [replacement_dict.get(val, val) for val in people_id]
mdata_All.obs["four_Groups"]=people_id

people_id=mdata_All.obs["Classification"]

replacement_dict = dict(zip(Hashtag_e1, Group_TA))

people_id = [replacement_dict.get(val, val) for val in people_id]

replacement_dict = dict(zip(Hashtag_e2, Group_TA))

people_id = [replacement_dict.get(val, val) for val in people_id]
mdata_All.obs["Group_TA"]=people_id

people_id=mdata_All.obs["Classification"]

replacement_dict = dict(zip(Hashtag_e1, Group_T))

people_id = [replacement_dict.get(val, val) for val in people_id]

replacement_dict = dict(zip(Hashtag_e2, Group_T))

people_id = [replacement_dict.get(val, val) for val in people_id]
mdata_All.obs["Group_T"]=people_id

#rna
mdata_All["rna"].var['mt'] = mdata_All["rna"].var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(mdata_All["rna"], qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
mdata_All["rna"].obs["log10GenesPerUMI"] = np.log10(mdata_All["rna"].obs["n_genes_by_counts"].values)/np.log10(mdata_All["rna"].obs["total_counts"])
mdata_All.obs["nCount_RNA"] = mdata_All["rna"].obs["total_counts"]
mdata_All.obs["nFeature_RNA"] = mdata_All["rna"].obs["n_genes_by_counts"]
mdata_All.obs["percent.mt"] = mdata_All["rna"].obs["pct_counts_mt"]
mdata_All.obs["log10GenesPerUMI"] = mdata_All["rna"].obs["log10GenesPerUMI"]

#prot
mdata_All["prot"].X = mdata_All["prot"].layers["counts"].copy()
sc.pp.calculate_qc_metrics(mdata_All["prot"],percent_top=None,var_type="antibodies",inplace=True)
mdata_All["prot"].obs["log10AntibodiesPerUMI"] = np.log10(mdata_All["prot"].obs["n_antibodies_by_counts"].values)/np.log10(mdata_All["prot"].obs["total_counts"])
mdata_All.obs["nCount_ADT"] = mdata_All["prot"].obs["total_counts"]
mdata_All.obs["nFeature_ADT"] = mdata_All["prot"].obs["n_antibodies_by_counts"]
mdata_All.obs["log10AntibodiesPerUMI"] = mdata_All["prot"].obs["log10AntibodiesPerUMI"]

# Filter cells using the original QC thresholds.
#nFeature_RNA：200-2500  percent.mt<10  log10GenesPerUMI>0.8
mdata_filtered= mdata_All[(mdata_All.obs["log10GenesPerUMI"]>=0.8)&(mdata_All.obs["percent.mt"]<=10)&(mdata_All.obs.nFeature_RNA>=200)&(mdata_All.obs["nFeature_ADT"]>=100), :]#&(mdata_All.obs.nFeature_RNA<=2500)
print(mdata_filtered)
# mdata_de.obs['lane'].value_counts()
gc.collect()
# Keep features detected in at least 50 cells; remove mitochondrial/ribosomal genes.
sc.pp.calculate_qc_metrics(mdata_filtered["rna"], qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.calculate_qc_metrics(mdata_filtered["prot"],percent_top=None,var_type="antibodies",inplace=True)
mdata_filtered.var['mt'] = mdata_filtered.var_names.str.startswith('MT-')
mdata_filtered.var['ribo'] = mdata_filtered.var_names.str.startswith('RPS')+mdata_filtered.var_names.str.startswith('RPL')

mdata_filtered.var["n_cells_by_counts"] = pd.concat([mdata_filtered["prot"].var["n_cells_by_counts"][0:130],mdata_filtered["rna"].var["n_cells_by_counts"][0:36601]],axis=0)

mu.pp.filter_var(mdata_filtered,'n_cells_by_counts', lambda x: x >= 50)#12819
print(mdata_filtered)
mu.pp.filter_var(mdata_filtered,'mt', lambda x: x == False)#13
print(mdata_filtered)
mu.pp.filter_var(mdata_filtered,'ribo', lambda x: x == False)#103
print(mdata_filtered)

mdata_filtered["rna"].layers["counts"]=mdata_filtered["rna"].X

read_dir = "data/concat_all_singlets/mdata_filtered.h5mu"
mdata_filtered.write(read_dir)

gc.collect()

# Compute modality-specific PCA and UMAP.
read_dir = "data/concat_all_singlets/mdata_filtered.h5mu"
mdata_filtered=mu.read(read_dir)

mdata_filtered.obs["lane_sample"] = mdata_filtered.obs["lane"].astype("str")+"_"+mdata_filtered.obs["Classification"].astype("str")
# #rna
rna = mdata_filtered.mod["rna"]
rna.obs["lane_sample"] = mdata_filtered.obs["lane_sample"]
rna.obs["lane"] = mdata_filtered.obs["lane"]
rna.obs["date"] = mdata_filtered.obs["date"]
rna.X = rna.layers["counts"].copy()
sc.pp.normalize_total(rna,target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna,n_top_genes=2000)
rna.layers["data"] = rna.X
sc.pp.scale(rna, max_value=10)
rna.layers["scale.data"] = rna.X
sc.pp.pca(rna,mask_var="highly_variable")
rna.obsm["X_pca_raw"]=rna.obsm["X_pca"]
sc.pp.neighbors(rna, n_neighbors=30)
sc.tl.umap(rna)
rna.obsm["X_umap_raw"]=rna.obsm["X_umap"]

# #protein
protein = mdata_filtered.mod["prot"]
protein.obs["lane_sample"] = mdata_filtered.obs["lane_sample"]
protein.obs["lane"] = mdata_filtered.obs["lane"]
protein.obs["date"] = mdata_filtered.obs["date"]
protein.X = protein.layers["dsb"].copy()
sc.pp.pca(protein)
protein.obsm["X_pca_raw"]=protein.obsm["X_pca"]
sc.pp.neighbors(protein, n_neighbors=30)   # why can't we just work with the default neighbors?
sc.tl.umap(protein)
protein.obsm["X_umap_raw"]=protein.obsm["X_umap"]

mdata_filtered.mod["rna"] = rna
mdata_filtered.mod["prot"] = protein
mdata_filtered.write(read_dir)

# Correct lane-sample batches with Harmony.

os.environ["OPENBLAS_NUM_THREADS"] = "8"
print("harmony...")
#harmony
#rna
sc.external.pp.harmony_integrate(rna,'lane_sample')
sc.pp.neighbors(rna,use_rep="X_pca_harmony")
rna.obsm["X_pca"] = rna.obsm["X_pca_harmony"]
sc.tl.umap(rna)
rna.obsm["X_umap_harmony"]=rna.obsm["X_umap"]

print("rna harmony done!")
#protein
sc.external.pp.harmony_integrate(protein,'lane_sample')
sc.pp.neighbors(protein,use_rep="X_pca_harmony")
protein.obsm["X_pca"] = protein.obsm["X_pca_harmony"]
sc.tl.umap(protein)
protein.obsm["X_umap_harmony"]=protein.obsm["X_umap"]
print("prot harmony done!")

mdata_filtered.mod["rna"] = rna
mdata_filtered.mod["prot"] = protein
mdata_filtered.write(read_dir)
# Integrate modalities with WNN.

#WNN

# Calculate weighted nearest neighbors
mu.pp.neighbors(mdata_filtered, key_added='wnn')
mu.tl.umap(mdata_filtered,neighbors_key='wnn',random_state=1)
sc.tl.leiden(mdata_filtered, resolution=0.5, neighbors_key='wnn', key_added='leiden_wnn_0.5')
mdata_filtered.obsm["X_umap_wnn"] = mdata_filtered.obsm["X_umap"]
mdata_filtered.update()
mdata_filtered.write(read_dir)
print("wnn done!")

print("All Done!")
