#!/usr/bin/python3
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import gc
import pandas as pd
import scanpy as sc
import muon as mu

# Run from the project root. Completed cell annotations are required.
# Sources: qwt_umap.ipynb, qwt_s0_umap_compute.py, qwt_m0_umap_compute.py.
input_file = "data/concat_all_singlets/mdata_filtered_light.h5mu"
annotation_file = "data/concat_all_singlets/mdata_filtered_metadata.pkl"
phenotype_file = "data/phenotype/Phenotype_All.csv"
output_file = "data/concat_all_singlets/wqq/mdata_subset_100.h5mu"

mdata = mu.read(input_file)
metadata = pd.read_pickle(annotation_file)
if not metadata.index.is_unique:
    raise ValueError("Annotation metadata contain duplicate cell barcodes.")
if not mdata.obs_names.isin(metadata.index).all():
    raise ValueError("Annotation metadata do not cover all input cell barcodes.")
mdata.obs["celltype"] = metadata["celltype"].reindex(mdata.obs_names)

# Both timepoints are excluded for five participants; reasons are not documented
# in the source code. Preserve the original exclusion list.
remove_samples = [
    "0611A.Hashtag5", "0611A.Hashtag12",
    "0628B.Hashtag4", "0628B.Hashtag10",
    "0611B.Hashtag1", "0611B.Hashtag7",
    "0611B.Hashtag6", "0611B.Hashtag13",
    "0609B.Hashtag1", "0609B.Hashtag7",
]
# Preserve the original spelling for matching stored labels.
remove_celltypes = ["RBC", "unkown"]
if "Classification" not in mdata.obs.columns:
    raise KeyError("Classification is missing from mdata.obs.")

mask_class = ~mdata.obs["Classification"].isin(remove_samples)
mask_celltype = ~mdata.obs["celltype"].isin(remove_celltypes)
n_cells_before = mdata.n_obs
mdata = mdata[mask_class & mask_celltype].copy()
print(f"Cells before filtering: {n_cells_before}; after exclusions: {mdata.n_obs}")
gc.collect()

# Keep the fine cell types specified in the original m0 workflow.
target_celltypes = [
    "CD4.Tn", "CD4.CD69+Tn", "CD4.Tcm", "Th1-like.Tem", "Th2-like.Tem",
    "Th1", "Th2", "Th17", "Treg",
    "CD8.Tn", "CD8.CD69+Tn", "CD8.Tcm", "CD8.Tem", "CD8.Teff",
    "cMono", "ncMono", "cDC", "pDCs", "CCL5+Myeloid",
    "CD56dim", "KIR2DL1+CD56dim", "KIR3DL1+CD56dim", "CD56bright",
    "γδT", "MAIT",
    "B.Naive", "B.SM", "B.MYC+SM", "B.IgD-IgM+UM", "Plasma",
]
mdata = mdata[mdata.obs["celltype"].isin(target_celltypes)].copy()
print(f"Cells retained in target cell types: {mdata.n_obs}")

# Apply the original main-cell-type mapping in the same order.
celltype_mapping = {
    "CD4": "CD4", "Th1": "CD4", "Th2": "CD4", "Th17": "CD4", "Treg": "CD4",
    "CD8": "CD8", "Mono": "Myeloid", "DC": "Myeloid", "CCL5+Myeloid": "Myeloid",
    "CD56": "NK", "KIR": "NK", "γδT": "γδT", "MAIT": "MAIT",
    "B.": "B", "Plasma": "B",
}
maincelltypes = []
for celltype in mdata.obs["celltype"]:
    maincelltype = celltype
    for key, value in celltype_mapping.items():
        if key in celltype:
            maincelltype = value
            break
    maincelltypes.append(maincelltype)
mdata.obs["maincelltype"] = maincelltypes

# Verify the retained 100-person cohort against the phenotype table.
phenotype = pd.read_csv(phenotype_file, dtype=str)
sample_d0 = phenotype["Hashtag_e1"].str.replace("_", ".", regex=False)
sample_d7 = phenotype["Hashtag_e2"].str.replace("_", ".", regex=False)
excluded = sample_d0.isin(remove_samples) | sample_d7.isin(remove_samples)
sample_to_participant = dict(zip(sample_d0, phenotype["unique_num"]))
sample_to_participant.update(dict(zip(sample_d7, phenotype["unique_num"])))
retained_samples = pd.Index(mdata.obs["Classification"].unique())
retained_participants = retained_samples.map(sample_to_participant)
if retained_participants.isna().any():
    raise ValueError("Some retained sample labels do not match the phenotype table.")
expected_participants = set(phenotype.loc[~excluded, "unique_num"])
if len(expected_participants) != 100 or set(retained_participants) != expected_participants:
    raise ValueError("Retained participants do not match the expected 100-person cohort.")
print(f"Participants excluded: {excluded.sum()}; retained: {len(expected_participants)}")

# Recompute RNA preprocessing on the retained cells.
rna = mdata.mod["rna"]
rna.obs["lane"] = mdata.obs["lane"]
rna.X = rna.layers["counts"].copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes=2000)
rna.layers["data"] = rna.X.copy()
sc.pp.scale(rna, max_value=10)
sc.pp.pca(rna, mask_var="highly_variable")
rna.obsm["X_pca_raw"] = rna.obsm["X_pca"].copy()
rna.X = rna.layers["data"].copy()
gc.collect()

# Recompute ADT PCA from the existing DSB-normalized layer.
protein = mdata.mod["protein"]
protein.obs["lane"] = mdata.obs["lane"]
protein.X = protein.layers["dsb"].copy()
sc.pp.pca(protein)
protein.obsm["X_pca_raw"] = protein.obsm["X_pca"].copy()

# Correct lane batches and recompute modality-specific UMAP embeddings.
print("Running RNA Harmony...")
sc.external.pp.harmony_integrate(rna, "lane")
sc.pp.neighbors(rna, use_rep="X_pca_harmony")
rna.obsm["X_pca"] = rna.obsm["X_pca_harmony"]
sc.tl.umap(rna)
rna.obsm["X_umap_harmony"] = rna.obsm["X_umap"].copy()
print("RNA Harmony and UMAP complete.")

print("Running ADT Harmony...")
sc.external.pp.harmony_integrate(protein, "lane")
sc.pp.neighbors(protein, use_rep="X_pca_harmony")
protein.obsm["X_pca"] = protein.obsm["X_pca_harmony"]
sc.tl.umap(protein)
protein.obsm["X_umap_harmony"] = protein.obsm["X_umap"].copy()
print("ADT Harmony and UMAP complete.")

mdata.mod["rna"] = rna
mdata.mod["protein"] = protein

# Rebuild WNN using the original RNA/ADT PC counts and UMAP parameters.
sc.pp.neighbors(mdata.mod["rna"], use_rep="X_pca_harmony", n_pcs=20)
sc.pp.neighbors(mdata.mod["protein"], use_rep="X_pca_harmony", n_pcs=15)
mu.pp.neighbors(mdata, key_added="wnn", n_neighbors=20)
mu.tl.umap(mdata, neighbors_key="wnn", random_state=1, min_dist=1.2, spread=1.5)
mdata.obsm["X_umap_wnn"] = mdata.obsm["X_umap"].copy()
mdata.update()
print("Subset WNN and UMAP complete.")

os.makedirs(os.path.dirname(output_file), exist_ok=True)
mdata.write(output_file)
print(f"Saved: {output_file}")
gc.collect()
