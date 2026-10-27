#UMAP visualization for all cell types
import muon as mu

mdata = mu.read_h5mu(
    "/public/users/jiangxinyu/Project/COVID19IV_CITE-seq_202408/data/concat_all_singlets/qwt/mdata_qwt_m0_l.h5mu"
)

print(mdata)

import pandas as pd
meta = pd.read_csv(
    "/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/mdata_filtered_metadata.tsv.gz",
    sep="\t"
)

keep_cells = mdata.obs_names.intersection(meta["rna:barcode" ])

mdata_sub = mdata[keep_cells, :].copy()

import scanpy as sc
import matplotlib as mpl

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
   "font.family": "serif",
    "font.serif": ["serif"],  
    "axes.unicode_minus": False,

    "font.size": 20,

    "axes.titlesize": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20
})
ax = mu.pl.embedding(
    mdata_sub,
    basis="umap_wnn",
    color="celltype_8",
    legend_fontsize=16,
    frameon=False,
    show=False
)

for coll in ax.collections:
    coll.set_rasterized(True)

ax.figure.savefig(
    "/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/umap/all_umap_celltype.pdf",
    bbox_inches="tight",
    dpi=300
)

#UMAP visualization for five major cell type subsets
import muon as mu
import pandas as pd
import os
import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib as mpl

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
   "font.family": "serif",
    "font.serif": ["serif"],  
    "axes.unicode_minus": False,

    "font.size": 20,

    "axes.titlesize": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20
})
# -----------------------------
# input files
# -----------------------------
files = {
    "CD4": "/public/users/jiangxinyu/Project/COVID19IV_CITE-seq_202408/data/concat_all_singlets/qwt/mdata_qwt_m0_CD4.h5mu",
    "CD8": "/public/users/jiangxinyu/Project/COVID19IV_CITE-seq_202408/data/concat_all_singlets/qwt/mdata_qwt_m0_CD8.h5mu",
    "NK": "/public/users/jiangxinyu/Project/COVID19IV_CITE-seq_202408/data/concat_all_singlets/qwt/mdata_qwt_m0_NK.h5mu",
    "Myeloid": "/public/users/jiangxinyu/Project/COVID19IV_CITE-seq_202408/data/concat_all_singlets/qwt/mdata_qwt_m0_Myeloid.h5mu",
    "B": "/public/users/jiangxinyu/Project/COVID19IV_CITE-seq_202408/data/concat_all_singlets/qwt/mdata_qwt_m0_B.h5mu"
}

# -----------------------------
# metadata
# -----------------------------
meta = pd.read_csv(
    "/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/mdata_filtered_metadata.tsv.gz",
    sep="\t"
)

outdir = "/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/umap"
os.makedirs(outdir, exist_ok=True)

# -----------------------------
# loop
# -----------------------------
for name, path in files.items():

    print(f"Processing: {name}")

    mdata = mu.read_h5mu(path)

    # keep cells
    keep_cells = mdata.obs_names.intersection(meta["rna:barcode"])
    mdata_sub = mdata[keep_cells, :].copy()

    # fix category
    if "celltype" in mdata_sub.obs.columns:
        mdata_sub.obs["celltype"] = (
            mdata_sub.obs["celltype"]
            .astype("category")
            .cat.remove_unused_categories()
        )

    # plot
    ax = mu.pl.embedding(
        mdata_sub,
        basis="umap_wnn",
        color="celltype",
        frameon=False,
        show=False
    )

    # rasterize points
    for coll in ax.collections:
        coll.set_rasterized(True)

    # save pdf
    outpdf = os.path.join(outdir, f"umap_celltype_{name}.pdf")

    ax.figure.savefig(
        outpdf,
        bbox_inches="tight",
        dpi=300
    )

    plt.close(ax.figure)

    print(f"Saved: {outpdf}")
