#!/usr/bin/python3
import muon as mu
import mudata as md
import gc
import os
import random
def retry_concat(mdata_merge):
    """Concatenate MuData objects, retrying once on failure."""
    try:
        mdata_All = md.concat(mdata_merge, merge="unique")
        return mdata_All
    except Exception as e:
        print(f"First concatenation failed: {e}; retrying once")
        try:
            mdata_All = md.concat(mdata_merge, merge="unique")
            return mdata_All
        except Exception as e2:
            print(f"Second concatenation failed: {e2}")
            return None

# Run from the project root.
random.seed(123)
names =   ["PBMC_01","PBMC_02",
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
dir = "data/hashsolo_dsb_single"
mdata_merge=[]

for i in range(len(names)):
  mdata=mu.read(os.path.join(dir, names[i], names[i]+"_mdata_singlet.h5mu"))
  print(mdata)
  mdata.obs["barcode"] = mdata.obs.index.values
  mdata.mod["rna"].obs=mdata.obs
  mdata.mod["prot"].obs=mdata.obs
  mdata_merge=mdata_merge+[mdata]
  print(names[i]+" merge done!")
  gc.collect()

mdata_All=retry_concat(mdata_merge)
print(mdata_All)

os.makedirs("data/concat_all_singlets", exist_ok=True)
mdata_All.write("data/concat_all_singlets/"+"mdata_All.h5mu")
print("save done!")
print("All Done!")
gc.collect()
