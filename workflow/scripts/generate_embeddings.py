""""
Generating Fused Network Embeddings
"""
import sys
import pandas as pd
import pyarrow.parquet as pq

from scipy import stats
from itertools import combinations
from damply import dirs
import pandas as pd
import tqdm 
import numpy as np
import argparse
import random
import torch 
from mvfusion.mv_integrator import MultiViewIntegrator
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import SimpleImputer, IterativeImputer
import time 
import yaml
# import snf2


with open(dirs.CONFIG / "embeddings.yaml", 'r') as file:
    # Load the contents safely
    data = yaml.safe_load(file)


hidden_dim = data['PARAMETERS']['hidden_dimension']
mu = data['PARAMETERS']['mu']
num_nbrs = data['PARAMETERS']['num_nbrs']
num_epochs = data['PARAMETERS']['num_epochs']
SEED = data['PARAMETERS']['seed']
metric = data['PARAMETERS']['metric']

## Initialize the seeds
SEED = 1234
np.random.seed(SEED)
torch.manual_seed(SEED)
random.seed(SEED)



colData = pd.read_csv( 
            dirs.RAWDATA / "HDD"/ "colData.tsv", 
            sep="\t",
            usecols= ['HDD.Compound.ID','Mechanism.of.Action',"Pubchem.CID",
                      "LINCS.CMap.Name",'GEOM.Source.SMILES',"JUMP.CP.ID"],
                      engine='python'   )



###
#
# Mapping sample names back to HDD Compound ID
####


LINCS_map = colData[["HDD.Compound.ID","LINCS.CMap.Name"]]
LINCS_map  = LINCS_map.dropna(subset=["LINCS.CMap.Name"])
GEOM_map = colData[["HDD.Compound.ID",'GEOM.Source.SMILES']]
VIABILITY_map = colData[["HDD.Compound.ID","Pubchem.CID"]]
VIABILITY_map.columns = ["HDD.Compound.ID","PubChem.CID"]
JUMP_map = colData[["HDD.Compound.ID","JUMP.CP.ID"]]
JUMP_map = JUMP_map.dropna(subset=['JUMP.CP.ID'])





####
#
# Prep JUMP
#
####
print("\nLoading in the JUMP data...\n")
JUMP_meta = pd.read_csv(
    dirs.RAWDATA / "JUMP" / "colData.tsv",
    sep="\t",
    usecols=["Sample.ID", "JUMP.CP.ID"],
)
# Apply the existing MoA eligibility filter before computing independent means.
moa_counts = colData["Mechanism.of.Action"].value_counts()
eligible_ids = colData.loc[
    colData["Mechanism.of.Action"].isin(moa_counts[moa_counts > 1].index),
    "HDD.Compound.ID",
]
JUMP_map = JUMP_map[JUMP_map["HDD.Compound.ID"].isin(eligible_ids)]
JUMP_mapping = JUMP_map.merge(JUMP_meta, on="JUMP.CP.ID")

# Read batches so the full JUMP assay does not need to fit in memory.
JUMP_sums = None
JUMP_counts = None
for batch in pq.ParquetFile(dirs.RAWDATA / "JUMP" / "cpcnn.parquet").iter_batches(
    batch_size=4096
):
    frame = JUMP_mapping.merge(batch.to_pandas(), on="Sample.ID")
    if frame.empty:
        continue
    frame = frame.drop(columns=["Sample.ID", "JUMP.CP.ID"])
    grouped = frame.set_index("HDD.Compound.ID").astype("float64").groupby(level=0)
    batch_sums = grouped.sum()
    batch_counts = grouped.count()
    JUMP_sums = (
        batch_sums if JUMP_sums is None else JUMP_sums.add(batch_sums, fill_value=0)
    )
    JUMP_counts = (
        batch_counts
        if JUMP_counts is None
        else JUMP_counts.add(batch_counts, fill_value=0)
    )
if JUMP_sums is None:
    raise ValueError("No JUMP samples mapped to HDD compounds")
JUMP_data = JUMP_sums.div(JUMP_counts.where(JUMP_counts > 0)).sort_index()
JUMP_mols = set(JUMP_data.index)


####
#
# Prep LINCS
#
####


print("\nLoading in the LINCS data...\n")
LINCS_data = pd.read_parquet(dirs.RAWDATA / "LINCS"/ "signatures.parquet")


LINCS_data = LINCS_map.merge(LINCS_data,on= "LINCS.CMap.Name")
LINCS_data = LINCS_data.drop(labels=["LINCS.CMap.Name"],axis=1)


LINCS_data = LINCS_data.set_index('HDD.Compound.ID')
LINCS_mols = set(LINCS_data.index)



####
#
# Prep GEOM
#
####

print("\nLoading in the GEOM data..")
GEOM_data = pd.read_parquet(dirs.RAWDATA / "GEOM"/ "WHIM_hp.parquet")

GEOM_data = GEOM_map.merge(GEOM_data,on='GEOM.Source.SMILES')
GEOM_data = GEOM_data.drop(labels=['GEOM.Source.SMILES'],axis=1)

GEOM_mols = set(GEOM_data['HDD.Compound.ID'])
GEOM_data = GEOM_data.set_index('HDD.Compound.ID')



###
#
#  Prep Viability
#
###
print("\nLoading in the NCI60 data..")


VIABILITY = pd.read_csv(dirs.PROCDATA / "NCI60_2026.csv",index_col=0)

VIABILITY = VIABILITY.dropna(subset=['PubChem.CID'])
# print(VIABILITY.head())
# sys.exit()
VIABILITY = VIABILITY_map.merge(VIABILITY,on="PubChem.CID")
VIABILITY = VIABILITY.drop(labels=['PubChem.CID','treatmentid'],axis=1)
VIABILITY = VIABILITY.set_index('HDD.Compound.ID')


X = VIABILITY.values[:,:20]
print("\nStarting imputation of missing values...\n")
imp = IterativeImputer(missing_values=np.nan,random_state=SEED)
# imp = SimpleImputer(missing_values=np.nan)
s= time.time()
X = imp.fit_transform(X)

e = time.time()
print(f"\t... imputation took {e-s} seconds.\n")
VIABILITY = pd.DataFrame(X, index=VIABILITY.index,columns = VIABILITY.columns[:20])
viability_mols = set(VIABILITY.index)


###
#
# Set Up Data Labels
#
####


colData = colData.dropna(subset=['Mechanism.of.Action'])
count = pd.DataFrame(colData['Mechanism.of.Action'].value_counts()).reset_index()
count = count[count['count']>1]
colData = colData[colData['Mechanism.of.Action'].isin(count['Mechanism.of.Action'])]



keep_mols = set(colData['HDD.Compound.ID'])

keep_mols = keep_mols.intersection(GEOM_mols)
keep_mols = keep_mols.intersection(LINCS_mols)
keep_mols = keep_mols.intersection(JUMP_mols)
keep_mols = keep_mols.intersection(viability_mols)
keep_mols = sorted(list(keep_mols))
print(len(keep_mols))
LINCS = LINCS_data.loc[keep_mols,]

print(LINCS.head())
print(LINCS.shape)

GEOM = GEOM_data.loc[keep_mols,]
print(GEOM.head())
print(GEOM.shape)

JUMP = JUMP_data.loc[keep_mols,]
print(JUMP.head())
print(JUMP.shape)



RESP = VIABILITY.loc[keep_mols,]
RESP = RESP.groupby(level=0).mean()

views = [GEOM, LINCS, JUMP, RESP]


integrator = MultiViewIntegrator(
            views = views,
            view_names=['GEOM','LINCS','JUMP','NCI60'],
                metrics = [metric]*len(views),
                neighborhood_size= num_nbrs,
                mu= mu,
                alignment_epochs=num_epochs,
                emb_dim =hidden_dim ,
                seed = SEED)



embeds_final, S_final, model = integrator.neural_integration()

embeds_final.to_csv(dirs.RESULTS / "molecule_embeddings.csv")
S_final = pd.DataFrame(S_final, index= embeds_final.index, columns= embeds_final.index)

print("All Done!")