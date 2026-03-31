import pandas as pd
import numpy as np

num_mols = max([int(i.split(",")[0]) if i.rstrip() != "FAIL" else 0 for i in open("docked/affins_all.csv").readlines()])+1

df_affins = pd.read_csv("docked/affins_all.csv", names=[f"pose_{i}" for i in range(0,num_mols+1)], engine="python").drop(columns=["pose_0"])
df_affins.index = range(1, len(df_affins)+1)

df_affins["avg_score"] = df_affins.apply(np.mean, axis=1)

df_affins.to_csv("all_affins.csv")

df_rmsds = pd.read_csv("docked/rmsds_all.csv", names=[f"pose_{i}" for i in range(0,num_mols+1)], engine="python").drop(columns=["pose_0"])
df_rmsds.index = range(1, len(df_rmsds)+1)

df_rmsds.to_csv("all_rmsds.csv")

df_smina = pd.read_csv("smina_breakdown/smina_all.csv", names=["1","2","3","4","5"], engine="python")
df_smina.index = range(1, len(df_smina)+1)

df_smina.to_csv("all_smina.csv")
