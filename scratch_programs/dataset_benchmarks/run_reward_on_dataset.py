import pandas as pd
from rdkit import Chem
from importlib import import_module
import argparse
import yaml
import sys
import os
import traceback
from joblib import Parrallel, delayed
import traceback

parser = argparse.ArgumentParser(
    description="Get reward values for a benchmark dataset"
)
parser.add_argument("-s", "--smi", type=str, help="file with SMILES strings")
parser.add_argument("-c", "--config", type=str, help="config file (same as MCTS)")
args = parser.parse_args()

with open(args.config, "r") as f:
    conf = yaml.load(f, Loader=yaml.SafeLoader)
    
outdir = conf["output_dir"]
outfile = f"{outdir}/out.csv"
conf["debug"] = True

sys.path.insert(0, os.getcwd())
rs = conf["reward_setting"]
print(rs["reward_module"])
reward_calculator = getattr(import_module(rs["reward_module"]), rs["reward_class"])
smis = [i.rstrip() for i in open(args.smi, "r").readlines()]

def process_smi(count, smi, reward_calculator, conf, outfile):
    conf["gid"] = count
    try:
        mol = Chem.MolFromSmiles(smi)
        # this will also write out vina poses to outdir/3D_pose
        objective_values = [f(mol) for f in reward_calculator.get_objective_functions(conf)]
        print(objective_values)
        vinascore = objective_values[0]
        reward = reward_calculator.calc_reward_from_objective_values(values=objective_values, conf=conf)
    except Exception as e:
        traceback.print_exc()
        reward = -1.0
        objective_values = [-1.0]
        vinascore = -1.0

    print(reward, objective_values[0])
    with open(outfile, "a") as fw:
        fw.write(f"{count},{smi},{reward},{objective_values[0]}\n")

    return reward, vinascore

# Parallelize the loop using joblib
rewards, vinascores = zip(*Parallel(n_jobs=conf["n_jobs"])(delayed(process_smi)(count, smi, reward_calculator, conf, outfile)
                                                for count, smi in enumerate(smis)))
