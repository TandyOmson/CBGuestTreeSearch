""" Testing for reward module
"""
import os
import sys
import yaml
import argparse
import joblib
import copy
from importlib import import_module

from rdkit import Chem
import numpy as np
from rdkit.Chem.MolStandardize import rdMolStandardize

# Load conf
def get_parser():
    parser = argparse.ArgumentParser(
        description="", usage=f"python {os.path.basename(__file__)} -c CONFIG_FILE"
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=True,
        help="path to a config file",
    )
    parser.add_argument(
        "--input_smiles",
        type=str,
        required=True,
        help="SMILES string (Need to put the atom you want to extend at the end of the string)",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="debug mode",
    )
    return parser.parse_args()

args = get_parser()
with open(args.config, "r") as f:
    conf = yaml.load(f, Loader=yaml.SafeLoader)
    conf["debug"] = args.debug

# Set working directory to chemtsv2 base directory (this is where we would be running from)
os.chdir("./../")
# modify sys path for imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

conf["mol_embed_pipeline"] = os.getcwd() + "/" + conf["mol_embed_pipeline"]
conf["kde"] = os.getcwd() + "/" + conf["kde"]
conf["vina_boltz_utils"] = os.getcwd() + "/" + conf["vina_boltz_utils"]

# Setup environement for call of reward module
# This is in the /chemtsv2/utils.py evaulate_mode method
rs = conf["reward_setting"]
reward_calculator = getattr(import_module(rs["reward_module"]), rs["reward_class"])
ps = conf["policy_setting"]
policy_evaluator = getattr(import_module(ps["policy_module"]), ps["policy_class"])

# Load the list of smiles to run
smi = args.input_smiles
print("Running SMILES:", smi)

_conf = copy.deepcopy(conf)
_conf["gid"] = 1

mol = Chem.MolFromSmiles(smi)

valid_mol_list = [mol]
valid_conf_list = [_conf]

# Run the reward function
def _get_objective_values(mol, conf):
    return [f(mol) for f in reward_calculator.get_objective_functions(conf)]

if conf["leaf_parallel"]:
    # the parallel is setup in mcts.py
    parallel = joblib.Parallel(
                    n_jobs=conf["leaf_parallel_num"],
                    prefer="processes",
                )
    
    values_list = parallel(
        joblib.delayed(_get_objective_values)(m, c)
        for m, c in zip(valid_mol_list, valid_conf_list)
    )
elif conf["batch_reward_calculation"]:
    values_list = [
        f(valid_mol_list, valid_conf_list)
        for f in reward_calculator.get_batch_objective_functions()
    ]
    values_list = np.array(values_list).T.tolist()
    
else:
    values_list = [_get_objective_values(m, c) for m, c in zip(valid_mol_list, valid_conf_list)]

print("output:", values_list)
