""" Random molecule generator using the Chemtsv2 RNN
    Playing around adding some restrictions on generated molecules, perhaps for dataset augmentation
"""

from logging import getLogger
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops 
from rdkit import RDLogger
import numpy as np
from numpy.random import default_rng
import pickle
import tensorflow as tf
tf.compat.v1.disable_eager_execution()

RDLogger.DisableLog('rdApp.*')

from chemtsv2.utils import load_tensorflow_model, get_model_structure_info
from chemtsv2.preprocessing import smi_tokenizer, selfies_tokenizer_from_smiles

def framework_filter(smi):
    """ Checking the criteria based on SMILES string:
        at least two rings intersecting
    """
    is_framework = False

    # check the number of rings
    rings = set([int(i) for i in smi if i.isnumeric()])
    ring_idxs = {i:[] for i in rings}
    if len(rings) < 3:
        return is_framework

    # get the indicies of ring start and terminations
    for i in rings:
        for count, c in enumerate(smi):
            if c.isnumeric():
                if i == int(c):
                    ring_idxs[i].append(count)
                    
        if i > 1:
            if ring_idxs[i][0] > ring_idxs[i-1][0] and ring_idxs[i][0] < ring_idxs[i-1][1] and ring_idxs[i][1] > ring_idxs[i-1][1] and ring_idxs[i-1][0] + 1 != ring_idxs[i][0]:
                print(smi, ring_idxs)
                is_framework = True
                return is_framework
            else:
                return is_framework

    return is_framework

def bridged_bicycles(smi):
    """ Checking the criteria based on graph
    """
    mol = Chem.MolFromSmiles(smi)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    num_paths = 0
    for i in range(len(atom_rings)):
        for j in range(i+1, len(atom_rings)):
            shared = set(atom_rings[i]).intersection(atom_rings[j])
            if len(shared) > 1:
                for path in Chem.rdmolops.FindAllPathsOfLengthN(mol, 3, useBonds=False):
                    if set(path) == set(shared):
                        num_paths += 1
                if num_paths >= 2:
                    return True
    return False

def ringsize_filter(smi):
    """ Returns True is a molecule has any 3 or 4 atom rings
    """
    has_small_rings = False
    
    mol = Chem.MolFromSmiles(smi)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    for i in atom_rings:
        if len(i) == 3 or len(i) == 4:
            has_small_rings = True
            return has_small_rings

    return has_small_rings


def generate_smiles_as_token_index(model, state, tokens, conf):
    end = "\n"
    position = []
    position.extend(state)
    get_int = [tokens.index(position[j]) for j in range(len(position))]
    x = np.reshape(get_int, (1, len(get_int)))
    model.reset_states()

    while not get_int[-1] == tokens.index(end):
        preds = model.predict_on_batch(x)
        state_pred = np.squeeze(preds)
        next_int = conf["random_generator"].choice(range(len(state_pred)), p=state_pred)
        get_int.append(next_int)
        x = np.reshape([next_int], (1, 1))
        if len(get_int) > conf["max_len"]:
            break
    return get_int

def build_smiles_from_token_index(generated_token_indexes, tokens, use_selfies=False):
    generate_tokens = [
        tokens[generated_token_indexes[j]] for j in range(len(generated_token_indexes) - 1)
    ]
    generate_tokens.remove("&")
    concat_tokens = "".join(generate_tokens)
    if use_selfies:
        # "[*]" is replaced with [Lr] because SELFIES (v2.1.0) currently does not support a wildcard representation.
        if "[Lr]" in concat_tokens:
            concat_tokens = sf.decoder(concat_tokens).replace("[Lr]", "[*]")
        else:
            concat_tokens = sf.decoder(concat_tokens)
    return concat_tokens

if __name__ == "__main__":

    model_dir = "/home/andyt/DProjects/DMCTS/VINA_ChemTSv2/model/HCs"

    logger = getLogger()

    conf = {
        "model_setting" :
        {
            "model_json" : f"{model_dir}/model.tf25.json",
            "model_weight" : f"{model_dir}/model.tf25.best.ckpt.h5",
        },
        "token" : f"{model_dir}/tokens.pkl",
        "random_generator" : default_rng(),
    }

    # If you want a stem molecule to produce molecules from
    input_smiles = None
    
    (
        conf["max_len"],
        conf["rnn_vocab_size"],
        conf["rnn_output_size"],
        conf["num_gru_units"],
    ) = get_model_structure_info(conf["model_setting"]["model_json"], logger)
    
    model = load_tensorflow_model(conf["model_setting"]["model_weight"], logger, conf)
    with open(conf["token"], "rb") as f:
        tokens = pickle.load(f)
    
    if input_smiles is not None:
        conf["input_smiles"] = args.input_smiles
        conf["tokenized_smiles"] = (
            selfies_tokenizer_from_smiles(conf["input_smiles"])
            if conf["use_selfies"]
            else smi_tokenizer(conf["input_smiles"])
        )

    # generate random compounds
    random_compounds = []
    num_random = 10000
    valid_count = 0

    # for testing purposes, specify input file
    #candidate_smi = [i.rstrip() for i in open("hydrophobe_hydrocarbons.smi", "r").readlines()]
    #candidate_smi = [i.rstrip() for i in open("hydros8992.smi", "r").readlines()]
    
    for i in range(num_random):
        print(f"running {i}", end="\r")
        new_compound = []
        generated_token_indexes = generate_smiles_as_token_index(
                        model, ["&"], tokens, conf
                    )
        
        new_compound.append(
            build_smiles_from_token_index(
                generated_token_indexes,
                tokens,
            )
        )

        # SMILES validity is done by checking rdkit sanitization, then filters are applied
        # see utils.evaluate_node
        if Chem.MolFromSmiles(new_compound[0]) == None:
            #print("invalid SMILES", new_compound)
            continue

        # Apply desired filters
        try:
            is_framework = bridged_bicycles(new_compound[0])
        except:
            is_framework = False

        try:
            has_small_rings = ringsize_filter(new_compound[0])
        except:
            has_small_rings = ringsize_filter(new_compound[0])

        if not is_framework or has_small_rings:
            continue
        else:
            print("framework!", new_compound[0])
        
        random_compounds.append(new_compound[0])
        valid_count += 1
        
    print(f"valid count, {valid_count} of {num_random}, {(valid_count/num_random)*100}%")
    print("\n")
    print(random_compounds)

    writer = Chem.SDWriter("valid.sdf")
    for i in random_compounds:
        mol = Chem.MolFromSmiles(i)
        mol = Chem.AddHs(mol)
        
        res = AllChem.EmbedMolecule(mol)
        if res == -1:
            res = AllChem.EmbedMolecule(mol, useBasicKnowledge=False)
            if res == -1:
                raise Exception
            else:
                AllChem.MMFFOptimizeMolecule(mol)
                
        writer.write(mol)
    writer.close()
