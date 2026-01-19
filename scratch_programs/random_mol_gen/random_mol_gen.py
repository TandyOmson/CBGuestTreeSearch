""" Random molecule generator using the Chemtsv2 RNN
    Playing around adding some restrictions on generated molecules, perhaps for dataset augmentation
"""

import yaml
import argparse

from logging import getLogger
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops 
from rdkit import RDLogger
import selfies as sf
import numpy as np
from numpy.random import default_rng
import pickle
import tensorflow as tf
from tensorflow.keras.models import Sequential, model_from_json  # pyright: ignore[reportMissingImports]
from tensorflow.keras.layers import Dense, Embedding, GRU  # pyright: ignore[reportMissingImports]
tf.compat.v1.disable_eager_execution()

RDLogger.DisableLog('rdApp.*')

from chemtsv2.preprocessing import smi_tokenizer, selfies_tokenizer_from_smiles

def load_tensorflow_model(model_weight, logger, conf):
    model = Sequential()
    model.add(
        Embedding(
            input_dim=conf["rnn_vocab_size"],
            output_dim=conf["rnn_vocab_size"],
            mask_zero=False,
            batch_size=1,
        )
    )
    model.add(
        GRU(
            conf["units_GRU_1"],
            batch_input_shape=(1, None, conf["rnn_vocab_size"]),
            activation="tanh",
            return_sequences=True,
            stateful=True,
        )
    )
    model.add(
        GRU(
            conf["units_GRU_2"],
            activation="tanh",
            return_sequences=False,
            stateful=True,
        )
    )
    model.add(
        Dense(
            conf["rnn_output_size"],
            activation="softmax",
        )
    )
    model.load_weights(model_weight)
    logger.info(f"Model weights loaded from {model_weight}")

    return model

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

    parser = argparse.ArgumentParser(
        description="Generate random mols using a trained RNN",
        usage="python {os.path.basename(__file__)} -m MODEL_DIRECTORY -o OUTFILE"
    )
    parser.add_argument(
        "-m",
        "--modeldir",
        type=str,
        required=True,
        help="path to directory with model",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    outfile = args.outfile
    model_dir = args.modeldir

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

    # Get RNN config from the same input file as train_rnn.py
    with open(conf["model_setting"]["model_json"], "r") as fr:
        loaded_model_json = fr.read()
    loaded_model = model_from_json(loaded_model_json)

    for layer in loaded_model.get_config()["layers"]:
        config = layer.get("config")
        if layer.get("class_name") == "InputLayer":
            conf["max_len"] = config["batch_input_shape"][1]
        if layer.get("class_name") == "Embedding":
            conf["rnn_vocab_size"] = config["input_dim"]
        if layer.get("class_name") == "TimeDistributed":
            conf["rnn_output_size"] = config["layer"]["config"]["units"]

#    # To specify a config file, otherwise assumes rnn_model_setting.yaml is in the model dir
#    parser = argparse.ArgumentParser(description="", usage="python3 random_mol_gen.py -c RNN_CONFIG_FILE")
#    parser.add_argument("-c", "--config", type=str, required=True, help="path to RNN config file")
#    args = parser.parse_args()

#    with open(args.config, "r") as f:
#        rnn_conf = yaml.load(f, Loader=yaml.SafeLoader)
#        conf.update(rnn_conf)

    
    rnn_config_file = model_dir + "/rnn_model_setting.yaml"
    with open(rnn_config_file, "r") as f:
        rnn_conf = yaml.load(f, Loader=yaml.SafeLoader)
        conf.update(rnn_conf)

    # Reconstruct model
    model = load_tensorflow_model(conf["model_setting"]["model_weight"], logger, conf)

    # Get tokens
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
    num_random = 10000
    valid_count = 0
    count = 0

    gen_smis = []
    with open(outfile, "w") as fw:
        while valid_count < num_random:
            print(f"running: {count}, valid count: {valid_count}", end="\r")
            new_compound = []
            generated_token_indexes = generate_smiles_as_token_index(
                            model, ["&"], tokens, conf
                        )
            
            new_compound.append(
                build_smiles_from_token_index(
                    generated_token_indexes,
                    tokens,
                    use_selfies=conf["use_selfies"]
                )
            )
            
            # SMILES validity is done by checking rdkit sanitization
            # see utils.evaluate_node
            if Chem.MolFromSmiles(new_compound[0]) == None:
                #print("invalid SMILES", new_compound)
                count +=1
                continue
            
            fw.write(f"{new_compound[0]}\n")
            gen_smis.append(new_compound[0])
            count += 1
            valid_count += 1

    #check for duplicates
    duplicate_count = 0
    canon_gen_smis = []
    for smi in gen_smis:
        mol = Chem.MolFromSmiles(smi)
        canon_smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True, canonical=True)
        if canon_smi in canon_gen_smis:
            duplicate_count += 1
            
    print(f"valid count, {valid_count} of {count}, {(valid_count/count)*100}%")
    print(f"{duplicate_count} duplicate molecules")
