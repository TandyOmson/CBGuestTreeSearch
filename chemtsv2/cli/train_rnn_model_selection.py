""" Train RNN with hyperparameter tuning
"""

import argparse
import os
import pickle

import tensorflow as tf
from tensorflow.keras.models import Sequential  # pyright: ignore[reportMissingImports]
from tensorflow.keras.layers import Dense, Embedding, GRU, TimeDistributed  # pyright: ignore[reportMissingImports]
from tensorflow.keras.optimizers import Adam  # pyright: ignore[reportMissingImports]
from tensorflow.keras.utils import to_categorical  # pyright: ignore[reportMissingImports]
from tensorflow.keras.preprocessing import sequence  # pyright: ignore[reportMissingImports]
from tensorflow.keras.callbacks import CSVLogger, EarlyStopping, ModelCheckpoint  # pyright: ignore[reportMissingImports]
import numpy as np
import yaml
import keras_tuner

from chemtsv2.preprocessing import read_smiles_dataset, tokenize_smiles

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"  # Hide all GPUs
#tf.config.threading.set_intra_op_parallelism_threads(44)
#tf.config.threading.set_inter_op_parallelism_threads(2)

def get_parser():
    parser = argparse.ArgumentParser(description="", usage="chemtsv2-train-rnn -c CONFIG_FILE")
    parser.add_argument("-c", "--config", type=str, required=True, help="path to a config file")
    return parser.parse_args()

def prepare_data(smiles, all_smiles):
    """TODO: need to be refactored"""
    all_smiles_index = []
    for i in range(len(all_smiles)):
        smiles_index = []
        for j in range(len(all_smiles[i])):
            smiles_index.append(smiles.index(all_smiles[i][j]))
        all_smiles_index.append(smiles_index)
    X_train = all_smiles_index
    y_train = []
    for i in range(len(X_train)):
        x1 = X_train[i]
        x2 = x1[1 : len(x1)]
        x2.append(0)
        y_train.append(x2)
    return X_train, y_train

def save_model(model, output_dir, use_selfies=False):
    output_json = os.path.join(
        output_dir, "model_sf.tf25.json" if use_selfies else "model.tf25.json"
    )
    output_weight = os.path.join(output_dir, "model_sf.tf25.h5" if use_selfies else "model.tf25.h5")
    model_json = model.to_json()

    with open(output_json, "w") as json_file:
        json_file.write(model_json)
    print(f"[INFO] Save a model structure to {output_json}")

    model.save_weights(output_weight)
    print(f"[INFO] Save model weights to {output_weight}")

class RnnHyperModel(keras_tuner.HyperModel):
    def __init__(self, token_list, X_len):
        self.token_list = token_list
        self.X_len = X_len

    def build(self, hp):
        """ RNN for predicting the next token in a sequence
            Batch size N, sequence size T
            Sample n, sequence position t
    
            Defines layers, callbacks, compiles model (currently no callbacks)
        """
        model = Sequential()
        # Input embedding layer for token probabilities at given sequence position
        model.add(
            Embedding(
                input_dim=len(self.token_list),
                output_dim=len(self.token_list),
                input_length=self.X_len,
                mask_zero=False,
            )
        )
        # Gated reccurent units (hidden layers)
        model.add(
            GRU(
                units=hp.Choice("units_GRU_1", values=[64, 128, 256, 512]),
                input_shape=(self.X_len, len(self.token_list)),
                activation="tanh",
                dropout=hp.Float("dropout_rate_GRU_1", min_value=0.0, max_value=0.5, step=0.1),
                recurrent_dropout=hp.Float("rec_dropout_GRU_1", min_value=0.0, max_value=0.5, step=0.1),
                return_sequences=True,
            )
        )
        model.add(
            GRU(
                units=hp.Choice("units_GRU_2", values=[64, 128, 256, 512]),
                activation="tanh",
                dropout=hp.Float("dropout_rate_GRU_2", min_value=0.0, max_value=0.5, step=0.1),
                recurrent_dropout=hp.Float("rec_dropout_GRU_2", min_value=0.0, max_value=0.5, step=0.1),
                return_sequences=True,
            )
        )
        # Output softmax distribution for probability of each token at given sequence position
        model.add(
            TimeDistributed(
                Dense(
                    len(self.token_list),
                    activation="softmax",
                )
            )
        )
        model.build(input_shape=(self.X_len, len(self.token_list)))
        model.summary()
            
        # Compile model
        model.compile(
            loss="categorical_crossentropy",
            optimizer=Adam(
                #learning_rate=hp.Choice("learning_rate", values=[1e-4, 3e-4, 1e-3, 3e-3])
                learning_rate=1e-3,
            ),
            metrics=["accuracy"],
        )
    
        return model

    def fit(self, hp, model, *args, **kwargs):
        return model.fit(
            *args,
            # Could add 16 for very small datasets
            batch_size=hp.Choice("batch_size", values=[256, 128, 64, 32]),
            **kwargs,
        )

def main():
    args = get_parser()

    # Setup configuration
    with open(args.config, "r") as f:
        conf = yaml.load(f, Loader=yaml.SafeLoader)
    print("========== Configuration ==========")
    for k, v in conf.items():
        print(f"{k}: {v}")
    print("===================================")

    os.makedirs(conf["output_model_dir"], exist_ok=True)

    # Prepare training dataset
    original_smiles_list = read_smiles_dataset(conf["dataset"])
    token_list, tokenized_smiles_list = tokenize_smiles(original_smiles_list, use_selfies=conf["use_selfies"])
    assert len(original_smiles_list) == len(tokenized_smiles_list)
    print(f"[INFO] Size of training dataset: {len(original_smiles_list)}")
    if conf["use_selfies"]:
        base, ext = os.path.splitext(conf["output_token"])
        conf["output_token"] = f"{base}_sf{ext}"
    with open(conf["output_token"], "wb") as f:
        pickle.dump(token_list, f)
    print(f"[INFO] Generated tokens: {token_list}")
    if conf["use_selfies"]:
        print(
            f"[INFO] Save generated tokens to {conf['output_token']}. "
            "Note that the file name was modified because `use_selfies` was specified."
        )
    else:
        print(f"[INFO] Save generated tokens to {conf['output_token']}")

    # Validation set indicies
    if conf["val_set_idxfile"]:
        val_set_idxs = [int(i.rstrip()) for i in open(conf["val_set_idxfile"], "r").readlines()]
    else:
        val_set_idxs = []

    tokenized_smiles_list_train = [tokenized_smiles_list[i] for i in range(len(tokenized_smiles_list)) if i not in val_set_idxs]
    tokenized_smiles_list_val = [tokenized_smiles_list[i] for i in val_set_idxs]

    X_train, y_train = prepare_data(token_list, tokenized_smiles_list_train)
    X = sequence.pad_sequences(
        X_train,
        maxlen=conf["maxlen"],
        dtype="int32",
        padding="post",
        truncating="pre",
        value=0.0,
    )
    y = sequence.pad_sequences(
        y_train,
        maxlen=conf["maxlen"],
        dtype="int32",
        padding="post",
        truncating="pre",
        value=0.0,
    )
    y_train_one_hot = np.array([
        to_categorical(sent_label, num_classes=len(token_list)) for sent_label in y
    ])
    print(f"[DEBUG] Shape of y_train_one_hot: {y_train_one_hot.shape}")

    # Prepare validation dataset
    X_val_raw, y_val_raw = prepare_data(token_list, tokenized_smiles_list_val)
    X_val = sequence.pad_sequences(
        X_val_raw,
        maxlen=conf["maxlen"],
        dtype="int32",
        padding="post",
        truncating="pre",
        value=0.0,
    )
    y_val = sequence.pad_sequences(
        y_val_raw,
        maxlen=conf["maxlen"],
        dtype="int32",
        padding="post",
        truncating="pre",
        value=0.0,
    )
    y_val_one_hot = np.array([
        to_categorical(sent_label, num_classes=len(token_list)) for sent_label in y_val
    ])
    print(f"[DEBUG] Shape of y_val_one_hot: {y_val_one_hot.shape}")

    # Hyperparameter search
    tuner = keras_tuner.RandomSearch(
        hypermodel = RnnHyperModel(token_list, X.shape[1]),
        objective = 'val_loss',
        max_trials=conf["num_trials"],
        executions_per_trial=1,
        overwrite=True,
        directory="hyperparms"
        )

    # define callbacs
    log_path = os.path.join(conf["output_model_dir"], "training_log.csv")
    logger = CSVLogger(log_path)
    early_stopping = EarlyStopping(monitor="val_accuracy", patience=5)
    model_ckpt = ModelCheckpoint(
        filepath=os.path.join(
            conf["output_model_dir"],
            "model_sf.tf25.best.ckpt.h5" if conf["use_selfies"] else "model.tf25.best.ckpt.h5",
        ),
        monitor="val_accuracy",
        verbose=1,
        save_best_only=True,
        save_weights_only=False,
        mode="max",
        save_freq="epoch",
    )
    callbacks = [early_stopping]

    tuner.search(X,
                 y_train_one_hot,
                 epochs = conf["epoch"],
                 validation_data=(X_val, y_val_one_hot),
                 callbacks=callbacks,
                 shuffle=True,
                 )
        
    models = tuner.get_best_models(num_models=3)
    for count, i in enumerate(models):
        print("model {}:".format(count))
        print(i.summary())
        print("\n")

if __name__ == "__main__":
    main()
