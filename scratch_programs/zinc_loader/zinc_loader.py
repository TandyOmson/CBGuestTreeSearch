from datasets import load_dataset
from tqdm import tqdm
import pandas as pd
from rdkit import Chem
import argparse

def stream_and_collect(
    url,
    target_n,
    filter_fn,
    split="train",
):
    ds = load_dataset(
        "json",
        data_files=url,
        split=split,
        streaming=True
    )

    collected = []
    for example in tqdm(ds, desc=f"Collected {len(collected)} smis"):
        if filter_fn(example["smiles"]):
            collected.append(example)
            
            if len(collected) >= target_n:
                break

    return collected

url = "https://huggingface.co/datasets/zpn/zinc20/resolve/main/zinc_processed/smiles_all_*_clean.jsonl.gz"

def size_filter(smi, max_atoms=15):
    m = Chem.MolFromSmiles(smi)
    n_atoms = m.GetNumAtoms()
    if n_atoms > max_atoms:
        return False
    else:
        return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Load incrementally from the zinc database",
        usage="python {os.path.basename(__file__)} -n num_smiles -m max_atoms -o out_prefix"
    )
    parser.add_argument("-n", "--num_smiles", type=int, default=None)
    parser.add_argument("-m", "--max_atoms", type=int, default=None)
    parser.add_argument("-o", "--output", type=str, default=None)
    args = parser.parse_args()

    samples = stream_and_collect(
        url=url,
        target_n=args.num_smiles,
        filter_fn=lambda x: size_filter(x, args.max_atoms),
    )
    
    df = pd.DataFrame(samples)
    df.to_csv(f"{args.output}.csv")
    
    with open(f"{args.output}.smi", "w") as fw:
        for i in df["smiles"]:
            fw.write(f"{i}\n")
