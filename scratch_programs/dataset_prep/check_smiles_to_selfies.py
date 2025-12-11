""" Checks whether a SMILES file can be translated to SELFIES, optionally produce a trimmed file
"""

import argparse
import selfies as sf

def selfies_tokenizer_from_smiles(smi):
    if "[*]" in smi:
        # Because SELFIES (v2.1.0) currently does not support a wildcard (*) representation.
        smi = smi.replace("[*]", "[Lr]")
    slfs = sf.encoder(smi)
    tokens = list(sf.split_selfies(slfs))
    assert slfs == "".join(tokens)
    tokens.insert(0, "&")
    return tokens

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="checks SMILES are translated to SELFIES",
        usage="python {os.path.basename(__file__)} -s SMILES_FILE -o OUT_SMILES"
    )
    parser.add_argument(
        "-s",
        "--smilesfile",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default=False
    )
    args = parser.parse_args()

    smis = [i.rstrip() for i in open(args.smilesfile, 'r').readlines()]

    len_smis = len(smis)
    failed_smis = 0
    failed_idxs = []
    passed_smiles = []
    for count, s in enumerate(smis):
        print(f"Processing SMILES {count} of {len_smis}", end="\r")
        try:
            selfies_tokenizer_from_smiles(s)
            passed_smiles.append(s)
        except:
            failed_smis += 1
            failed_idxs.append(count)
            
    print(f"\n\n Removed {failed_smis} of {len_smis}\n")
    print(failed_idxs)
    if args.outfile:
        with open(args.outfile, "w") as fw:
            for s in passed_smiles:
                fw.write(f"{s}\n")
