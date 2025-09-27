""" DATASET AUGMENTATION
"""
import math
from rdkit import Chem
from sklearn.utils import class_weight
import sys

def augment_by_cluster(smis, cluster_labels, target_size):
    """ Augments existing dataset up to target_size by weight of frequency in cluster labels
        New dataset includes original
    """
    class_weights = class_weight.compute_class_weight(
        class_weight='balanced',
        classes=np.unique(cluster_labels),
        y=cluster_labels
    )
    
    # dictionary to pass to class_weight argument in keras model.fit
    class_weight_dict = dict(enumerate(class_weights, min(np.unique(cluster_labels))))

    new_smis = []
    for count, smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        
        # calculate number of smis to enumerate
        weight = class_weight_dict[cluster_labels[count]]
        enums_per_smi  = math.ceil((target_size/len(dataset))*weight - 1)
        
        enumerated_smis = [smi]
        for i in range(enums_per_smi):
            new_smi = Chem.MolToSmiles(mol, doRandom=True)
            if new_smi not in enumerated_smis:
                enumerated_smis.append(new_smi)
                
        new_smis.extend(enumerated_smis)

    return new_smis

def augment_evenly(smis, target_size):
    """ Augments existing dataset up to target_size evenly
        New dataset includes original
    """

    new_smis = []
    for count, smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        
        # calculate number of smis to enumerate
        enums_per_smi = math.ceil(target_size/len(dataset) - 1)
        
        enumerated_smis = [smi]
        for i in range(enums_per_smi):
            new_smi = Chem.MolToSmiles(mol, doRandom=True)
            if new_smi not in enumerated_smis:
                enumerated_smis.append(new_smi)
                
        new_smis.extend(enumerated_smis)

    return new_smis

dataset = [i.rstrip() for i in open(sys.argv[1], "r").readlines()]
new_dataset = augment_evenly(dataset, 50000)

print("new dataset lenght", len(new_dataset))

with open(sys.argv[2], "w") as fw:
    for i in new_dataset:
        fw.write(f"{i}\n")
