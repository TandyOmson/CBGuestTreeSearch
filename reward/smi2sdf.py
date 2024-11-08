#!/usr/bin/env python3

# generation of up to N low energy conformers
# from 2D input (smi) to 3D output (sdf)
# see Ebejer et. al.
# "Freely Available Conformer Generation Methods: How Good Are They?"
# JCIM, 2012, DOI: 10.1021/ci2004658 for technical details
#
# Copyright (C) 2018 Francois Berenger
# System Cohort Division,
# Medical Institute of Bioregulation,
# Kyushu University
# 3-1-1 Maidashi, Higashi-ku, Fukuoka 812-8582, Japan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

import sys
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = words[1]
            yield name, Chem.MolFromSmiles(smile)

# nb. conformers to generate prior to energy minimization
# as an empirical function of the molecule's flexibility
def how_many_conformers(mol):
    nb_rot_bonds = AllChem.CalcNumRotatableBonds(mol)
    if nb_rot_bonds <= 7:
        return 50
    elif nb_rot_bonds <= 12:
        return 200
    return 300  # This is more

# keep only conformers which are far enough from the reference conformer
def rmsd_filter(mol, ref_conf, conf_energies, threshold):
    # we use heavy atoms RMSD; not all atoms (Peter Gedeck's suggestion)
    mol_noH = Chem.RemoveHs(mol)
    ref_conf_id = ref_conf.GetId()
    res = []
    for e, curr_conf in conf_energies:
        curr_conf_id = curr_conf.GetId()
        rms = AllChem.GetConformerRMS(mol_noH, ref_conf_id, curr_conf_id)
        if rms > threshold:
            res.append((e, curr_conf))
    return res

def process_smi(mol, n_confs, rmsd_threshold):
    n = how_many_conformers(mol)
    #print(f"init pool size {n}", file = sys.stderr)
    mol_H = Chem.AddHs(mol)
    res = Chem.Mol(mol_H)
    res.RemoveAllConformers()
    #print("generating starting conformers ...", file = sys.stderr)
    conf_energies = []
    #print("FF minimization ...", file = sys.stderr)
    for cid in AllChem.EmbedMultipleConfs(mol_H, n):
        ff = AllChem.UFFGetMoleculeForceField(mol_H, confId = cid)
        # print("E before: %f" % ff.CalcEnergy())
        ff.Minimize()
        energy = ff.CalcEnergy()
        # print("E after: %f" % energy)
        conformer = mol_H.GetConformer(cid)
        # print("cid: %d e: %f" % (cid, energy))
        conf_energies.append((energy, conformer))
    # sort by increasing E
    conf_energies = sorted(conf_energies, key = lambda x: x[0])
    # output non neighbor conformers
    kept = 0
    #print("RMSD pruning ...", file = sys.stderr)
    while kept < n_confs and len(conf_energies) > 0:
        (e, conf) = conf_energies.pop(0)
        kept += 1
        cid = res.AddConformer(conf, assignId = True)
        # align conformers to the one of lowest energy
        if cid != 0:
            rdMolAlign.AlignMol(res, res, prbCid = cid, refCid = 0)
        # remove neighbors
        conf_energies = rmsd_filter(mol_H, conf, conf_energies, rmsd_threshold)
    #print("kept %d confs" % (kept), file = sys.stderr)

    return res