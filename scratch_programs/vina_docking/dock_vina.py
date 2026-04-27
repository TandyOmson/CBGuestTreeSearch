""" Module for docking using AutoDock Vina
Inputs: Molecule object, Molecule ID, Host topology file, Vina input file
Outputs: Docking output file, and file with best pose
Returns: Complex molecule object
"""

import subprocess as sp
import os
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdMolAlign

from scipy.spatial import Delaunay
from sklearn.decomposition import PCA
import math
import numpy as np
import vina
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate

import tempfile

from contextlib import closing
import argparse

import pandas as pd
import re

def add_nitrogen_charges(m):
    m.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(m)
    if not ps:
        Chem.SanitizeMol(m)
        return m
    for p in ps:
        if p.GetType()=='AtomValenceException':
            at = m.GetAtomWithIdx(p.GetAtomIdx())
            if at.GetAtomicNum()==7 and at.GetFormalCharge()==0 and at.GetExplicitValence()==4:
                at.SetFormalCharge(1)
            if at.GetAtomicNum()==7 and at.GetFormalCharge()==0:
                bondcount = 0
                for b in at.GetBonds():
                    bondcount += b.GetBondTypeAsDouble()
                if int(bondcount) > 3:
                    at.SetFormalCharge(1)
                
    Chem.SanitizeMol(m)
    return m

def write_out_confs(confs, writer, n_poses):
    n, i = 0, 0
    while n < n_poses:
        try:
            writer.write(confs[i])
            n += 1
            i += 1
        except:
            break
    return n

def parse_vina_output(filename, exo_ignore):
    affinities = []
    rmsd_lbs = []
    rmsd_ubs = []

    with open(filename, 'r') as file:
        lines = file.readlines()

    parsing = False
    count = 0
    for line in lines:
        # Start parsing after the table header
        if re.match(r'^\s*-+\+-+', line):
            parsing = True
            continue
        if parsing:
            if line.strip().split()[0] == "Writing":
                break  # Stop at first empty line after data
            parts = line.strip().split()
            if len(parts) == 4:
                if count not in exo_ignore:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lb = float(parts[2])
                    rmsd_ub = float(parts[3])

                    affinities.append({"affinity": affinity})
                    rmsd_lbs.append({"rmsd_lb": rmsd_lb})
                    rmsd_ubs.append({"rmsd_ub": rmsd_ub})
                count += 1

    return pd.DataFrame(affinities), pd.DataFrame(rmsd_lbs), pd.DataFrame(rmsd_ubs)

def parse_smina_output(posesfile, exo_ignore):
    affinities = []
    rmsds = []
    atom_data_dfs = []

    ms = Chem.SDMolSupplier(posesfile, removeHs=False, sanitize=False)
    mols = [add_nitrogen_charges(m) for m in ms]

    for count, mol in enumerate(mols):
        if count in exo_ignore:
            mols.pop(count)
        else:
            affinity = mol.GetDoubleProp("minimizedAffinity")
            affinities.append(affinity)

            rmsd = rdMolAlign.CalcRMS(mol, mols[0])
            rmsds.append(rmsd)
            
            atom_data = mol.GetProp("atomic_interaction_terms")
            atom_data_df = pd.DataFrame([x.split() for x in atom_data.split("\n")][1:-1], columns=[x for x in atom_data.split("\n")[0].split(" ")])
            atom_data_dfs.append(atom_data_df)

    return affinities, rmsds, atom_data_dfs

def is_exo(mol, confId=-1):
    """ Checks for exo complex
    """
    centroiddiffthreshold=4
    cavityatomsthreshold=6
    
    hostmol = Chem.GetMolFrags(mol, asMols=True)[0]

    Chem.RemoveHs(mol)
    Chem.RemoveHs(hostmol)

    # Separate host and guest, get their coordinates
    guest_coords = np.array([mol.GetConformer().GetAtomPosition(atm.GetIdx()) for count, atm in enumerate(mol.GetAtoms()) if count >= hostmol.GetNumAtoms()])
    host_coords = np.array([mol.GetConformer().GetAtomPosition(atm.GetIdx()) for count, atm in enumerate(mol.GetAtoms()) if not count >= hostmol.GetNumAtoms()])

    # Get host and guest centroid
    guest_centroid = np.array(guest_coords).mean(axis=0)
    host_centroid = np.array(host_coords).mean(axis=0)

    # Calculate distance between guest and host centroid
    centroid_diff = np.linalg.norm(guest_centroid - host_centroid)

    # Delauny defines the convex hull of the host atoms
    # Delauny is a triangulation such that none of the host atoms are inside the circumsphere of any tetrahedron in the traingulation
    hull = Delaunay(host_coords)
    
    # Calculate number of atoms in cavity
    cavity_atoms = 0
    for atm in guest_coords:
        # Points outside the triangulation return -1
        if hull.find_simplex(atm) >= 0:
            cavity_atoms += 1

    # Check if exclusion complex
    isExo = False
    if mol.GetNumAtoms() > 141:
        if centroid_diff > centroiddiffthreshold or cavity_atoms < cavityatomsthreshold:
            isExo = True
    else:
        if centroid_diff > centroiddiffthreshold:
            isExo = True

    return isExo

# ALIGNMENT DOCKING FUNCTIONS
def align_mol(mol):
    """ Align principal axis of a molecule along the z axis 
    """
    guest_atoms = [atm.GetSymbol() for atm in mol.GetAtoms()]
    guest_coords = np.array([mol.GetConformer().GetAtomPosition(atm.GetIdx()) for atm in mol.GetAtoms()])

    #add clouds of points around the atoms of the guest molecule, 3 sections in polar (theta) and 6 in azimuthal (phi)
    #This prepares the coordinates for the PCA
    #   x = r*sin(theta)*cos(phi)
    #	y = r*sin(theta)*sin(phi)
    #	z = r*cos(theta)
    atomic_radii = {'H' :1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'Cl':1.75, 'Br':1.85, 'I':1.98, 'P':1.8, 'S':1.8, 'As':1.85, 'B':2.13, 'Si':2.1, 'Se':1.9, 'Te':2.06}
    
    cloud_points = [] 
    for index,atm in enumerate(guest_atoms):
        cloud_radius = atomic_radii[atm]
        # Add a cloud point above and below the atom
        cloud_points.append(guest_coords[index] + np.array([0,0,cloud_radius]))
        cloud_points.append(guest_coords[index] + np.array([0,0,-cloud_radius]))
        # Add a cloud of points around the atom, 3 in theta and 6 in phi
        for theta in np.linspace(0,np.pi,3):
            for phi in np.linspace(0,2*np.pi,6):
                cloud_points.append(guest_coords[index] + np.array([cloud_radius*np.sin(theta)*np.cos(phi),cloud_radius*np.sin(theta)*np.sin(phi),cloud_radius*np.cos(theta)]))

    # Add cloud points to the guest coordinates
    cloud_points = np.array(cloud_points)
    guest_coords_with_clouds = np.concatenate((guest_coords,cloud_points),axis=0)

    # Initiliase PCA
    pca = PCA(n_components=3)
    # Fit PCA to guest coordinates and transform the coordinates
    pca.fit_transform(guest_coords_with_clouds)
    transform_coord = pca.transform(guest_coords)

    # Direct the principal axis of the guest molecule towards the z-axis (the axis pointing through the cavity of the host)
    theta = np.arctan2(transform_coord[0,0],transform_coord[0,2])
    rotation_matrix = np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]])
    transform_coord = np.matmul(rotation_matrix,transform_coord.T).T

    # Centre the transformed coordinates on the host centroid (0,0,0)
    transform_coord_centered = transform_coord.copy()
    transform_coord_centered[:,0] = transform_coord[:,0] - np.mean(transform_coord[:,0])
    transform_coord_centered[:,1] = transform_coord[:,1] - np.mean(transform_coord[:,1])
    transform_coord_centered[:,2] = transform_coord[:,2] - np.mean(transform_coord[:,2])

    return transform_coord_centered

def PCA_align_pose(hostmol, guest):
    """ Generates binding poses via alignment docking
    """
    aligned_coords = align_mol(guest)
    
    conf = Chem.Conformer(guest.GetConformer(0))
    for index, atm in enumerate(guest.GetAtoms()):
        conf.SetAtomPosition(atm.GetIdx(), aligned_coords[index])
    guest.RemoveAllConformers()
    guest.AddConformer(conf, assignId=True)

    binding_pose = Chem.CombineMols(hostmol, guest)
    binding_pose = add_nitrogen_charges(binding_pose)

    return binding_pose

def MMFF94_vina_opt(mol, rundir):
    """ Runs MMFF94, then optimisation on a binding pose
        Modifies the pose in the conformers and adds a vina score
    """
    hostmol, guestmol = Chem.GetMolFrags(mol, asMols=True)

    mmff_res = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, ignoreInterfragInteractions=False)
    
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(guestmol)
    
    for setup in mol_setups:
        pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup, bad_charge_ok=True)

    if not is_ok:
        print(err_msg)

    vinaobj = vina.Vina(verbosity=0)
    vinaobj.set_receptor(f"{rundir}/host.pdbqt")
    vinaobj.compute_vina_maps(center=[0.0,0.0,0.0], box_size=[20.0, 20.0, 20.0])
    
    vinaobj.set_ligand_from_string(pdbqt_string)
    vina_decomp = vinaobj.optimize()
    vina_en = -vina_decomp[0]

    vinaobj.write_pose(f"{rundir}/temp_align_pose.pdbqt", overwrite=True, remarks="")
    pdbqt_mol = PDBQTMolecule.from_file(f"{rundir}/temp_align_pose.pdbqt", skip_typing=True)
    complexmol = Chem.CombineMols(hostmol, RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0])

    # Returns an rdkit mol
    return [complexmol], [vina_en]
        
def dock(repsdffile, ligsdffile, inp, rundir):
    """ Dock a molecule into a host
    hostfile is a .pdbqt file, use command "obabel host.pdb -O host.pdbqt -xrh" to convert
    Sets the vina interaction energy as a MolDoubleProp
    """
    # Possible improvment: run this all in one shell process, with && between lines
    sp.run(["obabel",repsdffile,"-O","host.pdbqt","-xr", "-xh"], cwd=rundir, stderr=open(f"{rundir}/obabel.log","w"))
    sp.run(["obabel",ligsdffile,"-O","guest.pdbqt","-xh"], cwd=rundir, stderr=open(f"{rundir}/obabel.log","a"))
    
    # Run Vina
    sp.run(["vina",f"--config={inp}","--num_modes=100","--energy_range=25","--receptor=host.pdbqt",f"--ligand=guest.pdbqt",f"--out=poses.pdbqt"], cwd=rundir, stdout=open(f"{rundir}/vina.log","w"))

    # Revert output .pdbqt to .sdf
    sp.run(["obabel","poses.pdbqt","-O","poses.sdf"], cwd=rundir, stderr=open(f"{rundir}/obabel.log","a"))

    # Disable RDKit logging
    RDLogger.DisableLog('rdApp.*')

    # Read in all binding poses
    guestmols = Chem.SDMolSupplier(f"{rundir}/poses.sdf",removeHs=False,sanitize=False)
    hostmol = Chem.MolFromMolFile(repsdffile,removeHs=False,sanitize=False)

    complexmols = []
    exos = []
    for count, mol in enumerate(guestmols):
        comp = Chem.CombineMols(hostmol, mol)
        comp = add_nitrogen_charges(comp) 
        if is_exo(comp) == True:
            exos.append(count)
        else:
            complexmols.append(comp)

    if len(complexmols) != 0:
        affinities, rmsd_lbs, rmsd_ubs = parse_vina_output(f"{rundir}/vina.log", exo_ignore=exos)        
        affinities = affinities.rename(columns={"affinity":str(len(complexmols))}).T
        rmsd_lbs = rmsd_lbs.rename(columns={"rmsd_lb":str(len(complexmols))}).T
        rmsd_ubs = rmsd_ubs.rename(columns={"rmsd_ub":str(len(complexmols))}).T

    # use aligment docking as as last resort
    else:
        print("docking by alignment as a last resort...")
        one_guestmol = Chem.MolFromMolFile(ligsdffile, removeHs=False)
        one_guestmol = add_nitrogen_charges(one_guestmol)
        align_pose = PCA_align_pose(hostmol, one_guestmol)
        complexmols, vina_ens = MMFF94_vina_opt(align_pose, rundir)

        affinities = pd.DataFrame({"affinity": affinity for affinity in vina_ens}, index=[0]).rename(columns={"affinity":str(len(complexmols))}).T
    
    return complexmols, affinities

def dock_smina(repsdffile, ligsdffile, inp, rundir):
    """Dock a molecule into a host using smina"""

    # Convert to pdbqt
    sp.run(
        ["obabel", repsdffile, "-O", "host.pdbqt", "-xr", "-xh"],
        cwd=rundir,
        stderr=open(f"{rundir}/obabel.log", "w")
    )
    sp.run(
        ["obabel", ligsdffile, "-O", "guest.pdbqt", "-xh"],
        cwd=rundir,
        stderr=open(f"{rundir}/obabel.log", "a")
    )

    # Run smina (note: direct SDF output + better scoring)
    sp.run([
        "smina",
        f"--config={inp}",
        "--num_modes=100",
        "--energy_range=25",
        "--receptor=host.pdbqt",
        "--ligand=guest.pdbqt",
        "--out=poses.sdf",
        "--log=smina.log",
        "--scoring=vinardo",            
        "--atom_term_data"
    ],
    cwd=rundir,
    stdout=open(f"{rundir}/smina.stdout", "w"),
    stderr=open(f"{rundir}/smina.stderr", "w")
    )

    # Disable RDKit logging
    RDLogger.DisableLog('rdApp.*')

    # Read poses
    guestmols = Chem.SDMolSupplier(
        f"{rundir}/poses.sdf",
        removeHs=False,
        sanitize=False
    )
    hostmol = Chem.MolFromMolFile(
        repsdffile,
        removeHs=False,
        sanitize=False
    )

    complexmols = []
    exos = []

    for count, mol in enumerate(guestmols):
        if mol is None:
            continue

        comp = Chem.CombineMols(hostmol, mol)
        comp = add_nitrogen_charges(comp)

        if is_exo(comp):
            exos.append(count)
        else:
            complexmols.append(comp)

    if len(complexmols) != 0:
        affinities, rmsds, atom_data_dfs = parse_smina_output(
            f"{rundir}/poses.sdf",
            exo_ignore=exos
        )

    else:
        print("docking by alignment as a last resort...")

        one_guestmol = Chem.MolFromMolFile(
            ligsdffile,
            removeHs=False
        )
        one_guestmol = add_nitrogen_charges(one_guestmol)

        align_pose = PCA_align_pose(hostmol, one_guestmol)
        complexmols, affinities = MMFF94_vina_opt(align_pose, rundir)
        rmsds, atom_data_dfs = [], []

    return complexmols, affinities, rmsds, atom_data_dfs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "dock a .sdf file into a host, return up to n options")
    parser.add_argument("-n",metavar = "n_poses", type = int, default = 1,
                        dest = "n_poses",
                        help = "docking poses to generate (default = 1)")
    parser.add_argument("-r", metavar = "input_receptor_sdf", dest = "input_sdf_rep")
    parser.add_argument("-l", metavar = "input_ligand_sdf", dest = "input_sdf_lig")
    parser.add_argument("-o", metavar = "output_sdf", dest = "output_sdf")
    parser.add_argument("-i", metavar = "vina_input", dest = "vina_inp")
    parser.add_argument("-a", metavar = "affins_output", type = str, default = "affinities.csv", dest = "affins_output")
    parser.add_argument("-d", metavar = "rmsds_output", type = str, default = "rmsds.csv", dest = "rmsds_output")
    parser.add_argument("-s", metavar = "smina_output", type = str, default = "smina_breakdown.csv", dest = "smina_output")
    parser.add_argument("-b", metavar = "atom_terms_output", type = str, default = "atom_terms_output.csv", dest = "atom_terms_output")
    parser.set_defaults(rename=True)

    args = parser.parse_args()
    n_poses = args.n_poses
    #input_sdf_lig = os.getcwd() + "/" + args.input_sdf_lig
    #input_sdf_rep = os.getcwd() + "/" + args.input_sdf_rep
    input_sdf_lig = args.input_sdf_lig
    input_sdf_rep = args.input_sdf_rep
    output_sdf = args.output_sdf
    inp = args.vina_inp
    affins_output = args.affins_output
    rmsds_output = args.rmsds_output
    smina_output = args.smina_output
    atom_terms_output = args.atom_terms_output
    
    suppl = Chem.SDMolSupplier(input_sdf_lig, sanitize=False, removeHs=False)
    print("docking {} conformers".format(len(suppl)))
    print("max {} complexes".format(n_poses))
    
    with closing(Chem.SDWriter(output_sdf)) as writer:
        for mol in suppl:
            if mol is None:
                continue
            with tempfile.TemporaryDirectory() as tmpdir:
                poses, affinities, rmsds, atom_data_dfs = dock_smina(input_sdf_rep, input_sdf_lig, inp, tmpdir)
                if atom_data_dfs:
                    scoring_contributions = [df.drop(columns=["atomid", "el", "pos"]).astype(float).apply(sum, axis=0) for df in atom_data_dfs]
                    pd.DataFrame(scoring_contributions[0]).T.to_csv(smina_output, header=False)
                    atom_data_dfs[0].to_csv(atom_terms_output)
                    
            if rmsds:
                rmsds_df = pd.DataFrame(rmsds, columns=[str(len(rmsds))], index=[i for i in range(1,len(rmsds)+1)]).T
                rmsds_df.to_csv(rmsds_output, header=False)

            docked = write_out_confs(poses, writer, n_poses)
            print("Docked {} poses".format(docked))

            affinities_df = pd.DataFrame(affinities, columns=[str(len(affinities))], index=[i for i in range(1,len(affinities)+1)]).T
            affinities_df.to_csv(affins_output, header=False)


