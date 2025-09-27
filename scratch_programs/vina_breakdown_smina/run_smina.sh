#!/bin/bash

num_complexes=100
complex_dir=/home/andyt/DProjects/DMCTS/VINA_ChemTSv2/data/dataset_benchmarks/HCs/PCCP_top_100_benchmark/3D_pose
outdir=results/PCCP_top_100
desc_script=descs_hydrocarbons.py

rm $outdir/results.csv

get_smina_and_desc_contributions () {
    # SMINA to get VINA scoring contributions
    python3 get_receptor_ligand_for_smina.py $1
    obabel -isdf host.sdf -opdbqt -O host.pdbqt &> /dev/null
    obabel -isdf guest.sdf -opdbqt -O guest.pdbqt &> /dev/null
    smina -r host.pdbqt -l guest.pdbqt --score_only --log score.log &> /dev/null
    res_smina=$(python3 read_smina_out.py score.log)
    if [ -z "${res_smina}" ]; then
	res_smina=",,,,,"
    else
	rm guest.pdbqt host.sdf host.pdbqt score.log
    fi

    # descriptor contributions
    res_desc=$(python3 ${desc_script} guest.sdf)
    if [ -z "${res_desc}" ]; then
	res_desc=",,,,"
    else
	rm guest.sdf
    fi

    echo $res_smina,$res_desc
}

for i in $(seq 0 $((num_complexes-1))); do
    conts=$( get_smina_and_desc_contributions $complex_dir/mol_${i}.sdf)
    echo $i,$conts >> $outdir/results.csv
done

