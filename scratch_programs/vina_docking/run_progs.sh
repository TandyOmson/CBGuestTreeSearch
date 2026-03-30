#!/bin/bash

hostfile=$PWD/data/cb7.sdf
vinainp=$PWD/vinadock.inp

outdir=/home/andyt/DProjects/Dsingle_system_binding/Ddft_descriptor_model/dft/Dsig_mols_3/docked

for i in {1..4}; do
    guestfile=/home/andyt/DProjects/Dsingle_system_binding/Ddft_descriptor_model/dft/Dsig_mols_3/pubchem/mol_$i.sdf
    python3 dock_vina.py -n 1 -r $hostfile -l $guestfile -o $outdir/mol_$i.sdf -i $vinainp -a $outdir/affins.csv
done
