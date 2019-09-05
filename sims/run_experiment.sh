#!/bin/bash

for GAMMA_B in 0.01 0.05 0.1 0.15 0.2 0.25 0.4 0.5 0.7 1.0 1.5 2.0
do
    echo $GAMMA_B
    ./plot_1d_lineage_tree.py 20 3 big_small_1d.slim BURNIN=1000 NUMGENS=5000 W=300 GAMMA_B=$GAMMA_B 
done
