#!/bin/bash


USAGE="
    $0 [.trees file [.trees file]]
"

if [ $# -lt 1 ]
then
    echo "$USAGE"
    exit 1
fi

for treefile in "$@"
do
    SCRIPT=$(dirname $treefile).slim
    ARGS="$(echo "$treefile" | cut -f 4- -d '_' | sed -e 's/.trees//' | sed -e 's/_\([0-9.]\+\)/=\1 /g' | sed -e 's/ _/ /g')"
    ./plot_1d_lineage_tree.py 100 3 $SCRIPT $ARGS
done
