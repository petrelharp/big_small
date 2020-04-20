#!/bin/bash

USAGE="
    $0 (paramfile)
"

if [ $# -lt 1 ]
then
    echo "$USAGE"
    exit 1
fi

PARAMFILE=$1

eidos <(echo "source(\"${PARAMFILE}\"); source(\"make_landscape.eidos\");")
