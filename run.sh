#!/bin/bash

set -x

WORK_FOLDER=data
EXE_FOLDER=.
PREFIX=${WORK_FOLDER}/at_at
#MCSCAN_PARAMS="-A -u 5000"
MCSCAN_PARAMS=""

# Filter redundant blast hits
#./filter_blast.py ${PREFIX}.blast.m8 ${PREFIX}.blast

# Build mcl gene families
more ${PREFIX}.blast | ${EXE_FOLDER}/mcl - --abc --abc-neg-log -abc-tf 'mul(0.4343), ceil(200)' -o ${PREFIX}.mcl

# Run mcscan program
${EXE_FOLDER}/mcscan ${MCSCAN_PARAMS} ${PREFIX}

