#!/bin/sh
NORM_ELM2_DATA=/nscratch/xpto.elm2.gz
NORM_IMAGE_PREFIX=xpto_V2
VOXEL_SIZE=2
DISTANCE=194

B=~ricardo/lmrec-trunk

# Run it all in one go!
zcat ${NORM_ELM2_DATA} | \
${B}/elm2todkfz -o -  - | \
${B}/Norm_Total_Gen ${DISTANCE} ${VOXEL_SIZE} - $NORM_IMAGE_PREFIX
