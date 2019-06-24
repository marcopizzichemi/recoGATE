#!/bin/sh
DATA_FILE=$1
shift

B=$(dirname $0)
if [[ -z $B ]]; then B=. ; fi

zcat ${DATA_FILE} | ${B}/elm2todkfz -o -  - | ${B}/ClearPEM_LMRec $*

