#!/bin/bash

# check if there are args, otherwise give feedback
if [ $# -eq 0 ]
  then
    echo "USAGE:"
    echo "$0 [binary-folder-path] [input-file] [energy-low] [energy-high] [time-window] [scanner-diameter] [axial-size] [pixel-length] [threads] [sensitivity-file] [output-prefix]"
    exit 1
fi

#give names to args
binfolder=$1
inputFile=$2
eMin=$3
eMax=$4
dt=$5
distance=$6
axial=$7
pixel=$8
threads=$9
sensitivity=${10}
output=${11}

#create pipes
mkfifo cat_pipe
mkfifo lmf_pipe

# create empty scatter.txt file if it's not there
if [[ ! -e scatter.txt ]]; then
    touch scatter.txt
fi

#start chain
pigz -d -c ${inputFile} > cat_pipe &
${binfolder}/elm2todkfz -o lmf_pipe --energy-low ${eMin} --energy-high ${eMax} --time-window ${dt} --scatter-total-events scatter.txt - < cat_pipe &
${binfolder}/ClearPEM_LMRec -i lmf_pipe -d ${distance} --axial-size ${axial} --pixel-length ${pixel} --threads ${threads} --save-iterations 1 -n ${sensitivity} -o ${output}

#remove pipes
rm cat_pipe
rm lmf_pipe
