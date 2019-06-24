#!/bin/bash

export ITK_DIR=/usr/local/lib/cmake/ITK-4.2/
export VTK_DIR=/home/marco/Programs/VTK5.10.1/build/

#export ITK_DIR=/usr/local/lib/cmake/ITK-4.2/
#export VTK_DIR=/home/clearpem/VTK/build/

RUNDIR=$PWD

#echo ">>>>>> Compiling QtRecoPem"

#cd QtRecoPem && qmake . && make 

echo ">>>>>> Compiling LMRec with compiler $CC"
cd $RUNDIR
#make clean
cd lmrec-trunk
make clean 
cd ..
make -j5 -C lmrec-trunk

echo ">>>>>> Compiling corrections"

cd corrections
mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make

echo ">>>>> Deploying all relevant files..."
cd $RUNDIR
cp corrections/build/MeshCreator pem-sonic-tools/
cp corrections/build/CorrectAttenuation_SV pem-sonic-tools/
cp corrections/build/CorrectScatter_multiCore_LMSW pem-sonic-tools/
cp lmrec-trunk/ClearPEM_LMRec pem-sonic-tools/
cp lmrec-trunk/Norm_Total_Gen pem-sonic-tools/
cp lmrec-trunk/Random_Gen pem-sonic-tools/
cp lmrec-trunk/elm2todkfz pem-sonic-tools/
