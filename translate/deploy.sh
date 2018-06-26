#!/bin/bash
if [ -z "$1" ]
  then
    echo "No argument supplied, so compiling with machine compiler..."
    g++ -o ../bin/simToData simToData.cpp `root-config --cflags --libs`
    g++ -o ../bin/dataset dataset.cpp `root-config --cflags --libs`
    rm ../bin/structDictionary*
    cp structDictionary.C ../bin/
    cp struct.hh ../bin/
    cp start.translate ../bin/
else
  if [ $1 = "lxplus" ]; then
    echo "Sourcing compiler for lxplus..."
    source /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6-gcc49-opt/setup.sh
    source /afs/cern.ch/sw/lcg/external/geant4/10.2.p03/x86_64-slc6-gcc49-opt/bin/geant4.sh
    # source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
    source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
    g++ -o ../bin/simToData simToData.cpp `root-config --cflags --libs`
    g++ -o ../bin/dataset dataset.cpp `root-config --cflags --libs`
    rm ../bin/structDictionary*
    cp structDictionary.C ../bin/
    cp struct.hh ../bin/
    cp start.translate ../bin/
  else
    echo "Invalid argument $1"
  fi
fi
