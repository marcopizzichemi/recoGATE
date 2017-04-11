g++ -o ../build/simToData simToData.cpp `root-config --cflags --libs`
g++ -o ../build/dataset dataset.cpp `root-config --cflags --libs`
cp structDictionary.C ../build/
cp struct.hh ../build/
cp start.reco ../build/
