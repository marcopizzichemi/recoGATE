g++ -o ../bin/simToData simToData.cpp `root-config --cflags --libs`
g++ -o ../bin/dataset dataset.cpp `root-config --cflags --libs`
cp structDictionary.C ../bin/
cp struct.hh ../bin/
cp start.translate ../bin/
