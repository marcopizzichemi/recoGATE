CXXFLAGS ?=-O3 -static

all:  ClearPEM_LMRec Norm_Total_Gen elm2todkfz elm2todkfzA Random_Gen UnitImage

ClearPEM_LMRec: LMRec/ClearPEM_LMRec.cpp LMRec/lmrec_common.cpp LMRec/interfile.cpp LMRec/Metz.cpp lmFormats/EventStore.cpp
	$(CXX) -o $@ LMRec/ClearPEM_LMRec.cpp LMRec/lmrec_common.cpp LMRec/interfile.cpp LMRec/Metz.cpp lmFormats/EventStore.cpp $(CXXFLAGS) -lpthread -lfftw3

Norm_Total_Gen: LMRec/Norm_Total_Gen.cpp LMRec/lmrec_common.cpp
	$(CXX) -o $@ LMRec/Norm_Total_Gen.cpp LMRec/lmrec_common.cpp $(CXXFLAGS) $(LDFLAGS) -lpthread

elm2todkfzA: lmFormats/elm2todkfzA.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

elm2todkfz: lmFormats/elm2todkfz.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

UnitImage: var/UnitImage.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

Random_Gen: LMRec/Random_Gen.cpp LMRec/lmrec_common.cpp LMRec/interfile.cpp LMRec/Metz.cpp lmFormats/EventStore.cpp
	$(CXX) -o $@ LMRec/Random_Gen.cpp LMRec/lmrec_common.cpp LMRec/interfile.cpp LMRec/Metz.cpp lmFormats/EventStore.cpp $(CXXFLAGS) -lpthread -lfftw3

clean:
	rm -f ClearPEM_LMRec Norm_Total_Gen elm2todkfz elm2todkfzA Random_Gen UnitImage

.PHONY: clean all
