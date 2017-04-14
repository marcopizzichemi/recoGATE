// compile with
// g++ -o simToData simToData.cpp `root-config --cflags --libs`

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>    // std::sort
#include <getopt.h>
#include "TROOT.h"
#include "TRandom1.h"
#include "TSystem.h"
#include "TChain.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH1.h"
#include "TAxis.h"
#include "TDirectory.h"
#include "TList.h"
#include "Rtypes.h"
#include "TChainElement.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TRandom.h"




//struct of energy depositions events.
// where (crystal id and position in the crystal) was the energy deposited, how much and when
// local coordinates: remember that in GATE simulation,
struct enDep
{
  Float_t  localPosX;
  Float_t  localPosY;
  Float_t  localPosZ;
  Int_t    parentID;
  Int_t    trackID;
  Int_t    globalCryID;
  Int_t    eventID;
  Int_t    gantryID;
  Int_t    rsectorID;
  Int_t    moduleID;
  Int_t    submoduleID;
  Int_t    crystalID;
  Int_t    layerID;
  Int_t    sourceID;
  Double_t time;
  Float_t  edep;
};

struct avgCryEnergyDep
{
  Int_t eventID;
  Int_t parentID;
  Int_t trackID;
  int   globalCryID;
  Int_t gantryID;
  Int_t rsectorID;
  Int_t moduleID;
  Int_t submoduleID;
  Int_t crystalID;
  Int_t layerID;
  Int_t sourceID;
  float x;
  float y;
  float z;
  float sx;
  float sy;
  float sz;
  float energy;
  float time;
};

struct point
{
  Int_t eventID;
  Int_t parentID;
  Int_t trackID;
  Int_t sourceID;
  Int_t crystalIDforSTIR;
  Int_t ringIDforSTIR;
  float x;
  float y;
  float z;
  float energy;
  float time;
};

//function to compare deposition event struct vectors using the field time
bool compareByTime(const enDep &a,const enDep  &b)
{
  return a.time < b.time;
}

//smear the input value given a FWHM (abs value, i.e. NOT the percentage)
Float_t gaussianSmear(Double_t mean,Double_t fwhm)
{
  Double_t sigma = fwhm/2.355;
  TRandom *r1 = new TRandom1();
  Double_t result = r1->Gaus(mean,sigma);
  return result;
}

int main(int argc, char** argv)
{

  // if(argc<2) {
  //   std::cout<<"You need at least to provide a directory (where your .root files, output by GATE, are) "<<std::endl ;
  //   std::cout<<"USAGE:"<<std::endl ;
  //   std::cout<<"gateToString input.root [output.root]"<<std::endl ;
  //   return 1;
  // }

  //ARGUMENTS WITH DEFAULTS
  float rmin = 100;
  float crylength = 15;
  //rsector
  int repsec = 2;
  // float ang = 3.1415*2/repsec;
  //module
  int repmodx = 16;
  int repmody = 6;
  int repmodz = 1;
  //array
  float arraymodx = 13;
  float arraymody = 13;
  float arraymodz = 0;
  //submod
  //repetitions
  int repsubx = 4;
  int repsuby = 4;
  int repsubz = 1;
  //positions
  float arraysubx = 3.2;
  float arraysuby = 3.2;
  float arraysubz = 0;
  //crystal
  //repetitions
  int repcryx = 2;
  int repcryy = 2;
  int repcryz = 1;
  //positions
  float arraycryx = 1.6;
  float arraycryy = 1.6;
  float arraycryz = 0;
  //detector parameters
  int smearedDoi = 1;
  int smearedEnergy = 1;
  Double_t energyResolutionFWHM = 0.12; //12% FWHM @511, for the moment we take the same also for low energy...
  Double_t doiResolutionFWHM = 2.8; //2.8mm doi res measured on the module

  char* inputName = NULL;
  char* outputName = NULL;
  // float pixelLength;

  //----------------------------------------------------------------------//
  //                                                                      //
  //                        ROOT STUFF                                    //
  //                                                                      //
  //----------------------------------------------------------------------//
  gROOT->ProcessLine("#include <vector>"); //needed by ROOT to deal with standard vectors
  //HACK to use the dictionary easily
  // Code taken from: http://www.gamedev.net/community/forums/topic.asp?topic_id=459511
  std::string fullFileName = "";
  std::string path = "";
  pid_t pid = getpid();
  char buf[20] = {0};
  sprintf(buf,"%d",pid);
  std::string _link = "/proc/";
  _link.append( buf );
  _link.append( "/exe");
  char proc[512];
  int ch = readlink(_link.c_str(),proc,512);
  if (ch != -1) {
    proc[ch] = 0;
    path = proc;
    std::string::size_type t = path.find_last_of("/");
    path = path.substr(0,t);
  }
  fullFileName = path + std::string("/");
  std::string command = ".L " + fullFileName + "structDictionary.C+";
  gROOT->ProcessLine(command.c_str());

  static struct option longOptions[] =
  {
			{ "inner-radius", required_argument, 0, 0 },
			{ "crystal-length", required_argument, 0, 0 },
			{ "radial-sectors", required_argument, 0, 0 },
			{ "module-repetions-x", required_argument, 0, 0 },
			{ "module-repetions-y", required_argument, 0, 0 },
			{ "module-repetions-z", required_argument, 0, 0 },
      { "module-distance-x", required_argument, 0, 0 },
			{ "module-distance-y", required_argument, 0, 0 },
			{ "module-distance-z", required_argument, 0, 0 },
      { "submodule-repetions-x", required_argument, 0, 0 },
			{ "submodule-repetions-y", required_argument, 0, 0 },
			{ "submodule-repetions-z", required_argument, 0, 0 },
      { "submodule-distance-x", required_argument, 0, 0 },
			{ "submodule-distance-y", required_argument, 0, 0 },
			{ "submodule-distance-z", required_argument, 0, 0 },
      { "crystals-repetions-x", required_argument, 0, 0 },
			{ "crystals-repetions-y", required_argument, 0, 0 },
			{ "crystals-repetions-z", required_argument, 0, 0 },
      { "crystals-distance-x", required_argument, 0, 0 },
			{ "crystals-distance-y", required_argument, 0, 0 },
			{ "crystals-distance-z", required_argument, 0, 0 },
      { "energy-resolution-fwhm", required_argument, 0, 0 },
      { "doi-resolution-fwhm", required_argument, 0, 0 },
      { "smeared-energy", required_argument, 0, 0 },
      { "smeared-doi", required_argument, 0, 0},
			{ NULL, 0, 0, 0 }
	};

  while(1) {

		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:", longOptions, &optionIndex);

		if (c == -1) {
			break;
		}

		if (c == 'i'){
			inputName = (char *)optarg;
    }
		else if (c == 'o'){
      outputName = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      rmin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 1){
      crylength = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 2){
      repsec = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 3){
      repmodx = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      repmody = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 5){
      repmodz = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      arraymodx = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      arraymody = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      arraymodz = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      repsubx = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      repsuby = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      repsubz = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 12){
      arraysubx = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 13){
      arraysuby = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 14){
      arraysubz = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 15){
      repcryx = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 16){
      repcryy = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 17){
      repcryz = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 18){
      arraycryx = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 19){
      arraycryy = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 20){
      arraycryz = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 21){
      energyResolutionFWHM = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 22){
      doiResolutionFWHM = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 23){
      smearedEnergy = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 24){
      smearedDoi = atoi((char *)optarg);
    }

		else { //FIXME
			std::cout	<< "Usage: " << argv[0] << std::endl
			<< "\t\t" << "[ -i <input file> ] " << std::endl
			<< "\t\t" << "[ -o <output file> ] " << std::endl
			<< "\t\t" << "[ -d | --crystal-distance <crystal distance in mm> ] " << std::endl
			<< "\t\t" << "[ -n | --norm-prefix <norm files prefix> ] " << std::endl
			<< "\t\t" << "[ --pixel-length <voxel size in mm> ] " << std::endl
			<< "\t\t" << "[ --iterations <number of iterations> ] " << std::endl
			<< "\t\t" << "[ --osem <0|1> ] " << std::endl
			<< "\t\t" << "[ --filter-fwhm <0 to disable or FWHM in mm> ] " << std::endl
			<< "\t\t" << "[ --filter-power < Metz power, 0 for gaussian ] " << std::endl
			<< "\t\t" << "[ --filter-iter <filter each n iterations> ]" << std::endl
			<< "\t\t"       << "[ --no-interfile ]" << std::endl
			<< "\t\t" << "[ --partial <event fraction %>" << std::endl
			<< "\t\t" << "[ --rand-image <random file> ]" << std::endl
			<< "\t\t" << "[ --rand-cut <random cut value, 0 to 0.5> ]" << std::endl
			<< "\t\t" << "[ --save-iterations <N> ] " << std::endl
			<< "\t\t" << std::endl;

			return 1;
		}
	}

  //CHECK PARAMETERS
  if(inputName == NULL | outputName == NULL)
  {
    std::cout	<< "ERROR: you need to provide at least input and output file!!!" << std::endl;
    return 1;
  }

  std::cout << "SIM TO DATA -----------------------------------"<< std::endl;
  std::cout << "Input file name      = " << inputName           << std::endl
            << "Output file name     = " << outputName          << std::endl
            << "rmin                 = " << rmin                << std::endl
            << "crylength            = " << crylength           << std::endl
            << "repsec               = " << repsec              << std::endl
            << "repmodx              = " << repmodx             << std::endl
            << "repmody              = " << repmody             << std::endl
            << "repmodz              = " << repmodz             << std::endl
            << "arraymodx            = " << arraymodx           << std::endl
            << "arraymody            = " << arraymody           << std::endl
            << "arraymodz            = " << arraymodz           << std::endl
            << "repsubx              = " << repsubx             << std::endl
            << "repsuby              = " << repsuby             << std::endl
            << "repsubz              = " << repsubz             << std::endl
            << "arraysubx            = " << arraysubx           << std::endl
            << "arraysuby            = " << arraysuby           << std::endl
            << "arraysubz            = " << arraysubz           << std::endl
            << "repcryx              = " << repcryx             << std::endl
            << "repcryy              = " << repcryy             << std::endl
            << "repcryz              = " << repcryz             << std::endl
            << "arraycryx            = " << arraycryx           << std::endl
            << "arraycryy            = " << arraycryy           << std::endl
            << "arraycryz            = " << arraycryz           << std::endl
            << "smearedDoi           = " << smearedDoi          << std::endl
            << "smearedEnergy        = " << smearedEnergy       << std::endl
            << "energyResolutionFWHM = " << energyResolutionFWHM<< std::endl
            << "doiResolutionFWHM    = " << doiResolutionFWHM   << std::endl
            <<  std::endl;



  std::string filedir, inputfilename ;
  std::string filename, Moutputfilename, Poutputfilename ;

  inputfilename = inputName ;
  // inputfilename = filedir + "*.root" ;

  std::cout << "Input file name is " << inputfilename << std::endl;
  TChain *Hits = new TChain("Hits") ;
  TChain *Singles = new TChain("Singles") ;
  TChain *Coincidences = new TChain("Coincidences") ;
  Hits->Add(inputfilename.c_str());
  Singles->Add(inputfilename.c_str());
  Coincidences->Add(inputfilename.c_str());
  long int Trues = 0;
  long int Scatters = 0;
  long int Randoms = 0;
  // Int_t   nbytes = 0;

  std::string outputFileName = outputName;

  // if(argc > 2) outputFileName = argv[2];
  TFile* outputFile = new TFile(outputFileName.c_str(),"recreate");

  TTree* outputTree = new TTree("test","test");
  // avgCryEnergyDep outAvg;
  std::vector<avgCryEnergyDep> averageDepEvents;
  std::vector<point> points;
  outputTree->Branch("points","std::vector<point>",&points);


  //######################################################################################
  //#                        Variables                                                   #
  //######################################################################################

  //HITS
  Int_t HITSPDGEncoding;
  Int_t HITStrackID;
  Int_t HITSparentID;
  Double_t HITStime;
  Float_t HITSedep;
  Float_t HITSstepLength;
  Float_t HITSposX;
  Float_t HITSposY;
  Float_t HITSposZ;
  Float_t HITSlocalPosX;
  Float_t HITSlocalPosY;
  Float_t HITSlocalPosZ;
  Int_t HITSgantryID;
  Int_t HITSrsectorID;
  Int_t HITSmoduleID;
  Int_t HITSsubmoduleID;
  Int_t HITScrystalID;
  Int_t HITSlayerID;
  Int_t HITSphotonID;
  Int_t HITSnPhantomCompton;
  Int_t HITSnCrystalCompton;
  Int_t HITSnPhantomRayleigh;
  Int_t HITSnCrystalRayleigh;
  Int_t HITSprimaryID;
  Float_t HITSsourcePosX;
  Float_t HITSsourcePosY;
  Float_t HITSsourcePosZ;
  Int_t HITSsourceID;
  Int_t HITSeventID;
  Int_t HITSrunID;
  Float_t HITSaxialPos;
  Float_t HITSrotationAngle;
  Int_t HITSvolumeID;
  Char_t HITSprocessName[40];
  Char_t HITScomptVolName[40];
  Char_t HITSRayleighVolName[40];

  //SINGLES
  Int_t SINGLESrunID;
  Int_t SINGLESeventID;
  Int_t SINGLESsourceID;
  Float_t SINGLESsourcePosX;
  Float_t SINGLESsourcePosY;
  Float_t SINGLESsourcePosZ;
  Float_t SINGLEStime;
  Float_t SINGLESenergy;
  Float_t SINGLESglobalPosX;
  Float_t SINGLESglobalPosY;
  Float_t SINGLESglobalPosZ;
  Int_t SINGLESgantryID;
  Int_t SINGLESrsectorID;
  Int_t SINGLESmoduleID;
  Int_t SINGLESsubmoduleID;
  Int_t SINGLEScrystalID;
  Int_t SINGLESlayerID;
  Int_t SINGLEScomptonPhantom;
  Int_t SINGLEScomptonCrystal;
  Int_t SINGLESRayleighPhantom;
  Int_t SINGLESRayleighCrystal;
  Float_t SINGLESaxialPos;
  Float_t SINGLESrotationAngle;
  Char_t SINGLEScomptVolName[40];
  Char_t SINGLESRayleighVolName[40];

  //COINCIDENCES
  Float_t         CoincAxialPos, CoincRotationAngle, sinogramS, sinogramTheta;
  Char_t          comptVolName1[40], comptVolName2[40];
  Char_t          RayleighVolName1[40], RayleighVolName2[40];
  Int_t           comptonPhantom1, comptonPhantom2;
  Int_t           comptonCrystal1, comptonCrystal2;
  Int_t           RayleighPhantom1, RayleighPhantom2;
  Int_t           RayleighCrystal1, RayleighCrystal2;
  Int_t           CoincRunID, sourceID1, sourceID2, eventID1, eventID2;
  Int_t           layerID1, layerID2, crystalID1, crystalID2;
  Int_t           submoduleID1, submoduleID2, moduleID1, moduleID2, rsectorID1, rsectorID2, gantryID1, gantryID2;
  Float_t         energy1, energy2;
  Float_t         globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2;
  Float_t         sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2;
  Double_t        time1, time2;

  //######################################################################################
  //#                        Set branch addresses                                        #
  //######################################################################################

  //TTree Hits
  Hits->SetBranchStatus("*",0); //WHY?
  Hits->SetBranchAddress("PDGEncoding",&HITSPDGEncoding);
  Hits->SetBranchAddress("trackID",&HITStrackID);
  Hits->SetBranchAddress("parentID",&HITSparentID);
  Hits->SetBranchAddress("time",&HITStime);
  Hits->SetBranchAddress("edep",&HITSedep);
  Hits->SetBranchAddress("stepLength",&HITSstepLength);
  Hits->SetBranchAddress("posX",&HITSposX);
  Hits->SetBranchAddress("posY",&HITSposY);
  Hits->SetBranchAddress("posZ",&HITSposZ);
  Hits->SetBranchAddress("localPosX",&HITSlocalPosX);
  Hits->SetBranchAddress("localPosY",&HITSlocalPosY);
  Hits->SetBranchAddress("localPosZ",&HITSlocalPosZ);
  Hits->SetBranchAddress("gantryID",&HITSgantryID);
  Hits->SetBranchAddress("rsectorID",&HITSrsectorID);
  Hits->SetBranchAddress("moduleID",&HITSmoduleID);
  Hits->SetBranchAddress("submoduleID",&HITSsubmoduleID);
  Hits->SetBranchAddress("crystalID",&HITScrystalID);
  Hits->SetBranchAddress("layerID",&HITSlayerID);
  Hits->SetBranchAddress("photonID",&HITSphotonID);
  Hits->SetBranchAddress("nPhantomCompton",&HITSnPhantomCompton);
  Hits->SetBranchAddress("nCrystalCompton",&HITSnCrystalCompton);
  Hits->SetBranchAddress("nPhantomRayleigh",&HITSnPhantomRayleigh);
  Hits->SetBranchAddress("nCrystalRayleigh",&HITSnCrystalRayleigh);
  Hits->SetBranchAddress("primaryID",&HITSprimaryID);
  Hits->SetBranchAddress("sourcePosX",&HITSsourcePosX);
  Hits->SetBranchAddress("sourcePosY",&HITSsourcePosY);
  Hits->SetBranchAddress("sourcePosZ",&HITSsourcePosZ);
  Hits->SetBranchAddress("sourceID",&HITSsourceID);
  Hits->SetBranchAddress("eventID",&HITSeventID);
  Hits->SetBranchAddress("runID",&HITSrunID);
  Hits->SetBranchAddress("axialPos",&HITSaxialPos);
  Hits->SetBranchAddress("rotationAngle",&HITSrotationAngle);
  Hits->SetBranchAddress("volumeID",&HITSvolumeID);
  Hits->SetBranchAddress("processName",&HITSprocessName);
  Hits->SetBranchAddress("comptVolName",&HITScomptVolName);
  Hits->SetBranchAddress("RayleighVolName",&HITSRayleighVolName);

  //TTree Singles
  Singles->SetBranchStatus("*",0); //WHY?
  Singles->SetBranchAddress("runID",&SINGLESrunID);
  Singles->SetBranchAddress("eventID",&SINGLESeventID);
  Singles->SetBranchAddress("sourceID",&SINGLESsourceID);
  Singles->SetBranchAddress("sourcePosX",&SINGLESsourcePosX);
  Singles->SetBranchAddress("sourcePosY",&SINGLESsourcePosY);
  Singles->SetBranchAddress("sourcePosZ",&SINGLESsourcePosZ);
  Singles->SetBranchAddress("time",&SINGLEStime);
  Singles->SetBranchAddress("energy",&SINGLESenergy);
  Singles->SetBranchAddress("globalPosX",&SINGLESglobalPosX);
  Singles->SetBranchAddress("globalPosY",&SINGLESglobalPosY);
  Singles->SetBranchAddress("globalPosZ",&SINGLESglobalPosZ);
  Singles->SetBranchAddress("gantryID",&SINGLESgantryID);
  Singles->SetBranchAddress("rsectorID",&SINGLESrsectorID);
  Singles->SetBranchAddress("moduleID",&SINGLESmoduleID);
  Singles->SetBranchAddress("submoduleID",&SINGLESsubmoduleID);
  Singles->SetBranchAddress("crystalID",&SINGLEScrystalID);
  Singles->SetBranchAddress("layerID",&SINGLESlayerID);
  Singles->SetBranchAddress("comptonPhantom",&SINGLEScomptonPhantom);
  Singles->SetBranchAddress("comptonCrystal",&SINGLEScomptonCrystal);
  Singles->SetBranchAddress("RayleighPhantom",&SINGLESRayleighPhantom);
  Singles->SetBranchAddress("RayleighCrystal",&SINGLESRayleighCrystal);
  Singles->SetBranchAddress("axialPos",&SINGLESaxialPos);
  Singles->SetBranchAddress("rotationAngle",&SINGLESrotationAngle);
  Singles->SetBranchAddress("comptVolName",&SINGLEScomptVolName);
  Singles->SetBranchAddress("RayleighVolName",&SINGLESRayleighVolName);

  //TTree Coincidences
  Coincidences->SetBranchStatus("*",0); //WHY?
  Coincidences->SetBranchAddress("axialPos",&CoincAxialPos);
  Coincidences->SetBranchAddress("comptVolName1",&comptVolName1);
  Coincidences->SetBranchAddress("comptVolName2",&comptVolName2);
  Coincidences->SetBranchAddress("comptonCrystal1",&comptonCrystal1);
  Coincidences->SetBranchAddress("comptonCrystal2",&comptonCrystal2);
  Coincidences->SetBranchAddress("crystalID1",&crystalID1);
  Coincidences->SetBranchAddress("crystalID2",&crystalID2);
  Coincidences->SetBranchAddress("comptonPhantom1",&comptonPhantom1);
  Coincidences->SetBranchAddress("comptonPhantom2",&comptonPhantom2);
  Coincidences->SetBranchAddress("RayleighVolName1",&RayleighVolName1);
  Coincidences->SetBranchAddress("RayleighVolName2",&RayleighVolName2);
  Coincidences->SetBranchAddress("RayleighPhantom1",&RayleighPhantom1);
  Coincidences->SetBranchAddress("RayleighPhantom2",&RayleighPhantom2);
  Coincidences->SetBranchAddress("RayleighCrystal1",&RayleighCrystal1);
  Coincidences->SetBranchAddress("RayleighCrystal2",&RayleighCrystal2);
  Coincidences->SetBranchAddress("energy1",&energy1);
  Coincidences->SetBranchAddress("energy2",&energy2);
  Coincidences->SetBranchAddress("eventID1",&eventID1);
  Coincidences->SetBranchAddress("eventID2",&eventID2);
  Coincidences->SetBranchAddress("globalPosX1",&globalPosX1);
  Coincidences->SetBranchAddress("globalPosX2",&globalPosX2);
  Coincidences->SetBranchAddress("globalPosY1",&globalPosY1);
  Coincidences->SetBranchAddress("globalPosY2",&globalPosY2);
  Coincidences->SetBranchAddress("globalPosZ1",&globalPosZ1);
  Coincidences->SetBranchAddress("globalPosZ2",&globalPosZ2);
  Coincidences->SetBranchAddress("gantryID1",&gantryID1);
  Coincidences->SetBranchAddress("gantryID2",&gantryID2);
  Coincidences->SetBranchAddress("layerID1",&layerID1);
  Coincidences->SetBranchAddress("layerID2",&layerID2);
  Coincidences->SetBranchAddress("moduleID1",&moduleID1);
  Coincidences->SetBranchAddress("moduleID2",&moduleID2);
  Coincidences->SetBranchAddress("rotationAngle",&CoincRotationAngle);
  Coincidences->SetBranchAddress("rsectorID1",&rsectorID1);
  Coincidences->SetBranchAddress("rsectorID2",&rsectorID2);
  Coincidences->SetBranchAddress("runID",&CoincRunID);
  Coincidences->SetBranchAddress("sinogramS",&sinogramS);
  Coincidences->SetBranchAddress("sinogramTheta",&sinogramTheta);
  Coincidences->SetBranchAddress("sourceID1",&sourceID1);
  Coincidences->SetBranchAddress("sourceID2",&sourceID2);
  Coincidences->SetBranchAddress("sourcePosX1",&sourcePosX1);
  Coincidences->SetBranchAddress("sourcePosX2",&sourcePosX2);
  Coincidences->SetBranchAddress("sourcePosY1",&sourcePosY1);
  Coincidences->SetBranchAddress("sourcePosY2",&sourcePosY2);
  Coincidences->SetBranchAddress("sourcePosZ1",&sourcePosZ1);
  Coincidences->SetBranchAddress("sourcePosZ2",&sourcePosZ2);
  Coincidences->SetBranchAddress("submoduleID1",&submoduleID1);
  Coincidences->SetBranchAddress("submoduleID2",&submoduleID2);
  Coincidences->SetBranchAddress("time1",&time1);
  Coincidences->SetBranchAddress("time2",&time2);

  Int_t HITSnentries = (Int_t)(Hits->GetEntries());
  Int_t SINGLESCnentries = (Int_t)(Singles->GetEntries());
  Int_t COINCnentries = (Int_t)(Coincidences->GetEntries());
  std::cout << "Number of hits         = " << HITSnentries << std::endl;
  std::cout << "Number of singles      = " << SINGLESCnentries << std::endl;
  std::cout << "Number of coincidences = " << COINCnentries << std::endl;

  Int_t eventsCheck = -1;
  std::vector<enDep> energyDeposition;



  //######################################################################################
  //#                        Parameters                                                  #
  //######################################################################################
  // TODO find a way to make them external
  // gate geometry
  // float rmin = 100;
  // float crylength = 15;
  // //rsector
  // int repsec = 2;
  float ang = 3.1415*2/repsec;
  // //module
  // int repmodx = 16;
  // int repmody = 6;
  // int repmodz = 1;
  //
  //
  // float arraymodx = 13;
  // float arraymody = 13;
  // float arraymodz = 1;
  //
  // //submod
  // int repsubx = 4;
  // int repsuby = 4;
  // int repsubz = 1;
  //
  // float arraysubx = 3.2;
  // float arraysuby = 3.2;
  // float arraysubz = 0;
  // //crystal
  //
  // int repcryx = 2;
  // int repcryy = 2;
  // int repcryz = 0;
  //
  // float arraycryx = 1.6;
  // float arraycryy = 1.6;
  // float arraycryz = 0;
  //
  //variables for global id
  Int_t nRSectors          = repsec;
  Int_t nModulesXRSector   = repmodx*repmody*repmodz;
  Int_t nSubmodulesXModule = repsubx*repsubz*repsuby;
  Int_t nCrystalXSubmodule  = repcryx*repcryz*repcryy;
  Int_t nLayersXCrystal    = 1;
  //
  // //detector parameters
  // bool smeared = true;
  // Double_t energyResolutionFWHM = 0.12; //12% FWHM @511, for the moment we take the same also for low energy...
  // Double_t doiResolutionFWHM = 2.8; //2.8mm doi res measured on the module

  // for (Int_t i = 0 ; i < (Int_t)(Hits->GetEntries()) ; i++)
  // {
  //   Hits->GetEntry(i); // each i is a energyDeposition
  //   if(HITStrackID > 3) std::cout << i  <<  " "  << HITStrackID  << " " << HITSparentID <<  std::endl;
  // }
  // std::cout << std::endl;


  Int_t eventCounter = 0;
  int statuscounter = 0;
  for (Int_t i = 0 ; i < (Int_t)(Hits->GetEntries()) ; i++)
  // for (Int_t i = 1244500 ; i < 1244611 ; i++)
  // for (Int_t i = 0 ; i < 40 ; i++)
  {
    Hits->GetEntry(i); // each i is a energyDeposition
    // std::cout<< eventsCheck << " " << HITSeventID << " " << HITStrackID <<  std::endl;

    if(eventsCheck != HITSeventID)
    {
      int CurrentID = -1;
      // "average" what happened into each crystal
      // stackingaction in g4 is last-in-first-out so we need to sort energyDeposition before we check what crystal was hit first (needed in gATE?)
      // check crystal sequence, check for returning gammas, calculate an average xyz per crystal

      std::vector < enDep > CrystalEnDepCollection;
      std::vector<std::vector < enDep > > separatedEnDep;
      if(energyDeposition.size() != 0)
      {
        std::sort(energyDeposition.begin(), energyDeposition.end(), compareByTime);
        for(int eEvent = 0; eEvent < energyDeposition.size(); eEvent++) //run on energy depositions and find in how many crystals energy was deposited
        {
          // std::cout << energyDeposition[eEvent].eventID << " " << energyDeposition[eEvent].globalCryID << " " <<  energyDeposition[eEvent].edep << " " << energyDeposition[eEvent].localPosX << " " << energyDeposition[eEvent].time << std::endl;
          //this for cycles on all the endep events. it stores them (each one is a struct) in a std::vector of structs, until the gamma changes crystal
          //when we enter a new crystal, the std::vector of that crystal is pushed_back into (guess what?) a std::vector of std::vector and another std::vector is created
          // so to sumarize:
          // energy deposition event -> struct
          // all energy deposition events in a given crystal -> std::vector<struct> -> i call this a "collection" below
          // all crystals touched by the gamma -> std::vector< std:: vector <struct> > this would be the std::vector of "collections"


          //read the crystal where energy was deposited
          //discard if there was no energy deposition (if edep == 0)
          // if(energyDeposition[eEvent].edep!=0)
          // {
            int cry = energyDeposition[eEvent].globalCryID;
            if(eEvent == 0) CurrentID = cry; //needed for the first step
            //create temp enDep variable and copy this eEvent into it
            enDep tempCrystalEnDep;
            tempCrystalEnDep = energyDeposition[eEvent];
            if(cry != CurrentID) // if this endep event is happening in a new crystal wrt the one before
            {
              separatedEnDep.push_back(CrystalEnDepCollection); // save the collection of this crystal into the std::vector of collections
              CrystalEnDepCollection.clear(); //clear this collection
              CurrentID = cry; //change the current id
            }

            CrystalEnDepCollection.push_back(tempCrystalEnDep); // save this enDep event into the collection of this crystal
            if(eEvent == energyDeposition.size() -1)
            {
              separatedEnDep.push_back(CrystalEnDepCollection);
            }
          // }
        }


        // now for each crystal average what happened inside
        // what happened in a crystal cannot be seen in split details by the detector so we need some useful
        // variables to gain the chance to filter the simulation dataset, if we want to understand
        // what's going on
        // std::vector<avgCryEnergyDep> averageDepEvents;
        averageDepEvents.clear();
        points.clear();
        // now the en dep events are collected by crystal
        // run on each collection and find average
        for(int iColl = 0 ; iColl < separatedEnDep.size(); iColl++)
        {
          avgCryEnergyDep tempAvgEnDep;
          // std::cout << separatedEnDep.at(iColl).size() << std::endl;

          //initialize
          tempAvgEnDep.eventID = separatedEnDep.at(iColl).at(0).eventID;
          tempAvgEnDep.parentID = separatedEnDep.at(iColl).at(0).parentID;
          tempAvgEnDep.trackID = separatedEnDep.at(iColl).at(0).trackID;
          tempAvgEnDep.globalCryID = separatedEnDep.at(iColl).at(0).globalCryID;
          tempAvgEnDep.gantryID = separatedEnDep.at(iColl).at(0).gantryID;
          tempAvgEnDep.rsectorID = separatedEnDep.at(iColl).at(0).rsectorID;
          tempAvgEnDep.moduleID = separatedEnDep.at(iColl).at(0).moduleID;
          tempAvgEnDep.submoduleID = separatedEnDep.at(iColl).at(0).submoduleID;
          tempAvgEnDep.crystalID = separatedEnDep.at(iColl).at(0).crystalID;
          tempAvgEnDep.layerID = separatedEnDep.at(iColl).at(0).layerID;
          tempAvgEnDep.sourceID = separatedEnDep.at(iColl).at(0).sourceID;
          tempAvgEnDep.x = 0;
          tempAvgEnDep.y = 0;
          tempAvgEnDep.z = 0;
          tempAvgEnDep.time = 0;
          tempAvgEnDep.energy = 0;
          tempAvgEnDep.sx = 0;
          tempAvgEnDep.sy = 0;
          tempAvgEnDep.sz = 0;
          //now run on the energy deposition events in this crystal and calculate the averages
          // std::cout << "aaaaaaaaaaaaaa" << std::endl;

          for(int iEndep = 0; iEndep < separatedEnDep.at(iColl).size(); iEndep++)
          {
            tempAvgEnDep.energy += separatedEnDep.at(iColl).at(iEndep).edep;
          }
          if(tempAvgEnDep.energy != 0)
          {
            for(int iEndep = 0; iEndep < separatedEnDep.at(iColl).size(); iEndep++)
            {
              if(separatedEnDep.at(iColl).at(iEndep).edep != 0)
              {
                tempAvgEnDep.x += separatedEnDep.at(iColl).at(iEndep).localPosX * separatedEnDep.at(iColl).at(iEndep).edep;
                tempAvgEnDep.y += separatedEnDep.at(iColl).at(iEndep).localPosY * separatedEnDep.at(iColl).at(iEndep).edep;
                tempAvgEnDep.z += separatedEnDep.at(iColl).at(iEndep).localPosZ * separatedEnDep.at(iColl).at(iEndep).edep;
                tempAvgEnDep.time += separatedEnDep.at(iColl).at(iEndep).time * separatedEnDep.at(iColl).at(iEndep).edep;
              }
            }
            tempAvgEnDep.x = tempAvgEnDep.x / tempAvgEnDep.energy;
            tempAvgEnDep.y = tempAvgEnDep.y / tempAvgEnDep.energy;
            tempAvgEnDep.z = tempAvgEnDep.z / tempAvgEnDep.energy;
            tempAvgEnDep.time = tempAvgEnDep.time / tempAvgEnDep.energy;
            float varx = 0.0; // variance (needed for stdev afterwards)
            float vary = 0.0;
            float varz = 0.0;
            for(int iEndep = 0; iEndep < separatedEnDep.at(iColl).size(); iEndep++)
            {
              varx += (separatedEnDep.at(iColl).at(iEndep).edep * pow(separatedEnDep.at(iColl).at(iEndep).localPosX  - tempAvgEnDep.x,2)) / tempAvgEnDep.energy;
              vary += (separatedEnDep.at(iColl).at(iEndep).edep * pow(separatedEnDep.at(iColl).at(iEndep).localPosY  - tempAvgEnDep.y,2)) / tempAvgEnDep.energy;
              varz += (separatedEnDep.at(iColl).at(iEndep).edep * pow(separatedEnDep.at(iColl).at(iEndep).localPosZ  - tempAvgEnDep.z,2)) / tempAvgEnDep.energy;
            }
            tempAvgEnDep.sx = sqrt(varx);
            tempAvgEnDep.sy = sqrt(vary);
            tempAvgEnDep.sz = sqrt(varz);

            //save into the std::vector of averages

            averageDepEvents.push_back(tempAvgEnDep);
          }
        }

        // now output the reading of this event, like this
        //calculate hit position, for each crystal

        for(int iAvg = 0; iAvg < averageDepEvents.size(); iAvg++)
        {
          //coordinates from crytal hit
          // float ymod = ( (averageDepEvents[iAvg].moduleID % repmody) + 0.5 ) * (arraymody) - arraymody*(repmody/2.0);
          // float zmod = ( (averageDepEvents[iAvg].moduleID / repmody) + 0.5 ) * (arraymodz) - arraymodz*(repmodz/2.0);
          // float ysub = ( (averageDepEvents[iAvg].submoduleID % repsuby) + 0.5 ) * (arraysuby) - arraysuby*(repsuby/2.0);
          // float zsub = ( (averageDepEvents[iAvg].submoduleID / repsuby) + 0.5 ) * (arraysubz) - arraysubz*(repsubz/2.0);
          // float ycry = ( (averageDepEvents[iAvg].crystalID % repcryy) + 0.5 ) * (arraycryy) - arraycryy*(repcryy/2.0);
          // float zcry = ( (averageDepEvents[iAvg].crystalID / repcryy) + 0.5 ) * (arraycryz) - arraycryz*(repcryz/2.0);

          float ymod = ( (averageDepEvents[iAvg].moduleID % repmody) + 0.5 ) * (arraymody) - arraymody*(repmody/2.0);
          float xmod = ( (averageDepEvents[iAvg].moduleID / repmody) + 0.5 ) * (arraymodx) - arraymodx*(repmodx/2.0);
          float ysub = ( (averageDepEvents[iAvg].submoduleID % repsuby) + 0.5 ) * (arraysuby) - arraysuby*(repsuby/2.0);
          float xsub = ( (averageDepEvents[iAvg].submoduleID / repsuby) + 0.5 ) * (arraysubx) - arraysubx*(repsubx/2.0);
          float ycry = ( (averageDepEvents[iAvg].crystalID % repcryy) + 0.5 ) * (arraycryy) - arraycryy*(repcryy/2.0);
          float xcry = ( (averageDepEvents[iAvg].crystalID / repcryy) + 0.5 ) * (arraycryx) - arraycryx*(repcryx/2.0);
          //calculate x and y in FOV BEFORE the rsectorID rotation on the basis of their module id, submodule etc..
          float yabs = ymod+ysub+ycry;
          float xabs = xmod+xsub+xcry;
          // CAREFUL!!!!! while posX posY posZ reflect the real positions in the FOV, localPosX localPosY localPosZ are relative to the volume BEFORE "any" rotation.
          // so not just the rsectorID rotation which we reconstruct later, but even the rotation that rotates the entire FOV before the simulation. It's quite confusing...
          // now we want to calculate the x y z in real FOV coordinates, moduleID, submoduleID, crystalID give us the x y coordinates in real FOV, because the the IDs rotate with the
          // elements they identify, of course. the localPosX etc coordinate do not, so the localPosX is always the doi coordinate because the crystals are placed in a system that
          // that in principle rotates around the z axis, then the entire FOV is rotatated around the Y axis to go back to the original ClearPEM coordinates and use
          // their reco algo (with that choiceof coordinates? who knows...). Crystal main axis is axis, as you can check in the GATE macro at the line
          // /gate/crystal/geometry/setXLength 15. mm
          // therefore for the z coordinate in FOV BEFORE the rsectorID rotation, the coordinate to use is the localPosX, hence averageDepEvents[iAvg].x, and NOT the localPosZ
          float doiAbs = rmin + (averageDepEvents[iAvg].x + (crylength/2.0));

          float zdoi;
          if(smearedDoi)
            zdoi = (float) gaussianSmear(doiAbs,doiResolutionFWHM);
          else
            zdoi = doiAbs;

          //force the event to be inside the crystal...
          if(fabs(zdoi) > fabs(rmin+crylength))
            if( zdoi > 0 )
              zdoi = rmin+crylength;
            else
              zdoi = - (rmin+crylength);
          if(fabs(zdoi) < fabs(rmin))
            if( zdoi > 0 )
              zdoi = rmin;
            else
              zdoi = - (rmin);

          // now we have the xyz coordinates if the sector was not rotated, but each sector is rotated and after all the FOV rotations our rot axis is x, so
          // the real coordinates will be
          float x = xabs;
          float y = yabs*cos(ang*averageDepEvents[iAvg].rsectorID) - zdoi*sin(ang*averageDepEvents[iAvg].rsectorID);
          float z = yabs*sin(ang*averageDepEvents[iAvg].rsectorID) + zdoi*cos(ang*averageDepEvents[iAvg].rsectorID);



          //debug
          // if(averageDepEvents[iAvg].trackID > 4 && averageDepEvents.size() > 1)
            // std::cout << i  << " " << averageDepEvents.size() << " " << averageDepEvents[iAvg].trackID << std::endl;
          point tempPoint;

          tempPoint.eventID = averageDepEvents[iAvg].eventID;
          tempPoint.parentID = averageDepEvents[iAvg].parentID;
          tempPoint.trackID = averageDepEvents[iAvg].trackID;
          tempPoint.sourceID = averageDepEvents[iAvg].sourceID;
          // calculate also the crystal and ring ID for STIR...
          // everything is recalculated according to rings
          Int_t N_MOD_xy  =  repmody;   // number of modules in xy plane
          Int_t N_SMOD_xy =  repsuby;   // number of submodules in xy plane
          Int_t N_CRY_xy  =  repcryy;   // number of crystals in xy plane
          Int_t N_SMOD_z  =  repsubz;   // number of submodules in z direction
          Int_t N_CRY_z   =  repcryz;   // number of crystals in z direction

          tempPoint.crystalIDforSTIR = averageDepEvents[iAvg].rsectorID * N_MOD_xy * N_SMOD_xy * N_CRY_xy
                                       + (averageDepEvents[iAvg].moduleID%N_MOD_xy) * N_SMOD_xy * N_CRY_xy
                                       + (averageDepEvents[iAvg].submoduleID%N_SMOD_xy) * N_CRY_xy
                                       + (averageDepEvents[iAvg].crystalID%N_CRY_xy);
          tempPoint.ringIDforSTIR    = (Int_t)(averageDepEvents[iAvg].crystalID/N_CRY_xy)
                                       + (Int_t)(averageDepEvents[iAvg].submoduleID/N_SMOD_xy)*N_CRY_z
                                      + (Int_t)(averageDepEvents[iAvg].moduleID/N_MOD_xy)*N_SMOD_z*N_CRY_z;
          tempPoint.x = x;
          tempPoint.y = y;
          tempPoint.z = z;
          tempPoint.time = averageDepEvents[iAvg].time;
          if(smearedEnergy)
          {
            tempPoint.energy = (Float_t) gaussianSmear(averageDepEvents[iAvg].energy,averageDepEvents[iAvg].energy*energyResolutionFWHM);
          }
          else
          {
            tempPoint.energy = averageDepEvents[iAvg].energy;
          }

          points.push_back(tempPoint);
        }

        // outAvg = averageDepEvents;
        outputTree->Fill();
        // if(averageDepEvents.size()<3)
        // {
        //   // std::cout << eventCounter << " " << averageDepEvents[iAvg].globalCryID << " " << averageDepEvents[iAvg].time << " " << averageDepEvents[iAvg].energy << " " << averageDepEvents[iAvg].x << " " << averageDepEvents[iAvg].y << " " << averageDepEvents[iAvg].z << " " << averageDepEvents[iAvg].parentID << " " << averageDepEvents[iAvg].trackID <<   std::endl;
        // }

          //
          // eventCounter++;
          // for(int iAvg = 0; iAvg < averageDepEvents.size() ; iAvg++)
          // {
          //   std::cout << eventCounter << " " << averageDepEvents[iAvg].globalCryID << " " << averageDepEvents[iAvg].time << " " << averageDepEvents[iAvg].energy << " " << averageDepEvents[iAvg].x << " " << averageDepEvents[iAvg].y << " " << averageDepEvents[iAvg].z << " " << averageDepEvents[iAvg].parentID << " " << averageDepEvents[iAvg].trackID <<   std::endl;
          // }
        // }

      }
      // std::cout<< "----------------------" << std::endl;


      //move to next pair of gammas, reset
      eventsCheck = HITSeventID;
      energyDeposition.clear();
      // points.clear();
    }
    //continue (or start) adding to the event

    //first thing first, go back to the original parent. It can be
    // std::cout << HITSeventID << " "<< HITStrackID << " " << HITSparentID << " " << HITSedep <<  std::endl;

    enDep hitDeposition;       //create a new endep struct
    hitDeposition.localPosX    = HITSlocalPosX;
    hitDeposition.localPosY    = HITSlocalPosY;
    hitDeposition.localPosZ    = HITSlocalPosZ;
    hitDeposition.eventID      = HITSeventID;
    hitDeposition.parentID     = HITSparentID;
    hitDeposition.trackID      = HITStrackID;
    hitDeposition.globalCryID  = HITSrsectorID * (nModulesXRSector * nSubmodulesXModule * nCrystalXSubmodule * nLayersXCrystal) +
                                 HITSmoduleID  * (nSubmodulesXModule * nCrystalXSubmodule * nLayersXCrystal) +
                                 HITSsubmoduleID * (nCrystalXSubmodule * nLayersXCrystal) +
                                 HITScrystalID * (nLayersXCrystal) +
                                 HITSlayerID;
    hitDeposition.gantryID     = HITSgantryID;
    hitDeposition.rsectorID    = HITSrsectorID;
    hitDeposition.moduleID     = HITSmoduleID;
    hitDeposition.submoduleID  = HITSsubmoduleID;
    hitDeposition.crystalID    = HITScrystalID;
    hitDeposition.layerID      = HITSlayerID;
    hitDeposition.sourceID     = HITSsourceID;
    hitDeposition.time         = HITStime;
    hitDeposition.edep         = HITSedep;
    energyDeposition.push_back(hitDeposition); //push it to the vector
    // if(HITStrackID > 2)

    statuscounter++;
    int perc = ((100*statuscounter)/HITSnentries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      // std::cout << statuscounter << std::endl;
    }

  }
  std::cout << std::endl;

  // for (Int_t i = 0 ; i < COINCnentries ; i++)
  // {
  //   Coincidences->GetEntry(i);
  //   //
  //   if (eventID1 == eventID2)
	//   {
	//     if ((comptonPhantom1 == 0) && (comptonPhantom2 == 0)) Trues++;
	//     else Scatters++;
	//   }
  //   else Randoms++;
  //
  //
  // }
  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
  // std::cout << "Trues    = " << Trues << std::endl;
  // std::cout << "Scatters = " << Scatters << std::endl;
  // std::cout << "Randoms  = " << Randoms << std::endl;

  return 0;
}
