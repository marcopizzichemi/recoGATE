// compile with
// g++ -o dataset dataset.cpp `root-config --cflags --libs`

// This program takes and output of gateToString, a root file
// then produces a list of points in variuos text files
// x1 y1 z1 x2 y2 z2
// with different strategies:
// 1: only events where energy deposition was in 2 and only 2 crystals
// 2: events where energy dep was in 3 and only 3 crystals. this is split in more subcases
//    2a. the scatter events are averaged in one average event
//    2b. the first hit of the scatter is saved (100% efficient compton reconstruction)
//    2c. the first hit of the scatter is recognized with a given efficiency (parameter compton_eff below)
// data is saved in 4 different files


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

// DLS Parameters -- CONVERSION TO STIR PART
#define  N_RINGS              128 // DLS: Number of rings - these are the number of crystals in z direction (?)
#define  N_DET                384 // DLS: Detectors per ring - these the the number of crystals in a ring of 1 crystals height (?)
#define  S_WIDTH              128 // DLS: Number of radial sinogram bins - binning of radial sinogram. i guess the same of N_RINGS is fine?
#define  MAX_D_RING           128 // DLS: Maximum ring difference - in principle all rings accepted
#define  USE_OFFSET             0 // DLS: On/Off use of offset
#define  OFFSET               499 // DLS: Sets initial sinogram angle

struct Pairs // definition fo a pair of points
{
  float x1;
  float y1;
  float z1;
  float x2;
  float y2;
  float z2;
} __attribute__((__packed__));


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

struct EventFormat {
		double ts;
		u_int8_t random;
		float d;
		float yozRot;
		float x1;
		float y1;
		float z1;
		float e1;
	  u_int8_t n1;
		float x2;
		float y2;
		float z2;
		float e2;
	  u_int8_t n2;
		float dt;
	} __attribute__((__packed__));


//STIR PART
Float_t    Mich_r1r2fu[N_RINGS][N_RINGS][N_DET/2][S_WIDTH]={0};

int main(int argc, char** argv)
{
  if(argc<6) {
    std::cout<<"You need to provide an input file "<<std::endl ;
    std::cout<<"USAGE:"<<std::endl ;
    std::cout<<"dataset input.root output.root [binary <0|1>] [listmodeoutput <0|1>] [NoBackground <0|1>] [stir <0|1>] "<<std::endl ;
    return 1;
  }

  unsigned short ans;

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

  //variables with defaults
  bool binary = false;
  bool listmodeoutput = true;
  bool NoBackground = false;
  bool stir = false;
  float compton_eff = 0.65; // fraction efficiency of compton assignment. e.g. 0.7 means 70%
  float enMin = 100.0;
  float enMax = 700.0;
  float plates_distance = 200.0;
  float rotationAngle = 0.0;
  //INPUT FILE AND TTREE
  std::string inputfilename,outputfilename ;
  // inputfilename = argv[1] ;
  // outputfilename = argv[2];

  static struct option longOptions[] =
  {
      { "inner-radius",required_argument, 0, 0},
      { "compton-efficiency",required_argument, 0, 0},
      { "energy-low",required_argument, 0, 0},
      { "energy-high",required_argument, 0, 0},
      { "binary-output",no_argument, 0, 0},
      { "listmode-output",no_argument, 0, 0},
      { "stir-output",no_argument, 0, 0},
      { "no-background",no_argument, 0, 0},
      { "rotation-angle",required_argument, 0, 0},
			{ NULL, 0, 0, 0 }
	};

  while(1) {

		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:", longOptions, &optionIndex);

		if (c == -1) {
			break;
		}

		if (c == 'i'){
			inputfilename = (char *)optarg;
    }
		else if (c == 'o'){
      outputfilename = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      plates_distance = 2.0*atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 1){
      compton_eff = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 2){
      enMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 3){
      enMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      binary   = true;
      listmodeoutput = false;
      stir     = false;
    }
    else if (c == 0 && optionIndex == 5){
      binary   = false;
      listmodeoutput = true;
      stir     = false;
    }
    else if (c == 0 && optionIndex == 6){
      binary   = false;
      listmodeoutput = false;
      stir     = true;
    }
    else if (c == 0 && optionIndex == 7){
      NoBackground = true;
    }
    else if (c == 0 && optionIndex == 8){
      rotationAngle = atof((char *)optarg);
    }
		else { //FIXME
			std::cout	<< "Usage: " << argv[0] << std::endl
			<< "\t\t" << "[ -i <input file> ] " << std::endl
			<< "\t\t" << "[ -o <output file> ] " << std::endl
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


  std::cout << "Rotation angle is = " << rotationAngle << std::endl;
  std::cout << "Input file name is " << inputfilename << std::endl;
  TFile *inputFile = TFile::Open(inputfilename.c_str());
  TTree *tree = (TTree*) inputFile->Get("test");
  std::vector<point> *points = 0;
  tree->SetBranchAddress("points",&points);

  //OUTPUT FILES for STIR
  std::string Moutputfilename, Poutputfilename ;
  FILE  *Mich_r1r2fuFile, *Proj_File;


  std::string::size_type t_point = outputfilename.find_last_of(".");
  // path = path.substr(0,t);

  // std::stringstream ssBinary(argv[3]);
  //
  // std::stringstream sslistmodeoutput(argv[4]);
  // std::stringstream ssNoBackground(argv[5]);
  // std::stringstream ssstir(argv[6]);
  //
  //
  //
  // ssBinary >>  binary;
  // sslistmodeoutput >> listmodeoutput;
  // ssNoBackground >> NoBackground;
  // ssstir >> stir;


  std::cout << binary << listmodeoutput << NoBackground << stir << std::endl;

  std::string baseName = outputfilename.substr(0,t_point);
  if(NoBackground)
  {
    baseName += "_NOBG_";
  }
  //OUTPUT STREAMS
  std::string ofs2cryName                =  baseName + "_2cry";
  std::string ofs3cry_avgName            =  baseName + "_3cry-avg";
  std::string ofs3cry_magicalComptonName =  baseName + "_3cry-magicalCompton";
  std::string ofs3cry_effComptonName     =  baseName + "_3cry-effCompton";

  std::cout << ofs2cryName << " " << ofs3cry_avgName <<std::endl;

  Moutputfilename = "./Mich_"+ baseName + ".s" ;
  std::cout << "Michelogram file name is = " << Moutputfilename << std::endl ;
  Poutputfilename = "./Proj_"+ baseName + ".s" ;
  std::cout << "Projection file name is = " << Poutputfilename << std::endl ;


  // Pairs pair;

  std::ofstream ofs2cry;
  std::ofstream ofs3cry_avg;
  std::ofstream ofs3cry_magicalCompton;
  std::ofstream ofs3cry_effCompton;

  // bool listmodeoutput = true;
  // std::ofstream *listMode2cry = NULL;
  // std::ofstream *listMode3cry_avg = NULL;
  // std::ofstream *listMode3cry_magicalCompton = NULL;
  // std::ofstream *listMode3cry_effCompton = NULL;



  // CurrentOutput = new std::ofstream(OutStringStream.str().c_str(), std::ios::binary);

  if(binary)
  {
    ofs2cryName+= ".bin";
    ofs3cry_avgName+= ".bin";
    ofs3cry_magicalComptonName+= ".bin";
    ofs3cry_effComptonName+= ".bin";

    ofs2cry.open (ofs2cryName.c_str(), std::ios::binary);
    ofs3cry_avg.open (ofs3cry_avgName.c_str(), std::ios::binary);
    ofs3cry_magicalCompton.open (ofs3cry_magicalComptonName.c_str(), std::ios::binary);
    ofs3cry_effCompton.open (ofs3cry_effComptonName.c_str(), std::ios::binary);
  }
  else
  {
    if(listmodeoutput)
    {
      ofs2cryName+= ".elm2";
      ofs3cry_avgName+= ".elm2";
      ofs3cry_magicalComptonName+= ".elm2";
      ofs3cry_effComptonName+= ".elm2";

      ofs2cry.open (ofs2cryName.c_str(), std::ios::binary);
      ofs3cry_avg.open (ofs3cry_avgName.c_str(), std::ios::binary);
      ofs3cry_magicalCompton.open (ofs3cry_magicalComptonName.c_str(), std::ios::binary);
      ofs3cry_effCompton.open (ofs3cry_effComptonName.c_str(), std::ios::binary);
    }
    else
    {
      ofs2cryName+= ".txt";
      ofs3cry_avgName+= ".txt";
      ofs3cry_magicalComptonName+= ".txt";
      ofs3cry_effComptonName+= ".txt";

      ofs2cry.open (ofs2cryName.c_str(), std::ofstream::out);
      ofs3cry_avg.open (ofs3cry_avgName.c_str(), std::ofstream::out);
      ofs3cry_magicalCompton.open (ofs3cry_magicalComptonName.c_str(), std::ofstream::out);
      ofs3cry_effCompton.open (ofs3cry_effComptonName.c_str(), std::ofstream::out);
    }


  }
  if(stir)
  {
    Mich_r1r2fuFile = fopen(Moutputfilename.c_str(),"wb");
    Proj_File       = fopen(Poutputfilename.c_str(),"wb");
  }


  //COUNTERS
  Int_t nsamples = tree->GetEntries();
  std::cout << nsamples << std::endl;
  Int_t onlyOneCrystals = 0;
  Int_t onlyTwoCrystals = 0;
  Int_t onlyThreeCrystals = 0;
  Int_t onlyFourCrystals = 0;
  Int_t moreThanFourCrystals = 0;
  Int_t onlyTwoCrystalsInEnergyWindow = 0;
  Int_t onlyThreeCrystalsInEnergyWindow = 0;
  Int_t onlyFourCrystalsInEnergyWindow = 0;
  Int_t moreThanFourCrystalsInEnergyWindow = 0;

  //COMPTON efficiency parameter

  double timeCounter = 0.0; //fake time for the elm2 format

  for(int i = 0 ; i < nsamples ; i ++)
  {

    tree->GetEntry(i);
    timeCounter += 1e-6;
    Pairs pair;  //the output pair
    //--------------------------//
    // 1 crystal hit            //
    //--------------------------//
    // if(NoBackground)
    // {
    //   if(points->at(0).sourceID == 0)
    //   {
    //     continue;
    //   }
    // }
    if(points->size() > 0) //why points == 0????
    {
      // std::cout << points->at(0).sourceID << std::endl;
      Int_t EventSourceID = points->at(0).sourceID;

      if(NoBackground && EventSourceID == 0)
      {

      }
      else
      {
        if(points->size() < 2)  //single crystal hit, count and discard
        {
          onlyOneCrystals++;
        }
        //--------------------------//


        //selection of dataset. First, only when 2 and only 2 crystals are hit, and with the energy in the correct window -> more or less standard clearpem strategy
        //--------------------------//
        // 2 crystals hit           //
        //--------------------------//
        if(points->size() == 2)
        {
          if(points->at(0).z * points->at(1).z > 0)
          {
            //two hits in same head!
          }
          else
          {
            onlyTwoCrystals++;
            if(1000.0*points->at(0).energy > enMin && 1000.0*points->at(0).energy < enMax && 1000.0*points->at(1).energy > enMin && 1000.0*points->at(1).energy < enMax) //energy limits
            {
              onlyTwoCrystalsInEnergyWindow++; //baseline sensitivity
            }
            if(listmodeoutput)
            {
              EventFormat fe;
              fe.ts     = timeCounter;
              fe.random = 0;
              fe.d      = plates_distance;
              fe.yozRot = (rotationAngle / 180.0) * 3.1415;
              fe.x1     = points->at(0).x;
              fe.y1     = points->at(0).y;
              fe.z1     = points->at(0).z;
              fe.e1     = points->at(0).energy * 1000;  // clearpem reco wants KeV
              fe.n1     = 1;
              fe.x2     = points->at(1).x;
              fe.y2     = points->at(1).y;
              fe.z2     = points->at(1).z;
              fe.e2     = points->at(1).energy * 1000;// clearpem reco wants KeV
              fe.n2     = 1;
              // fe.dt     = (double) (points->at(0).time - points->at(1).time);
              fe.dt     = (double) 1e-9;
              ofs2cry.write((char*)&fe,sizeof(fe));
            }
          }
          // if(points->at(0).energy > enMin && points->at(0).energy < enMax && points->at(1).energy > enMin && points->at(1).energy < enMax) //energy limits
          // {
          //   onlyTwoCrystalsInEnergyWindow++; //baseline sensitivity
          //   if(binary)
          //   {
          //     pair.x1 = points->at(0).x;
          //     pair.y1 = points->at(0).y;
          //     pair.z1 = points->at(0).z;
          //     pair.x2 = points->at(1).x;
          //     pair.y2 = points->at(1).y;
          //     pair.z2 = points->at(1).z;
          //     ofs2cry.write((char*)&pair,sizeof(pair));
          //   }
          //   else
          //   {
          //     if(listmodeoutput)
          //     {
          //       EventFormat fe;
          //       fe.ts     = 0;
          //       fe.random = 0;
          //       fe.d      = 200.0; // for now hardcoded
          //       fe.yozRot = 0;
          //       //change axis to accomodate clearpem directions...
          //       fe.x1     = points->at(0).z;
          //       fe.y1     = points->at(0).y;
          //       fe.z1     = -points->at(0).x;
          //       fe.e1     = points->at(0).energy;
          //       fe.n1     = 1;
          //       fe.x2     = points->at(1).z;
          //       fe.y2     = points->at(1).y;
          //       fe.z2     = -points->at(1).x;
          //       fe.e2     = points->at(1).energy;
          //       fe.n2     = 1;
          //       fe.dt     = points->at(0).time - points->at(1).time;
          //       ofs2cry.write((char*)&fe,sizeof(fe));
          //     }
          //     else
          //     {
          //       ofs2cry   << points->at(0).x << " "
          //       << points->at(0).y << " "
          //       << points->at(0).z << " "
          //       << points->at(1).x << " "
          //       << points->at(1).y << " "
          //       << points->at(1).z
          //       << std::endl;
          //     }
          //   }
          //
          //   if(stir)
          //   {
          //     //rings
          //     Int_t ring1 = points->at(0).ringIDforSTIR;
          //     Int_t ring2 = points->at(1).ringIDforSTIR;
          //     Int_t crystal1 = points->at(0).crystalIDforSTIR;
          //     Int_t crystal2 = points->at(1).crystalIDforSTIR;
          //
          //     //
          //     //--------------------------------
          //     //  Bin the crystal ring pairs into Michelograms
          //     //  u - radial sinogram component
          //     //  phi - azimuthal sinogram component
          //     //  ring pairs are sorted according to c1 < c2 else flip
          //     //  where c1 and c2 are crystals at phi(u = S_WIDTH/2)
          //     //--------------------------------
          //     Int_t    phi, u;
          //     int flip, swap, zi, c1, c2;
          //     phi = ((crystal1 + crystal2 + N_DET/2)%N_DET)/2;
          //
          //     if (((crystal1 + crystal2) < (3*N_DET/2)) && ((crystal1 + crystal2) >= (N_DET/2)))
          //     u    =  abs(crystal1 - crystal2) -  N_DET/2 + S_WIDTH/2;
          //     else u = -abs(crystal1 - crystal2) +  N_DET/2 + S_WIDTH/2;
          //
          //     if ( u >= S_WIDTH || u < 0 ) continue;
          //
          //     if (u%2 == 0)
          //     {
          //       zi = (N_DET/2 - (crystal1 - crystal2) - 1)/2;
          //       if (zi >=  N_DET/4) zi = zi - N_DET/2 + 1;
          //       if (zi <= -N_DET/4) zi = zi + N_DET/2 - 1;
          //     }
          //     else
          //     {
          //       zi = (N_DET/2 - (crystal1 - crystal2))/2;
          //       if (zi >=  N_DET/4) zi = zi - N_DET/2;
          //       if (zi <= -N_DET/4) zi = zi + N_DET/2;
          //     }
          //
          //     c1 = crystal1 + zi;
          //     c2 = crystal2 - zi;
          //     if (c1 >= N_DET) c1 = c1 - N_DET;
          //     if (c1 < 0)      c1 = c1 + N_DET;
          //     if (c2 >= N_DET) c2 = c2 - N_DET;
          //     if (c2 < 0)      c2 = c2 + N_DET;
          //
          //     if (c1 < c2) flip = 0;
          //     else         flip = 1;
          //
          //     if (flip)
          //     {
          //       swap  = ring1;
          //       ring1 = ring2;
          //       ring2 = swap;
          //     }
          //
          //     // Update the different arrays...
          //     //-------------------------------
          //     //***ALL EVENTS
          //     Mich_r1r2fu[ring2][ring1][phi][u] += 1.;
          //
          //   }
          //   //DEBUG
          //   // std::cout   << points->at(0).x << " "
          //   //           << points->at(0).y << " "
          //   //           << points->at(0).z << " "
          //   //           << points->at(1).x << " "
          //   //           << points->at(1).y << " "
          //   //           << points->at(1).z
          //   //           << std::endl;
          //
          // }

          //DEBUG
          // std::cout << points->at(0).x << " "
          //           << points->at(0).y << " "
          //           << points->at(0).z << " "
          //           << points->at(1).x << " "
          //           << points->at(1).y << " "
          //           << points->at(1).z
          //           << std::endl;
        }
        //--------------------------//


        //--------------------------//
        // 3 crystals hit           //
        //--------------------------//
        if(points->size() == 3)
        {
          //indentification of the pair. at least one must be in the "511 window", i.e. should have experience photoelectric deposition by one of the 2 gammas
          bool foundSingle511 = false;
          int SingleID = -1;
          for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          {
            if(1000.0*points->at(pCount).energy > enMin && 1000.0*points->at(pCount).energy < enMax)
            {
              foundSingle511 = true;
              SingleID = pCount;    //take the ID of the single 511 deposition
            }
          }
          if(foundSingle511) //now take the other two
          {
            float sum_energy = 0;
            for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
            {
              if(pCount != SingleID) //only the others
              {
                //sum the energies
                sum_energy += points->at(pCount).energy;
              }
            }
            //check if sum is in window
            if(1000.0*sum_energy > enMin && 1000.0*sum_energy < enMax)
            {
              //count this event as 3 crystals only, in energy window
              onlyThreeCrystalsInEnergyWindow++;
              //---------------------------------------------//
              //strategy a: calculate an average event
              //---------------------------------------------//
              float temp_x = 0;
              float temp_y = 0;
              float temp_z = 0;
              for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
              {
                if(pCount != SingleID)
                {
                  //average x,y,z
                  temp_x += (points->at(pCount).x * points->at(pCount).energy) / sum_energy;
                  temp_y += (points->at(pCount).y * points->at(pCount).energy) / sum_energy;
                  temp_z += (points->at(pCount).z * points->at(pCount).energy) / sum_energy;
                }
              }
              //save in another file
              if(binary)
              {
                pair.x1 = points->at(SingleID).x;
                pair.y1 = points->at(SingleID).y;
                pair.z1 = points->at(SingleID).z;
                pair.x2 = temp_x;
                pair.y2 = temp_y;
                pair.z2 = temp_z;
                ofs3cry_avg.write((char*)&pair,sizeof(pair));
              }
              else
              {
                ofs3cry_avg << points->at(SingleID).x << " "
                << points->at(SingleID).y << " "
                << points->at(SingleID).z << " "
                << temp_x << " "
                << temp_y << " "
                << temp_z
                << std::endl;
              }


              //---------------------------------------------//
              //strategy b: magical compton
              //---------------------------------------------//
              //check which crystal was hit first
              Int_t firstCrystalID;
              float timeOfFirst = INFINITY;
              for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
              {
                if(pCount != SingleID)
                {
                  if(points->at(pCount).time <  timeOfFirst)
                  {
                    timeOfFirst = points->at(pCount).time;
                    firstCrystalID = pCount;
                  }
                }
              }

              //now save only the single 511 and the position of first between the two comptons (100% efficient compton discrimination)
              if(binary)
              {
                pair.x1 = points->at(SingleID).x;
                pair.y1 = points->at(SingleID).y;
                pair.z1 = points->at(SingleID).z;
                pair.x2 = points->at(firstCrystalID).x;
                pair.y2 = points->at(firstCrystalID).y;
                pair.z2 = points->at(firstCrystalID).z;
                ofs3cry_magicalCompton.write((char*)&pair,sizeof(pair));
              }
              else
              {
                ofs3cry_magicalCompton << points->at(SingleID).x << " "
                << points->at(SingleID).y << " "
                << points->at(SingleID).z << " "
                << points->at(firstCrystalID).x << " "
                << points->at(firstCrystalID).y << " "
                << points->at(firstCrystalID).z
                << std::endl;
              }


              //---------------------------------------------//
              //strategy c: compton with efficiency
              //---------------------------------------------//
              //get the remaing id among 0,1,2
              Int_t otherID = -1;
              for(int fastI = 0 ; fastI < 3 ; fastI++)
              {
                if(fastI != SingleID && fastI != firstCrystalID)
                otherID = fastI;
              }
              if(otherID == -1) std::cout << "ERROR!!!!!!!!!!!!" << std::endl;
              // choose the right one with an efficiency
              double dice = (double)rand() / (double)RAND_MAX;
              Int_t chosenCryID;
              if(dice <= compton_eff)
              chosenCryID = firstCrystalID;
              else
              chosenCryID = otherID;

              if(binary)
              {
                pair.x1 = points->at(SingleID).x;
                pair.y1 = points->at(SingleID).y;
                pair.z1 = points->at(SingleID).z;
                pair.x2 = points->at(chosenCryID).x;
                pair.y2 = points->at(chosenCryID).y;
                pair.z2 = points->at(chosenCryID).z;
                ofs3cry_effCompton.write((char*)&pair,sizeof(pair));
              }
              else
              {
                ofs3cry_effCompton << points->at(SingleID).x    << " "
                << points->at(SingleID).y    << " "
                << points->at(SingleID).z    << " "
                << points->at(chosenCryID).x << " "
                << points->at(chosenCryID).y << " "
                << points->at(chosenCryID).z
                << std::endl;
              }

            }
            else // the two scatter don't sum up to enter 511 window
            {
              //don't save the event for "reconstruction"
            }
          }
          else // no crystal among the 3 has a deposition that enters the 511 window
          {
            //don't save the event for "reconstruction"
          }
          onlyThreeCrystals++;
        }
        //--------------------------//


        //--------------------------//
        // 4 crystals hit           //
        //--------------------------//
        if(points->size() == 4)
        {
          onlyFourCrystals++;
          float sum_energy = 0;
          for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          {
            sum_energy += points->at(pCount).energy;
          }
          if(1000.0*sum_energy > enMin && 1000.0*sum_energy < enMax)
          {
            onlyFourCrystalsInEnergyWindow++;
          }

        }
        //--------------------------//


        //--------------------------//
        // more than 4 crystals hit //
        //--------------------------//
        if(points->size() > 4)
        {
          moreThanFourCrystals++;
          float sum_energy = 0;
          for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          {
            sum_energy += points->at(pCount).energy;
          }
          if(1000.0*sum_energy > enMin && 1000.0*sum_energy < enMax)
          {
            moreThanFourCrystalsInEnergyWindow++;
          }
        }
        //--------------------------//


        //FEEDBACK TO USER
        int perc = ((100*i)/nsamples); //should strictly have not decimal part, written like this...
        if( (perc % 10) == 0 )
        {
          std::cout << "\r";
          std::cout << perc << "% done... ";
          // std::cout << statuscounter << std::endl;
        }
      }
    }
  }
  std::cout << std::endl;

  //SENSITIVITY COUNTS
  std::cout << "Total dataset                       = " << nsamples << std::endl;
  std::cout << "Only 1 crystal                      = " << onlyOneCrystals << std::endl;
  std::cout << "Only 2 crystals                     = " << onlyTwoCrystals << "\t , @511 = " << onlyTwoCrystalsInEnergyWindow<< std::endl;
  std::cout << "Only 3 crystals                     = " << onlyThreeCrystals << "\t , @511 = "<< onlyThreeCrystalsInEnergyWindow << std::endl;
  std::cout << "Only 4 crystals                     = " << onlyFourCrystals << "\t , @511 = "<< onlyFourCrystalsInEnergyWindow <<std::endl;
  std::cout << "More than 4 crystals                = " << moreThanFourCrystals << "\t , @511 = "<< moreThanFourCrystalsInEnergyWindow <<std::endl;

  // Write the data to disk, and then close Michelogram file...
  //----------------------------------------------------------------
  if(stir)
  {
    ans = fwrite(Mich_r1r2fu,4,(N_RINGS*N_RINGS*N_DET/2*S_WIDTH),Mich_r1r2fuFile);
    fclose(Mich_r1r2fuFile);
  }
  // Generate projection files...
  //-----------------------------
  // From Segment number -MAX_D_RING to +MAX_D_RING
  // After phi
  // After z
  // After r
  Int_t    ring1, ring2;
  Int_t S_NUM;
  if(stir)
  {
    for (int i = 0 ; i < 2*MAX_D_RING + 1 ; i++)
    {
      if (i <= MAX_D_RING) S_NUM = N_RINGS - MAX_D_RING + i;
      else                 S_NUM = N_RINGS + MAX_D_RING - i;

      for (int j = 0 ; j < N_DET/2 ; j++)
      {
        float Proj[S_NUM][S_WIDTH];
        for (int k = 0 ; k < S_NUM ; k++)
        {
          if (i <= MAX_D_RING) ring1 = k;
          else                 ring2 = k;

          if (i <= MAX_D_RING) ring2 = ring1 + MAX_D_RING - i;
          else                 ring1 = ring2 - MAX_D_RING + i;

          for (int l = 0 ; l < S_WIDTH ; l++)
          Proj[k][l] = Mich_r1r2fu[ring2][ring1][j][l];
        }
        ans = fwrite(Proj,4,(S_NUM*S_WIDTH),Proj_File);
      }
    }
    fclose(Proj_File);
  }



  //CLOSE FILES
  inputFile->Close();
  ofs2cry.close();
  ofs3cry_avg.close();
  ofs3cry_magicalCompton.close();
  ofs3cry_effCompton.close();
  return 0;
}
