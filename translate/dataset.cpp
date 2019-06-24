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
//    2d. the first hit of the scatter is assigned to the crystal of maximum energy deposition
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
  Int_t primaryID;
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


struct primary_t
{
  int id;
  int num;
};


  //distance between 2 points in 3d space
Float_t distance3D(Float_t ax, Float_t ay, Float_t az, Float_t bx, Float_t by, Float_t bz)
{
 Float_t v[3] = {bx-ax, by-ay, bz-az};
 Float_t vMod = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
 return fabs(vMod);
}

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
  float energyFwhm = 0.12;
  //INPUT FILE AND TTREE
  std::string inputfilename,outputfilename ;
  // inputfilename = argv[1] ;
  // outputfilename = argv[2];
  int effSteps = 11;
  float min_efficiency = 0.5;
  float max_efficiency = 1.0;
  float gammaEnergy = 511.0;

  static struct option longOptions[] =
  {
      { "inner-radius",required_argument, 0, 0},
      { "eff-steps",required_argument, 0, 0},
      { "energy-fwhm",required_argument, 0, 0},
      // { "energy-high",required_argument, 0, 0},
      { "binary-output",no_argument, 0, 0},
      { "listmode-output",no_argument, 0, 0},
      { "stir-output",no_argument, 0, 0},
      { "no-background",no_argument, 0, 0},
      { "rotation-angle",required_argument, 0, 0},
      { "gamma-energy",required_argument, 0, 0},
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
      effSteps = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 2){
      energyFwhm = atof((char *)optarg);
    }
    // else if (c == 0 && optionIndex == 3){
    //   enMax = atof((char *)optarg);
    // }
    else if (c == 0 && optionIndex == 3){
      binary   = true;
      listmodeoutput = false;
      stir     = false;
    }
    else if (c == 0 && optionIndex == 4){
      binary   = false;
      listmodeoutput = true;
      stir     = false;
    }
    else if (c == 0 && optionIndex == 5){
      binary   = false;
      listmodeoutput = false;
      stir     = true;
    }
    else if (c == 0 && optionIndex == 6){
      NoBackground = true;
    }
    else if (c == 0 && optionIndex == 7){
      rotationAngle = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      gammaEnergy = atof((char *)optarg);
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
			<< "\t\t" << "[ --no-interfile ]" << std::endl
			<< "\t\t" << "[ --partial <event fraction %>" << std::endl
			<< "\t\t" << "[ --rand-image <random file> ]" << std::endl
			<< "\t\t" << "[ --rand-cut <random cut value, 0 to 0.5> ]" << std::endl
			<< "\t\t" << "[ --save-iterations <N> ] " << std::endl
			<< "\t\t" << std::endl;

			return 1;
		}
	}



  enMin =gammaEnergy - (gammaEnergy * energyFwhm);
  enMax = gammaEnergy + 2.0*(gammaEnergy * energyFwhm);
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
  // std::string ofs3cry_magicalComptonName =  baseName + "_3cry-magicalCompton";
  // std::string ofs3cry_effComptonName     =  baseName + "_3cry-effCompton";
  std::string ofs3cry_maxEnergyName      =  baseName + "_3cry-maxEnergy";
  std::string summary_out_name           = baseName + "_summary";
  std::vector<std::string>  ofs3cry_effComptonName;
  float efficiency_step = (float) (max_efficiency - min_efficiency)/(effSteps-1);
  for(int i = 0 ; i < effSteps ; i++)
  {
    float efficiency = min_efficiency + i*efficiency_step;
    std::stringstream sname;
    sname << baseName << "_3cry-eff" << efficiency;
    ofs3cry_effComptonName.push_back(sname.str());
  }

  std::cout << ofs2cryName << " " << ofs3cry_avgName <<std::endl;

  Moutputfilename = "./Mich_"+ baseName + ".s" ;
  std::cout << "Michelogram file name is = " << Moutputfilename << std::endl ;
  Poutputfilename = "./Proj_"+ baseName + ".s" ;
  std::cout << "Projection file name is = " << Poutputfilename << std::endl ;


  // Pairs pair;

  std::ofstream ofs2cry;
  std::ofstream ofs3cry_avg;
  // std::ofstream ofs3cry_magicalCompton;
  // std::ofstream ofs3cry_effCompton;
  std::ofstream ofs3cry_maxEnergy;
  std::vector<std::ofstream*> ofs3cry_effCompton;
  std::ofstream summary_out;

  summary_out_name += ".txt";
  summary_out.open (summary_out_name.c_str(),  std::ofstream::out);
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
    // ofs3cry_magicalComptonName+= ".bin";
    // ofs3cry_effComptonName+= ".bin";
    ofs3cry_maxEnergyName+= ".bin";

    ofs2cry.open (ofs2cryName.c_str(), std::ios::binary);
    ofs3cry_avg.open (ofs3cry_avgName.c_str(), std::ios::binary);
    // ofs3cry_magicalCompton.open (ofs3cry_magicalComptonName.c_str(), std::ios::binary);
    // ofs3cry_effCompton.open (ofs3cry_effComptonName.c_str(), std::ios::binary);
    ofs3cry_maxEnergy.open (ofs3cry_maxEnergyName.c_str(), std::ios::binary);
    for(int i = 0 ; i < effSteps ; i++)
    {
      std::stringstream sname;
      sname << ofs3cry_effComptonName[i] << ".bin";
      std::ofstream* temp_ofs;
      temp_ofs = new std::ofstream(sname.str().c_str(), std::ios::binary);
      // temp_ofs->open (sname.str(), std::ios::binary);
      ofs3cry_effCompton.push_back(temp_ofs);
    }


  }
  else
  {
    if(listmodeoutput)
    {
      ofs2cryName+= ".elm2";
      ofs3cry_avgName+= ".elm2";
      ofs3cry_maxEnergyName+= ".elm2";

      ofs2cry.open (ofs2cryName.c_str(), std::ios::binary);
      ofs3cry_avg.open (ofs3cry_avgName.c_str(), std::ios::binary);
      ofs3cry_maxEnergy.open (ofs3cry_maxEnergyName.c_str(), std::ios::binary);

      for(int i = 0 ; i < effSteps ; i++)
      {
        std::stringstream sname;
        sname << ofs3cry_effComptonName[i] << ".elm2";
        // std::cout << sname.str() << std::endl;
        std::ofstream* temp_ofs;
        temp_ofs = new std::ofstream(sname.str().c_str(), std::ios::binary);
        // temp_ofs->open (sname.str().c_str(), std::ios::binary);
        // std::cout << sname.str() << std::endl;
        ofs3cry_effCompton.push_back(temp_ofs);
        // std::cout << sname.str() << std::endl;
      }
    }
    else
    {
      ofs2cryName+= ".txt";
      ofs3cry_avgName+= ".txt";
      // ofs3cry_magicalComptonName+= ".txt";
      // ofs3cry_effComptonName+= ".txt";
      ofs3cry_maxEnergyName+= ".txt";

      ofs2cry.open (ofs2cryName.c_str(), std::ofstream::out);
      ofs3cry_avg.open (ofs3cry_avgName.c_str(), std::ofstream::out);
      // ofs3cry_magicalCompton.open (ofs3cry_magicalComptonName.c_str(), std::ofstream::out);
      // ofs3cry_effCompton.open (ofs3cry_effComptonName.c_str(), std::ofstream::out);
      ofs3cry_maxEnergy.open (ofs3cry_maxEnergyName.c_str(),  std::ofstream::out);
      for(int i = 0 ; i < effSteps ; i++)
      {
        std::stringstream sname;
        sname << ofs3cry_effComptonName[i] << ".txt";
        std::ofstream* temp_ofs;
        temp_ofs = new std::ofstream(sname.str().c_str(), std::ofstream::out);
        // temp_ofs->open (sname.str().c_str(), std::ofstream::out);
        ofs3cry_effCompton.push_back(temp_ofs);
      }
    }


  }

  std::cout << "qui" << std::endl;

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
  Int_t NoCrystals = 0;

  Int_t moreThanFourCrystals = 0;
  Int_t onlyOneCrystalsInEnergyWindow = 0;
  Int_t onlyTwoCrystalsInEnergyWindow = 0;
  Int_t onlyThreeCrystalsInEnergyWindow = 0;
  Int_t onlyFourCrystalsInEnergyWindow = 0;
  Int_t moreThanFourCrystalsInEnergyWindow = 0;
  Int_t onlyTwoCrystalsHitBySameGammaInEnergyWindow = 0;

  Int_t twoCrystalsSamePrimary = 0;
  Int_t twoCrystalsDiffPrimary = 0;
  Int_t moreThanThreeCrystals = 0;
  Int_t moreThanThreeCrystalsInEnergyWindow = 0;

  //COMPTON efficiency parameter

  double timeCounter = 0.0; //fake time for the elm2 format

  for(int i = 0 ; i < nsamples ; i ++)
  {

    tree->GetEntry(i);
    timeCounter += 1e-6;
    Pairs pair;  //the output pair

    if(points->size() == 0)
    {
      NoCrystals++;
    }

    if(points->size() > 0) //why points == 0????
    {
      // std::cout << points->at(0).sourceID << std::endl;
      Int_t EventSourceID = points->at(0).sourceID;

      if(NoBackground && EventSourceID == 0)
      {

      }
      else
      {
        if(points->size() == 1)  //single crystal hit, count and discard
        {
          onlyOneCrystals++;
          if(1000.0*points->at(0).energy > enMin &&
             1000.0*points->at(0).energy < enMax ) //energy limits
          {
            onlyOneCrystalsInEnergyWindow++;
          }
        }
        //--------------------------//

        //selection of dataset. First, only when 2 and only 2 crystals are hit, and with the energy in the correct window -> more or less standard clearpem strategy
        //--------------------------//
        // 2 crystals hit           //
        //--------------------------//
        if(points->size() == 2)
        {
          onlyTwoCrystals++;
          if(points->at(0).primaryID == points->at(1).primaryID)
          {
            twoCrystalsSamePrimary++;
            // onlyTwoCrystalsHitBySameGamma++;
            float sum_energy_for_compton_camera = 1000.0*points->at(0).energy + 1000.0*points->at(1).energy;
            if(sum_energy_for_compton_camera > enMin &&
               sum_energy_for_compton_camera < enMax) //energy limits
            {
              onlyTwoCrystalsHitBySameGammaInEnergyWindow++; //baseline sensitivity
            }
          }
          else
          {
            twoCrystalsDiffPrimary++;
            // classic PET event, two crystals hit, 511 KeV in each crystal
            if(1000.0*points->at(0).energy > enMin &&
               1000.0*points->at(0).energy < enMax &&
               1000.0*points->at(1).energy > enMin &&
               1000.0*points->at(1).energy < enMax) //energy limits
            {
              onlyTwoCrystalsInEnergyWindow++; //baseline sensitivity
            }
          }






          if(listmodeoutput) //save anyway all the coincidences, even with wrong energy, they are filtered by reco
          {
            EventFormat fe;
            fe.ts     = timeCounter;
            fe.random = 0;
            fe.d      = plates_distance;
            fe.yozRot = 0;
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
            fe.dt     = (double) (points->at(0).time - points->at(1).time);
            ofs2cry.write((char*)&fe,sizeof(fe));
          }
        }
        //--------------------------//

        //--------------------------//
        // 3 crystals hit           //
        //--------------------------//
        if(points->size() == 3)
        {
          onlyThreeCrystals++;

          //using primaryID
          // find single event with primaryID


          // std::vector<primary_t> primaryCounter;
          // //look for different primaries
          // for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          // {
          //   if(pCount == 0)
          //   {
          //     primary_t tempPrimary;
          //     tempPrimary.id = points->at(pCount).primaryID;
          //     tempPrimary.num = 1;
          //     primaryCounter.push_back(tempPrimary);
          //   }
          //   else
          //   {
          //     for(Int_t prI = 0 ; prI < primaryCounter.size() ; prI++)
          //     {
          //       if(primaryCounter[prI].id == points->at(pCount).primaryID)
          //       {
          //         primaryCounter[prI].num++;
          //       }
          //       else
          //       {
          //         primary_t tempPrimary;
          //         tempPrimary.id = points->at(pCount).primaryID;
          //         tempPrimary.num = 1;
          //         primaryCounter.push_back(tempPrimary);
          //       }
          //     }
          //   }
          // }

          float distance[3];
          Int_t SingleID;
          distance[0] = distance3D(points->at(0).x,
                                   points->at(0).y,
                                   points->at(0).z,
                                   points->at(1).x,
                                   points->at(1).y,
                                   points->at(1).z);
          distance[1] = distance3D(points->at(0).x,
                                   points->at(0).y,
                                   points->at(0).z,
                                   points->at(2).x,
                                   points->at(2).y,
                                   points->at(2).z);
          distance[2] = distance3D(points->at(2).x,
                                   points->at(2).y,
                                   points->at(2).z,
                                   points->at(1).x,
                                   points->at(1).y,
                                   points->at(1).z);
          float minDistance = INFINITY;
          int distanceID;
          for(int iMin = 0 ; iMin < 3 ; iMin++)
          {
            if(distance[iMin] < minDistance)
            {
              minDistance = distance[iMin];
              distanceID = iMin;
            }
          }

          if(distanceID == 0) // 2 is single, 0 and 1 summed
          {
            SingleID = 2;
          }
          if(distanceID == 1) // 1 is single, 0 and 2 summed
          {
            SingleID = 1;
          }
          if(distanceID == 2) // 0 is single, 1 and 2 summed
          {
            SingleID = 0;
          }

          // for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          // {
          //
          // }
          //
          // // if(primaryCounter.size() == 2) // only if events from both primaries
          // // {
          //
          // Int_t SinglePrimaryID;
          //
          // if(primaryCounter[0].num == 1)
          // {
          //   SinglePrimaryID = primaryCounter[0].id;
          // }
          // else
          // {
          //   SinglePrimaryID = primaryCounter[1].id;
          // }

          if(listmodeoutput)
          {
            EventFormat fe;
            fe.ts     = timeCounter;
            fe.random = 0;
            fe.d      = plates_distance;
            fe.yozRot = 0;

            float sum_energy = 0;
            float first_time = INFINITY;
            Int_t firstCrystalID;
            for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
            {
              if(pCount == SingleID)
              {
                // SingleID = pCount;
                //save first gamma
                fe.x1     = points->at(pCount).x;
                fe.y1     = points->at(pCount).y;
                fe.z1     = points->at(pCount).z;
                fe.e1     = points->at(pCount).energy * 1000;  // clearpem reco wants KeV
                fe.n1     = 1;
              }
              else
              {
                //sum energy and find first point
                sum_energy += points->at(pCount).energy;
                if(points->at(pCount).time < first_time)
                {
                  first_time = points->at(pCount).time;
                  firstCrystalID = pCount;
                }
              }
            }

            if(1000.0*points->at(SingleID).energy > enMin && 1000.0*points->at(SingleID).energy < enMax && sum_energy*1000.0 > enMin && 1000.0*sum_energy < enMax) //energy limits
            {
              onlyThreeCrystalsInEnergyWindow++; //baseline sensitivity
            }

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
            fe.x2     = temp_x;
            fe.y2     = temp_y;
            fe.z2     = temp_z;
            fe.e2     = sum_energy * 1000.0;// clearpem reco wants KeV
            fe.n2     = 1;
            // fe.dt     = (double) (points->at(0).time - points->at(1).time);
            fe.dt     = (double) (points->at(SingleID).time - first_time);
            ofs3cry_avg.write((char*)&fe,sizeof(fe));
            //---------------------------------------------//


            //---------------------------------------------//
            //strategy b: assign first hit to greater charge
            //---------------------------------------------//
            Double_t maxEnergy = 0.0;
            Int_t maxEnergyID = -1;
            for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
            {
              if(pCount != SingleID)
              {
                if(points->at(pCount).energy >  maxEnergy)
                {
                  maxEnergy = points->at(pCount).energy;
                  maxEnergyID = pCount;
                }
              }
            }
            fe.x2     = points->at(maxEnergyID).x;
            fe.y2     = points->at(maxEnergyID).y;
            fe.z2     = points->at(maxEnergyID).z;
            //fe.e2     = sum_energy * 1000.0;// clearpem reco wants KeV
            // fe.n2     = 1;
            // fe.dt     = (double) (points->at(0).time - points->at(1).time);
            fe.dt     = (double) (points->at(SingleID).time - points->at(maxEnergyID).time);
            ofs3cry_maxEnergy.write((char*)&fe,sizeof(fe));
            //---------------------------------------------//


            // ofs3cry_effCompton
            //---------------------------------------------//
            //strategy c: compton with efficiency
            //---------------------------------------------//
            Int_t otherID = -1;
            for(int fastI = 0 ; fastI < 3 ; fastI++)
            {
              if(fastI != SingleID && fastI != firstCrystalID)
              otherID = fastI;
            }
            for(int iEff = 0 ; iEff < effSteps ; iEff++)
            {

              float efficiency = min_efficiency + iEff*efficiency_step;
              if(otherID == -1) std::cout << "ERROR!!!!!!!!!!!!" << std::endl;
              // choose the right one with an efficiency
              double dice = (double)rand() / (double)RAND_MAX;
              Int_t chosenCryID;
              if(dice <= efficiency)
              {
                chosenCryID = firstCrystalID;
              }
              else
              {
                chosenCryID = otherID;
              }

              fe.x2     = points->at(chosenCryID).x;
              fe.y2     = points->at(chosenCryID).y;
              fe.z2     = points->at(chosenCryID).z;
              fe.dt     = (double) (points->at(SingleID).time - points->at(chosenCryID).time);
              ofs3cry_effCompton[iEff]->write((char*)&fe,sizeof(fe));
            }
          }






          // indentification of the pair. at least one must be in the "511 window", i.e. should have experienced photoelectric deposition by one of the 2 gammas
          // bool foundSingle511 = false;
          // int SingleID = -1;
          // for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          // {
          //   if(1000.0*points->at(pCount).energy > enMin && 1000.0*points->at(pCount).energy < enMax)
          //   {
          //     foundSingle511 = true;
          //     SingleID = pCount;    //take the ID of the single 511 deposition
          //   }
          // }
          // if(foundSingle511) //now take the other two
          // {
          //   float sum_energy = 0;
          //   for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          //   {
          //     if(pCount != SingleID) //only the others
          //     {
          //       //sum the energies
          //       sum_energy += points->at(pCount).energy;
          //     }
          //   }
          //   //check if sum is in window
          //   if(1000.0*sum_energy > enMin && 1000.0*sum_energy < enMax)
          //   {
          //     //count this event as 3 crystals only, in energy window
          //     onlyThreeCrystalsInEnergyWindow++;
          //     //---------------------------------------------//
          //     //strategy a: calculate an average event
          //     //---------------------------------------------//
          //     float temp_x = 0;
          //     float temp_y = 0;
          //     float temp_z = 0;
          //     for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          //     {
          //       if(pCount != SingleID)
          //       {
          //         //average x,y,z
          //         temp_x += (points->at(pCount).x * points->at(pCount).energy) / sum_energy;
          //         temp_y += (points->at(pCount).y * points->at(pCount).energy) / sum_energy;
          //         temp_z += (points->at(pCount).z * points->at(pCount).energy) / sum_energy;
          //       }
          //     }
          //     //save in another file
          //     if(binary)
          //     {
          //       pair.x1 = points->at(SingleID).x;
          //       pair.y1 = points->at(SingleID).y;
          //       pair.z1 = points->at(SingleID).z;
          //       pair.x2 = temp_x;
          //       pair.y2 = temp_y;
          //       pair.z2 = temp_z;
          //       ofs3cry_avg.write((char*)&pair,sizeof(pair));
          //     }
          //     else
          //     {
          //       ofs3cry_avg << points->at(SingleID).x << " "
          //       << points->at(SingleID).y << " "
          //       << points->at(SingleID).z << " "
          //       << temp_x << " "
          //       << temp_y << " "
          //       << temp_z
          //       << std::endl;
          //     }
          //
          //
          //     //---------------------------------------------//
          //     //strategy b: magical compton
          //     //---------------------------------------------//
          //     //check which crystal was hit first
          //     Int_t firstCrystalID;
          //     float timeOfFirst = INFINITY;
          //     for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          //     {
          //       if(pCount != SingleID)
          //       {
          //         if(points->at(pCount).time <  timeOfFirst)
          //         {
          //           timeOfFirst = points->at(pCount).time;
          //           firstCrystalID = pCount;
          //         }
          //       }
          //     }
          //
          //     //now save only the single 511 and the position of first between the two comptons (100% efficient compton discrimination)
          //     if(binary)
          //     {
          //       pair.x1 = points->at(SingleID).x;
          //       pair.y1 = points->at(SingleID).y;
          //       pair.z1 = points->at(SingleID).z;
          //       pair.x2 = points->at(firstCrystalID).x;
          //       pair.y2 = points->at(firstCrystalID).y;
          //       pair.z2 = points->at(firstCrystalID).z;
          //       ofs3cry_magicalCompton.write((char*)&pair,sizeof(pair));
          //     }
          //     else
          //     {
          //       ofs3cry_magicalCompton << points->at(SingleID).x << " "
          //       << points->at(SingleID).y << " "
          //       << points->at(SingleID).z << " "
          //       << points->at(firstCrystalID).x << " "
          //       << points->at(firstCrystalID).y << " "
          //       << points->at(firstCrystalID).z
          //       << std::endl;
          //     }
          //
          //
          //     //---------------------------------------------//
          //     //strategy c: compton with efficiency
          //     //---------------------------------------------//
          //     //get the remaing id among 0,1,2
          //
          //     Int_t otherID = -1;
          //     for(int fastI = 0 ; fastI < 3 ; fastI++)
          //     {
          //       if(fastI != SingleID && fastI != firstCrystalID)
          //       otherID = fastI;
          //     }
          //     if(otherID == -1) std::cout << "ERROR!!!!!!!!!!!!" << std::endl;
          //
          //
          //     // choose the right one with an efficiency
          //     double dice = (double)rand() / (double)RAND_MAX;
          //     Int_t chosenCryID;
          //
          //     if(dice <= efficiency)
          //     {
          //       chosenCryID = firstCrystalID;
          //     }
          //     else
          //     {
          //       chosenCryID = otherID;
          //     }
          //
          //
          //     if(binary)
          //     {
          //       pair.x1 = points->at(SingleID).x;
          //       pair.y1 = points->at(SingleID).y;
          //       pair.z1 = points->at(SingleID).z;
          //       pair.x2 = points->at(chosenCryID).x;
          //       pair.y2 = points->at(chosenCryID).y;
          //       pair.z2 = points->at(chosenCryID).z;
          //       ofs3cry_effCompton.write((char*)&pair,sizeof(pair));
          //     }
          //     else
          //     {
          //       ofs3cry_effCompton << points->at(SingleID).x    << " "
          //       << points->at(SingleID).y    << " "
          //       << points->at(SingleID).z    << " "
          //       << points->at(chosenCryID).x << " "
          //       << points->at(chosenCryID).y << " "
          //       << points->at(chosenCryID).z
          //       << std::endl;
          //     }
          //
          //     //---------------------------------------------//
          //     //strategy d: assign first hit to greater charge
          //     //---------------------------------------------//
          //     Double_t maxEnergy = 0.0;
          //     Int_t maxEnergyID = -1;
          //     for(Int_t pCount = 0 ; pCount < points->size() ; pCount++)
          //     {
          //       if(pCount != SingleID)
          //       {
          //         if(points->at(pCount).energy >  maxEnergy)
          //         {
          //           maxEnergy = points->at(pCount).energy;
          //           maxEnergyID = pCount;
          //         }
          //       }
          //     }
          //
          //     if(binary)
          //     {
          //       pair.x1 = points->at(SingleID).x;
          //       pair.y1 = points->at(SingleID).y;
          //       pair.z1 = points->at(SingleID).z;
          //       pair.x2 = points->at(maxEnergyID).x;
          //       pair.y2 = points->at(maxEnergyID).y;
          //       pair.z2 = points->at(maxEnergyID).z;
          //       ofs3cry_maxEnergy.write((char*)&pair,sizeof(pair));
          //     }
          //     else
          //     {
          //       ofs3cry_maxEnergy << points->at(SingleID).x    << " "
          //       << points->at(SingleID).y    << " "
          //       << points->at(SingleID).z    << " "
          //       << points->at(maxEnergyID).x << " "
          //       << points->at(maxEnergyID).y << " "
          //       << points->at(maxEnergyID).z
          //       << std::endl;
          //     }
          //   }
          //   else // the two scatter don't sum up to enter 511 window
          //   {
          //     //don't save the event for "reconstruction"
          //   }
          // }
          // else // no crystal among the 3 has a deposition that enters the 511 window
          // {
          //   //don't save the event for "reconstruction"
          // }

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
          if(1000.0*sum_energy > 2.0*enMin && 1000.0*sum_energy < 2.0*enMax)
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
          if(1000.0*sum_energy > 2.0*enMin && 1000.0*sum_energy < 2.0*enMax)
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
  std::cout << "Total interactions with det.           = " << nsamples << std::endl;
  std::cout << "0 crystals hit                         = " << NoCrystals <<  std::endl;
  std::cout << "1 crystal hit                          = " << onlyOneCrystals <<  std::endl;
  std::cout << "1 crystal hit - in energy. window      = " << onlyOneCrystalsInEnergyWindow <<  std::endl;
  std::cout << "2 crystals hit                         = " << onlyTwoCrystals <<  std::endl;
  std::cout << "2 crystals same primary                = " << twoCrystalsSamePrimary <<  std::endl;
  std::cout << "2 cry same primary - in energy. window = " << onlyTwoCrystalsHitBySameGammaInEnergyWindow <<  std::endl;
  std::cout << "2 crystals diff primary                = " << twoCrystalsDiffPrimary <<  std::endl;
  std::cout << "2 cry diff primary - in energy. window = " << onlyTwoCrystalsInEnergyWindow <<  std::endl;
  std::cout << "3 crystals hit                         = " << onlyThreeCrystals <<  std::endl;
  std::cout << "3 crystals hit - in energy. window     = " << onlyThreeCrystalsInEnergyWindow <<  std::endl;
  std::cout << "4 crystals hit                         = " << onlyFourCrystals <<  std::endl;
  std::cout << "4 crystals hit - in energy. window     = " << onlyFourCrystalsInEnergyWindow <<  std::endl;
  std::cout << "More than 4 crystals hit               = " << moreThanFourCrystals <<  std::endl;
  std::cout << "More than 4 - in energy. window        = " << moreThanFourCrystalsInEnergyWindow <<  std::endl;


  //SENSITIVITY COUNTS
  summary_out << "Total interactions with det.           = " << nsamples << std::endl;
  summary_out << "0 crystals hit                         = " << NoCrystals <<  std::endl;
  summary_out << "1 crystal hit                          = " << onlyOneCrystals <<  std::endl;
  summary_out << "1 crystal hit - in energy. window      = " << onlyOneCrystalsInEnergyWindow <<  std::endl;
  summary_out << "2 crystals hit                         = " << onlyTwoCrystals <<  std::endl;
  summary_out << "2 crystals same primary                = " << twoCrystalsSamePrimary <<  std::endl;
  summary_out << "2 cry same primary - in energy. window = " << onlyTwoCrystalsHitBySameGammaInEnergyWindow <<  summary_out;
  summary_out << "2 crystals diff primary                = " << twoCrystalsDiffPrimary <<  std::endl;
  summary_out << "2 cry diff primary - in energy. window = " << onlyTwoCrystalsInEnergyWindow <<  std::endl;
  summary_out << "3 crystals hit                         = " << onlyThreeCrystals <<  std::endl;
  summary_out << "3 crystals hit - in energy. window     = " << onlyThreeCrystalsInEnergyWindow <<  std::endl;
  summary_out << "4 crystals hit                         = " << onlyFourCrystals <<  std::endl;
  summary_out << "4 crystals hit - in energy. window     = " << onlyFourCrystalsInEnergyWindow <<  std::endl;
  summary_out << "More than 4 crystals hit               = " << moreThanFourCrystals <<  std::endl;
  summary_out << "More than 4 - in energy. window        = " << moreThanFourCrystalsInEnergyWindow <<  std::endl;
  // std::cout << "Only 1 crystal   [all - in energy window]   = " << onlyOneCrystals   << " - " << onlyOneCrystalsInEnergyWindow   << std::endl;
  // std::cout << "Only 2 crystals  [all - in energy window]   = " << onlyTwoCrystals      << " - " << onlyTwoCrystalsInEnergyWindow   <<  std::endl;
  //
  // std::cout <<
  //
  // onlyTwoCrystalsHitBySameGamma
  // onlyTwoCrystalsHitBySameGammaInEnergyWindow
  //
  // "onlyTwoCrystalsFromSameGammaInEnergyWindow   = " << onlyTwoCrystalsFromSameGammaInEnergyWindow  <<  std::endl;
  //
  //
  // std::cout << "Only 3 crystals  [all - in energy window]   = " << onlyThreeCrystals    << " - " << onlyThreeCrystalsInEnergyWindow <<  std::endl;
  // std::cout << "Only 4 crystals                             = " << onlyFourCrystals     << std::endl;
  // std::cout << "More than 4 crystals                        = " << moreThanFourCrystals << std::endl;
  //
  //
  // summary_out << "Total dataset                               = " << nsamples << std::endl;
  // summary_out << "Only 1 crystal   [all - in energy window]   = " << onlyOneCrystals   << " - " << onlyOneCrystalsInEnergyWindow   << std::endl;
  // summary_out << "Only 2 crystals  [all - in energy window]   = " << onlyTwoCrystals      << " - " << onlyTwoCrystalsInEnergyWindow   <<  std::endl;
  // summary_out << "Only 3 crystals  [all - in energy window]   = " << onlyThreeCrystals    << " - " << onlyThreeCrystalsInEnergyWindow <<  std::endl;
  // summary_out << "Only 4 crystals                             = " << onlyFourCrystals     << std::endl;
  // summary_out << "More than 4 crystals                        = " << moreThanFourCrystals << std::endl;

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
  summary_out.close();
  // ofs3cry_magicalCompton.close();
  for(int i = 0 ; i < effSteps ; i++)
  {
    ofs3cry_effCompton[i]->close();
  }
  ofs3cry_maxEnergy.close();
  return 0;
}
