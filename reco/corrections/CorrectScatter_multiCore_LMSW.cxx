/*======================================================================
  Simulation for Scatter Correction for the ClearPEM scanner
  
  Inputs: Input segmented image, outputPrefix, min. number of detected events,
          Crystals distance, lower limit for energy window, total number of 
          non-random events on the real acquisition. 
       
  Outputs: LMF of simulated events (scattered events scaled to real acq.)

  Author: C. S. Ferreira (claudiaf@lip.pt)

  Modification History: April 2012 (final version)
                        February 2013 (adapted to ListMode reconstruction)
			April 2013 speed related improvements

  ======================================================================*/
/* This code implementation was developed in the context of the 
   ClearPEM project and the Crystal Clear collaboration
// Copyright Claudia Ferreira
// 2013 IBEB, Lisbon and LIP, Lisbon
========================================================================*/



#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <string>

#include <boost/random.hpp>
#include <boost/nondet_random.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;

#include <pthread.h>
#define nThreads 8

//ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"



 //----------SEGMENTED IMAGE OPTION--------------//
typedef float                  PixelType;
const unsigned int             Dimension = 3;

typedef itk::Image< PixelType, Dimension >   ImageType;

typedef itk::ImageFileReader< ImageType >  ReaderType;
ReaderType::Pointer reader = ReaderType::New();


struct InfoFile{ 
    float c0;
    float c1; 
    float c2;
    float c3;
    float c4;
    float c5;
    float c6;
    float c7; } __attribute__((__packed__));


struct ThreadData{
  double nFinalEvents;
  InfoFile * crossSecInfo;
  // FILE* filePrompt;
  FILE* fileCompt;
  //  FILE* fileDirect;
  int nphoton;
  float xDim;
  float yDim; 
  float axialDim; 
  ImageType::SpacingType imageSpacing; 
  ImageType::SizeType size;
  float half_plates_distance;
  float Min_detect_threshold;
  double nEvent_45;
  double nEvent0;
  double nEvent45;
  double nEvent90;
};
 

ThreadData threadDataArray[nThreads]; 

void* doJob(void*targ);

InfoFile * readDataFile();


float plane(float p0, float z0, float w);


struct FinalCoordDetectAxis{float x; float y; float z; float u; float v; float w;}__attribute__((__packed__));
struct PhotonCoord{float x; float y; float z;}__attribute__((__packed__));

PhotonCoord detector(float half_plates_distance, float axialLength, float transaxialLength, FinalCoordDetectAxis finalCoordDetect);


struct ScatterCoord{float x; float y; float z; float u; float v; float w; float E; int nCompton;}__attribute__((__packed__));


struct Coord{float x; float y; float z; float u; float v; float w; float E; float THETA; float PHI;}__attribute__((__packed__));

ScatterCoord scatter(Coord coord, InfoFile * crossSecInfo, float cutOff, boost::mt19937 *generator, int nphoton, float xDim, float yDim, float axialDim, ImageType::SpacingType imageSpacing, ImageType::SizeType size, float Min_Ethresh);


std::vector<float> compton(float E0, boost::mt19937 *generator);
std::vector<float> crossSection(float E, InfoFile * crossSecInfo);
float cylinder(float rCylin, float x0, float y0, float u, float v);

//VARS
float E = 511;           //(keV) Initial Energy of the photon
float cutOff = 70;      //(keV) cut-off energy
float axialLength =  173.5;        //(mm)  detector plates length
float transaxialLength = 152.2;    //(mm)  detector plates width (transaxial dimension)
float Max_detect_threshold = 650;      //(keV) Energy window (max) 
double nSimulatedEvents = 0;



struct Pos{  
    float x1; 
    float y1;
    float z1;
    float x2;
    float y2;
    float z2; 
  }  __attribute__((__packed__));
      
  Pos pos;


struct LmEvent{    //DKFZFormat
  float x1; 
  float y1;
  float z1;
  float x2;
  float y2;
  float z2;
  float angle;
  float weight; }  __attribute__((__packed__));





int main(int argc, char *argv[])
{
  
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  InputSegmentedImage.nii  outputBinaryFile  nMinEvents  CRYSTALS_distance  low-energy  nRealPromptEv" << std::endl;
    return EXIT_FAILURE;
    }
  
  
  srand ( time(NULL) );
 
  //TIME STAMP
  time_t time1,time2;
  time(&time1);


  //----------SEGMENTED IMAGE OPTION--------------//
  const char * fileName = argv[1];
  reader->SetFileName( fileName );
  reader->Update();


  // get size of the whole 3D image
  ImageType::RegionType inputSize = reader->GetOutput()->GetLargestPossibleRegion();
  ImageType::SizeType size = inputSize.GetSize();
  std::cerr << "Input size (voxels): " << size << std::endl;
  float imageSizeX = size[0];
  float imageSizeY = size[1];
  float imageSizeZ = size[2];
 

  
  ImageType::SpacingType imageSpacing = reader->GetOutput()->GetSpacing();
  std::cerr << "Voxel Size: " << imageSpacing << std::endl;

  float axialDim = imageSizeZ*imageSpacing[2];  //mm 
  float xDim = imageSizeX*imageSpacing[0];  //mm
  float yDim = imageSizeY*imageSpacing[1];  //mm
  std::cerr<< "Dimensions (mm):  " <<xDim<<", "<<yDim<<", "<<axialDim<<std::endl;
  


  //----- FILES OUT --------//
  /*FILE * out = NULL;
  std::string outFile = std::string(argv[2]) + ".lm";
  const char* fileOut = outFile.c_str();
  if(strcmp(argv[2], "-") == 0) 
    out = stdout;
  else
    out = fopen(fileOut, "w");
  if (out == NULL) 
    {
      fprintf(stderr, "Can't open output file \n");
      return -1;
    }
  */
  std::string comptonEvFile = std::string(argv[2]) + "_comptonEv.lm";
  //std::string directEvFile = std::string(argv[2]) + "_directEv.lm";
  const char* comptonFileName = comptonEvFile.c_str();
  //const char* directFileName = directEvFile.c_str();
  FILE * out2 = fopen(comptonFileName, "w");
  //FILE * out3 = fopen(directFileName,"w");
  std::string comptonScaledEvFile = std::string(argv[2]) + "_comptonEvSc.lm";
  const char* comptonFileNameSc = comptonScaledEvFile.c_str();
  FILE * outScatterScaled = fopen(comptonFileNameSc, "w");
  

 
  //Vars
  double nFinalEvents = atol(argv[3]);
  // This number shall never be less than 10e8
  std::cerr << "Minimum total number of detected events:  " << nFinalEvents<< std::endl;
  // user shall give the distance with the offset included (mm) Sonic (+59.06); ClearPEM (+64.3) //**********
  float half_plates_distance = (atof(argv[4]) )/2;
  std::cerr << "half plates distance= "<<  half_plates_distance << std::endl;
  float Min_detect_threshold = atof(argv[5]);     //(keV) Energy window (min)
  std::cerr << "Lower threshold for Energy window: " << Min_detect_threshold << std::endl;
  int nphoton = 0; //check purposes 


  InfoFile * crossSecInfo = readDataFile();


  
  //DO THE JOB
  double nTotalDetectedEv = 0;
  double nEventsDetected__45deg = 0;
  double nEventsDetected_0deg = 0;
  double nEventsDetected_45deg = 0;
  double nEventsDetected_90deg = 0;


  pthread_t threads[nThreads];
  for(int i = 0; i<nThreads;i++) 
    {
      threadDataArray[i].nFinalEvents = nFinalEvents/nThreads;
      threadDataArray[i].crossSecInfo = crossSecInfo;
      //threadDataArray[i].filePrompt = out;
      threadDataArray[i].fileCompt = out2;
      //  threadDataArray[i].fileDirect = out3;
      threadDataArray[i].nphoton = nphoton;
      threadDataArray[i].xDim = xDim;
      threadDataArray[i].yDim = yDim;
      threadDataArray[i].axialDim = axialDim;
      threadDataArray[i].imageSpacing = imageSpacing; 
      threadDataArray[i].size = size;
      threadDataArray[i].half_plates_distance = half_plates_distance;
      threadDataArray[i].Min_detect_threshold = Min_detect_threshold;
      threadDataArray[i].nEvent_45 = 0;
      threadDataArray[i].nEvent0 = 0;
      threadDataArray[i].nEvent45 = 0;
      threadDataArray[i].nEvent90 = 0;

      int rc = pthread_create(&threads[i], NULL, doJob, (void *) &threadDataArray[i]);     
    }

  for(int t=0; t<nThreads; t++) 
    {
      int rc = pthread_join(threads[t], NULL);
      nEventsDetected__45deg = nEventsDetected__45deg + threadDataArray[t].nEvent_45;
      nEventsDetected_0deg = nEventsDetected_0deg + threadDataArray[t].nEvent0;
      nEventsDetected_45deg = nEventsDetected_45deg + threadDataArray[t].nEvent45;
      nEventsDetected_90deg = nEventsDetected_90deg + threadDataArray[t].nEvent90;

      if (rc) 
	{
	  printf("ERROR; return code from pthread_join() is %d\n", rc);
	  exit(-1);
	}
      
    }
  
  
  nTotalDetectedEv = nEventsDetected__45deg + nEventsDetected_0deg + nEventsDetected_45deg + nEventsDetected_90deg;
  std::cerr << "\n" << "Number of simulated events inside object: " << nSimulatedEvents << std::endl;
  std::cerr << "Number of detected events for -45 deg: " << nEventsDetected__45deg << std::endl;
  std::cerr << "Number of detected events for 0 deg: " << nEventsDetected_0deg << std::endl;
  std::cerr << "Number of detected events for 45 deg: " << nEventsDetected_45deg << std::endl;
  std::cerr << "Number of detected events for 90 deg: " << nEventsDetected_90deg << std::endl;
  std::cerr << "Total number of detected events: " <<  nTotalDetectedEv << std::endl;


  //-----TIME STAMP----//
  time(&time2);
  double diff_sec = difftime (time2,time1);
  std::cout << "Time required for execution (simulation): " << diff_sec << " seconds."<<std::endl;
  //------------------//
  

  //fclose(out);
  fclose(out2);
  //fclose(out3);


  FILE * outScatter = fopen(comptonFileName, "r");
  LmEvent eventLMsc;
  double nPromptEv;
  sscanf(argv[6],"%lf",&nPromptEv);
  
  int sWeight = 0;
  while(fread((void*)&eventLMsc,sizeof(eventLMsc),1,outScatter)==1)
    {
      eventLMsc.weight = (float)(-nPromptEv/nTotalDetectedEv);
      if(sWeight == 0){
	std::cout<<"Scaling weight: "<<eventLMsc.weight<< std::endl;
	sWeight=1;
      }
      fwrite((void*) &eventLMsc, sizeof(eventLMsc), 1, outScatterScaled);
    }
  
  fclose(outScatter);
  fclose(outScatterScaled);

  return 0;
  
}

void* doJob(void*targ)
{
  ThreadData *threadData = (ThreadData *)targ;

  double nFinalEvents = threadData->nFinalEvents;
  InfoFile *crossSecInfo = threadData->crossSecInfo;
  // FILE*out = threadData->filePrompt;
  FILE*out2 = threadData->fileCompt;
  // FILE*out3 = threadData->fileDirect;
  int nphoton = threadData->nphoton;
  float xDim = threadData->xDim;
  float yDim = threadData->yDim;
  float axialDim = threadData->axialDim;
  ImageType::SpacingType imageSpacing = threadData->imageSpacing;
  ImageType::SizeType size = threadData->size;
  float half_plates_distance = threadData->half_plates_distance;
  float Min_detect_threshold = threadData->Min_detect_threshold;


  int nThreadDetectedEv = 0;

  //random numbers
  size_t randResult;
  FILE * fileRandomSeed = fopen("/dev/urandom","r");
  time_t seed;
  randResult = fread((void*)&seed, sizeof(time_t), 1, fileRandomSeed); // FIXME unused value.
  
  if (randResult != 1){std::cerr << "error reading random value" << std::endl;}
  mt19937 generator;
  generator.seed(time_t(seed));
  boost::uniform_real<> range(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > next(generator, range);
  


  double nevents = 1e6; 
  while (nThreadDetectedEv < nFinalEvents)
    {
      for(int j=-45; j<91; j=j+45)
	{
	  
	  double nEventsDetected = 0;
	  
	  for(double i=0; i<nevents; i++)
	    {
	      
	      Coord coord;
	      FinalCoordDetectAxis finalCoordDetectAxis;
	      struct DetectCoord{float x1; float y1; float z1; float x2; float y2; float z2;}detectCoord;
	      
		
	      //GENERATE PARTICLES
	      
	      //Generation inside FOV's dimensions
	      float x0 = (next() * xDim) - (xDim/2);
	      float y0 = (next() * yDim) - (yDim/2);
	      float z0 = (next() * axialDim) - (axialDim/2);
		
		
	      //Check if event is inside object's boundaries
	      ImageType::IndexType pixIndex; 	  
	      pixIndex[0] = (round(x0) + ((xDim)/2))/(imageSpacing[0]);
	      pixIndex[1] = (round(y0) + ((yDim)/2))/(imageSpacing[1]);
	      pixIndex[2] = (round(z0) + ((axialDim)/2))/(imageSpacing[2]);
		

	      float imageSizeX = size[0];
	      float imageSizeY = size[1];
	      float imageSizeZ = size[2];
	      if (pixIndex[0]> imageSizeX || pixIndex[1]> imageSizeY || pixIndex[2]> imageSizeZ){std::cerr<<"Check events' generation" << std::endl;}
		
		
	      ImageType::PixelType pixValue = reader->GetOutput()->GetPixel(pixIndex); 
		
	      if (pixValue <=0) {continue;}
	      nSimulatedEvents++;
		
		
	      //Spherical coordinates - particles can travel along any direction 
	      float thetaSphericCoord = M_PI * next();
	      float phi = 2 * M_PI * next();
		
	      float u0 = sin(thetaSphericCoord) * cos(phi);  //x-axis direction cosine of the particle
	      float v0 = sin(thetaSphericCoord) * sin(phi);  //y-axis direction cosine of the particle
	      float w0 = cos(thetaSphericCoord);             //z-axis direction cosine of the particle 
		
	      coord.x = x0;
	      coord.y = y0;
	      coord.z = z0;
	      coord.u = u0;
	      coord.v = v0;
	      coord.w = w0;
	      coord.E = E;
	      coord.THETA = thetaSphericCoord;
	      coord.PHI = phi;
		

		
	      //--------final coordinates-----//
	      threadData->nphoton++;
	      ScatterCoord photon1;
	      //photon1.x=x0; photon1.y=y0; photon1.z=z0; photon1.u=u0; photon1.v=v0; photon1.w=w0; photon1.E=511; photon1.nCompton=0;  //for 511 keV only
	      photon1 = scatter(coord,crossSecInfo, cutOff, &generator, nphoton, xDim, yDim, axialDim, imageSpacing, size, Min_detect_threshold);
		
	      //Energy window for the detector
	      if(photon1.E >= Max_detect_threshold || photon1.E <= Min_detect_threshold) {continue;}
		
		
	      //conversion of coordinates: simulation (z is the axial) to the detector referential (x is the axial)
	      // coordinates -> ClearPEM coordinates
	      // vertical -> YY = -YYcp
	      // horizontal -> XX = ZZcp
	      // scanner axis -> ZZ = -XXcp
	      finalCoordDetectAxis.x = -photon1.z;
	      finalCoordDetectAxis.y = -photon1.y;
	      finalCoordDetectAxis.z = photon1.x;
		
	      finalCoordDetectAxis.u = -photon1.w;
	      finalCoordDetectAxis.v = -photon1.v;
	      finalCoordDetectAxis.w = photon1.u;
		
		
	      //Rotation of coordinates for detection
	      /*--------------------------------------------
		Rotation Matrix
		1       0             0
		0   cos(rotation_deg)   -sin(rotation_deg)
		0   sin(rotation_deg)   cos(rotation_deg)
		--------------------------------------------*/   	  
	      float newFinalCoordDetectAxis_y = finalCoordDetectAxis.y * cos(j*((2*M_PI)/360)) + finalCoordDetectAxis.z * (-sin (j*((2*M_PI)/360))); //radians
	      float newFinalCoordDetectAxis_z = finalCoordDetectAxis.y * sin(j*((2*M_PI)/360)) + finalCoordDetectAxis.z * cos(j*((2*M_PI)/360));
	      float newFinalCoordDetectAxis_v = finalCoordDetectAxis.v * cos(j*((2*M_PI)/360)) + finalCoordDetectAxis.w * (-sin (j*((2*M_PI)/360)));
	      float newFinalCoordDetectAxis_w = finalCoordDetectAxis.v * sin(j*((2*M_PI)/360)) + finalCoordDetectAxis.w * cos(j*((2*M_PI)/360));
		
	      finalCoordDetectAxis.y = newFinalCoordDetectAxis_y;
	      finalCoordDetectAxis.z = newFinalCoordDetectAxis_z;
	      finalCoordDetectAxis.v = newFinalCoordDetectAxis_v;
	      finalCoordDetectAxis.w = newFinalCoordDetectAxis_w;
		
	      //final coordinates (detector)
	      PhotonCoord detected1 = detector(half_plates_distance, axialLength, transaxialLength, finalCoordDetectAxis);
	      detectCoord.x1 = detected1.x;
	      detectCoord.y1 = detected1.y;
	      detectCoord.z1 = detected1.z;
		
	      if (detectCoord.z1 == 0){continue;}
	
	      //+++++++++++++++++++++++++++++++
	      //++         photon 2          ++
	      //+++++++++++++++++++++++++++++++
		
	      //Same parameters as photon1 except for the direction
	      coord.u = -u0;
	      coord.v = -v0;
	      coord.w = -w0;
		
	      nphoton++;
	      ScatterCoord photon2;
	      // photon2.x=x0; photon2.y=y0; photon2.z=z0; photon2.u=-u0; photon2.v=-v0; photon2.w=-w0; photon2.E=511; photon2.nCompton=0;  //for 511 keV only
	      photon2 = scatter(coord,crossSecInfo, cutOff, &generator, nphoton, xDim, yDim, axialDim, imageSpacing, size, Min_detect_threshold);
		
		
		
	      //Energy window for the detector
	      if(photon2.E >= Max_detect_threshold || photon2.E <= Min_detect_threshold) {continue;}
		
	      //conversion of coordinates to the detector axis
	      finalCoordDetectAxis.x = -photon2.z;
	      finalCoordDetectAxis.y = -photon2.y;
	      finalCoordDetectAxis.z = photon2.x;
	      finalCoordDetectAxis.u = -photon2.w;
	      finalCoordDetectAxis.v = -photon2.v;
	      finalCoordDetectAxis.w = photon2.u;
		
		
	      //Rotation of coordinates for detection 	  
	      newFinalCoordDetectAxis_y = finalCoordDetectAxis.y * cos(j*((2*M_PI)/360)) + finalCoordDetectAxis.z * (-sin (j*((2*M_PI)/360))); //radians
	      newFinalCoordDetectAxis_z = finalCoordDetectAxis.y * sin(j*((2*M_PI)/360)) + finalCoordDetectAxis.z * cos(j*((2*M_PI)/360));
	      newFinalCoordDetectAxis_v = finalCoordDetectAxis.v * cos(j*((2*M_PI)/360)) + finalCoordDetectAxis.w * (-sin (j*((2*M_PI)/360)));
	      newFinalCoordDetectAxis_w = finalCoordDetectAxis.v * sin(j*((2*M_PI)/360)) + finalCoordDetectAxis.w * cos(j*((2*M_PI)/360));
		
	      finalCoordDetectAxis.y = newFinalCoordDetectAxis_y;
	      finalCoordDetectAxis.z = newFinalCoordDetectAxis_z;
	      finalCoordDetectAxis.v = newFinalCoordDetectAxis_v;
	      finalCoordDetectAxis.w = newFinalCoordDetectAxis_w;
		
		
	      //final coordinates (detector)
	      PhotonCoord detected2 = detector(half_plates_distance, axialLength, transaxialLength, finalCoordDetectAxis);
	      detectCoord.x2 = detected2.x;
	      detectCoord.y2 = detected2.y;
	      detectCoord.z2 = detected2.z;
	       
	      if (detectCoord.z2 == 0){continue;}
		
	      if (((detectCoord.z1 <= 0) && (detectCoord.z2 <= 0)) || ((detectCoord.z1 >= 0) && (detectCoord.z2 >= 0))){continue;}
		
	      nEventsDetected ++;
		
	      if(detectCoord.x1 > ((axialLength/2) - 2.7109375/2.)) {detectCoord.x1 = (axialLength/2)  -2.7109375/2.;}
	      if(detectCoord.x1 < ((-axialLength/2) + 2.7109375/2.)) {detectCoord.x1 = (-axialLength/2) +2.7109375/2.;}
	      if(detectCoord.x2 > ((axialLength/2) - 2.7109375/2.)) {detectCoord.x2 = (axialLength/2) -2.7109375/2.;}
	      if(detectCoord.x2 < ((-axialLength/2) + 2.7109375/2.)) {detectCoord.x2 = (-axialLength/2) +2.7109375/2.;}
		
		
	      if(detectCoord.y1 > ((transaxialLength/2) -3.1708333/2.)) {detectCoord.y1 = (transaxialLength/2) -3.1708333/2.;}
	      if(detectCoord.y1 < ((-transaxialLength/2) +3.1708333/2.)) {detectCoord.y1 = (-transaxialLength/2) +3.1708333/2.;}
	      if(detectCoord.y2 > ((transaxialLength/2) -3.1708333/2.)) {detectCoord.y2 = (transaxialLength/2) -3.1708333/2.;}
	      if(detectCoord.y2 < ((-transaxialLength/2) +3.1708333/2.)) {detectCoord.y2 = (-transaxialLength/2) +3.1708333/2.;}
		
		
		
	      /*Control for strange values*/
	      if((detectCoord.x1 < (-axialLength/2)) || (detectCoord.x1 > (axialLength/2)) || (detectCoord.y1 < (-transaxialLength/2)) || (detectCoord.y1 > (transaxialLength/2)) || 
		 (detectCoord.z1 < (-half_plates_distance-0.001)) || (detectCoord.z1 > (half_plates_distance+0.001)) || (detectCoord.x2 < (-axialLength/2)) || (detectCoord.x2 > (axialLength/2)) || 
		 (detectCoord.y2 < (-transaxialLength/2)) || (detectCoord.y2 > (transaxialLength)) || (detectCoord.z2 < (-half_plates_distance-0.001)) || (detectCoord.z2 > (half_plates_distance+0.001))) 
		{
		  std::cerr<< "error!!! "<<  detectCoord.x1 << " ; "<< detectCoord.y1<<" ; "<< detectCoord.z1<<" ; "<< detectCoord.x2<<" ; "<< detectCoord.y2<<" ; "<< detectCoord.z2<<std::endl;
		}
		
		
		
	      /*--------------------------------------------
		Rotation for N angles
		--------------------------------------------
		1       0             0
		0   cos(rotation_deg)   sin(rotation_deg)
		0   -sin(rotation_deg)   cos(rotation_deg)
		--------------------------------------------*/   
		
	      float newDetectCoord_y1 = detectCoord.y1 * cos(j*((2*M_PI)/360)) + detectCoord.z1 * sin (j*((2*M_PI)/360)); //radians
	      float newDetectCoord_z1 = -detectCoord.y1 * sin(j*((2*M_PI)/360)) + detectCoord.z1 * cos(j*((2*M_PI)/360));
	      float newDetectCoord_y2 = detectCoord.y2 * cos(j*((2*M_PI)/360)) + detectCoord.z2 * sin (j*((2*M_PI)/360)); 
	      float newDetectCoord_z2 = -detectCoord.y2 * sin(j*((2*M_PI)/360)) + detectCoord.z2 * cos(j*((2*M_PI)/360));
		
		
	      LmEvent eventLM;
		
	      eventLM.x1 = detectCoord.x1;
	      eventLM.y1 = newDetectCoord_y1;
	      eventLM.z1 = newDetectCoord_z1;
	      eventLM.x2 = detectCoord.x2;
	      eventLM.y2 = newDetectCoord_y2;
	      eventLM.z2 = newDetectCoord_z2;
	      eventLM.angle = (j*((2*M_PI)/360));
	      eventLM.weight = 1;

		
	      //  if((photon1.e >= Min_detect_threshold) && (photon2.e >= Min_detect_threshold))
	      //fwrite((void*) &eventLM, sizeof(eventLM), 1, out); 
		
	      //Saves events where at least one of the photons suffered at least 1 scatter.
	      if ((photon1.nCompton >= 1) || (photon2.nCompton >= 1)) 
		{ fwrite((void*) &eventLM, sizeof(eventLM), 1, out2);}
		
	      //Saves only direct unscattered events.
	      //if ((photon1.nCompton == 0) && (photon2.nCompton == 0)) 
	      //{ fwrite((void*) &eventLM, sizeof(eventLM), 1, out3); }
		
		
	    }

	    
	
	  if (j == -45) {threadData->nEvent_45 = threadData->nEvent_45 + nEventsDetected;}
	  if (j == 0) {threadData->nEvent0 = threadData->nEvent0 + nEventsDetected;}
	  if (j == 45) {threadData->nEvent45 = threadData->nEvent45 + nEventsDetected;}
	  if (j == 90) {threadData->nEvent90 = threadData->nEvent90 + nEventsDetected;}
	    
// 	  std::cerr << ">";
	    
	  nThreadDetectedEv += nEventsDetected;
	}
    }
  fclose(fileRandomSeed);

  return NULL;
}


InfoFile * readDataFile()
{
  //--------------------------------
  //File cross_sectionH2O (ASCII adapted from: http://physics.nist.gov/PhysRefData/Xcom/Text/XCOM.html)
  //--------------------------------
  //column0: Photon Energy (MeV)
  //column1: Scattering Coherent (cm2/g)
  //column2: Scattering Incoherent (cm2/g)
  //column3: Photo-Electric Absorption (cm2/g)
  //column4: Pair in Nuclear Field (cm2/g)
  //column5: Production in electron field (cm2/g)
  //column6: Total Attenuation with coherent scatt. (cm2/g)
  //column7: Total Attenuation without coherent scatt. (cm2/g)
  //----------------------------------------------------------
 
  FILE * infile = NULL;
  infile = fopen("cross_sectionH2O_2.txt", "r");
  if (infile == NULL) 
    {
      fprintf(stderr, "Can't open Cross-Section file. File must be on current directory. \n");
    }
  InfoFile infoCrossSection;
  InfoFile * infoCrossSectionArray = new InfoFile[48];
  int line = 0;


  while(fscanf (infile, "%f %f %f %f %f %f %f %f", &infoCrossSection.c0, &infoCrossSection.c1, &infoCrossSection.c2, &infoCrossSection.c3, &infoCrossSection.c4, &infoCrossSection.c5, &infoCrossSection.c6, &infoCrossSection.c7) == 8)
    {
 
      infoCrossSectionArray[line] = infoCrossSection;
      line++;
    }

  fclose(infile);
  return infoCrossSectionArray;
}

float plane(float p0, float z0, float w)
{
  float s;
  if(((z0 < p0) && (w > 0)) || ((z0 > p0) && (w < 0)))
    {
      s = 0 > ((p0-z0)/w) ? 0 : ((p0-z0)/w); //the distance has to be positive
    }
  else
    {
      s = INFINITY;
    }

  return s;
} 


PhotonCoord detector(float half_plates_distance, float axialLength, float transaxialLength, FinalCoordDetectAxis finalCoordDetect)
{

  PhotonCoord photonCoord;

  //Photon projection to the detector (assuming 0deg as rotation angle)
  
  if (finalCoordDetect.w != 0)
    {
      float detect0 = plane(-half_plates_distance, finalCoordDetect.z, finalCoordDetect.w); //distance to DH
      float detect1 = plane(half_plates_distance, finalCoordDetect.z, finalCoordDetect.w);

     
      if (detect0 != INFINITY || detect1 != INFINITY)
	{
	  float s = detect0 < detect1 ? detect0 : detect1;
	  float x = finalCoordDetect.x + finalCoordDetect.u * s;
	  float y = finalCoordDetect.y + finalCoordDetect.v * s;
	  float z = finalCoordDetect.z + finalCoordDetect.w * s;
	  
	  if((x >= -axialLength/2) && (x <= axialLength/2) && (y >= -transaxialLength/2) && (y <= transaxialLength/2))
	    {
	      photonCoord.x = x;
	      photonCoord.y = y;
	      photonCoord.z = z;
	    } 
	  else
	    {
	      photonCoord.x = 0;
	      photonCoord.y = 0;
	      photonCoord.z = 0;
	    }	  
	}
      else
	{
	  photonCoord.x = 0;
	  photonCoord.y = 0;
	  photonCoord.z = 0;
	}
    }
  else
    {
      photonCoord.x = 0;
      photonCoord.y = 0;
      photonCoord.z = 0;
    }

  
  return photonCoord;

}


std::vector<float> crossSection(float E, InfoFile * crossSecInfo)
{
  float E_Mev = E *0.001;      //(MeV)
  
  //INTERPOLATE
  int i=0;
  float compt1 = 0;
  float photo1 = 1;
  float e1 = 0;
  while (crossSecInfo[i].c0 < E_Mev)
    {
      e1 = crossSecInfo[i].c0;
      compt1 = crossSecInfo[i].c2;
      photo1 = crossSecInfo[i].c3;
      i++;
    }

  //at the end of the while loop i= i(used for compt1) +1

  float e2 = crossSecInfo[i].c0;
  float compt2 = crossSecInfo[i].c2;
  float photo2 = crossSecInfo[i].c3;
 
  
  float interpolCompt = compt1 + ((compt2-compt1)/(e2 - e1)) * (E_Mev - e1);  //(cm2/g)
  float interpolPhoto = photo1 + ((photo2-photo1)/(e2 - e1)) * (E_Mev - e1);  //(cm2/g)

  
  std::vector<float> crossSec(2);
  crossSec[0] = interpolCompt;
  crossSec[1] = interpolPhoto;

  return crossSec; 
}



std::vector<float> compton(float E0, boost::mt19937 *generator)
{

  //random number generator
  boost::uniform_real<> range(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > next(*generator, range);


  float mec2 = 511.0; //keV

  float Eps0 = mec2/(mec2 + 2.0 * E0);
  float Eps0sq = pow(Eps0,2);

  float alpha1 = log(1/Eps0);
  float alpha2 = (1.0 - Eps0sq)/2.0;


  float Eps;
  float cosTHETA;
  while(1)
    {
      float r1 = next();                     
      float r2 = next();                  
      float r3 = next();           
      

      if (r1 < (alpha1/(alpha1+alpha2)))
	{
	  Eps = exp(-r2 * alpha1);
	}
      else
	{
	  Eps = sqrt(Eps0sq + ((1.0 - Eps0sq) * r2));
	}
      
      float t = mec2 * (1.0 - Eps)/(E0 * Eps);  // =(1-cos(theta))

      float sin2THETA = t * (2.0 - t);
      cosTHETA = 1.0 - t;

      float g = (1.0 - (Eps / (1 + pow(Eps,2))) * sin2THETA); //rejection function

      if (g < r3){continue;} else {break;}
      
    }

  float THETA = acos(cosTHETA);  //(Rad)
  float E1 = Eps * E0;  //(keV)

  std::vector<float> result(2); 
  result[0] = E1;
  result[1] = THETA;

  return result;

}



ScatterCoord scatter(Coord coord, InfoFile * crossSecInfo, float cutOff, boost::mt19937 *generator, int nphoton, float xDim, float yDim, float axialDim, ImageType::SpacingType imageSpacing, ImageType::SizeType size, float Min_Ethresh)
{

  //random number generator
  boost::uniform_real<> range(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > next(*generator, range);


  float x0 = coord.x;
  float y0 = coord.y;
  float z0 = coord.z;
  float u0 = coord.u;
  float v0 = coord.v;
  float w0 = coord.w;
  float E  = coord.E;
  float theta0 = coord.THETA;
  float phi0 = coord.PHI;



  int nCompton = 0;
  

  ScatterCoord finalCoord;

  float densH2O = 1; //(g/cm3) water density
  std::vector<float> sigma(2);
  sigma = crossSection(E, crossSecInfo); //(cm2/g)  //returns crossSec[compt, photo]
  float sigma_compt = sigma[0];
  float sigma_photo = sigma[1];
  float niu_compt = densH2O * sigma_compt;          //(cm-1)
  niu_compt = niu_compt * (1/10.);            //(mm-1)
  float niu_photo = densH2O * sigma_photo;          //(cm-1)
  niu_photo = niu_photo * (1/10.);            //(mm-1)
  float niu_total = niu_compt + niu_photo;          //(mm-1)
  float P_compt = niu_compt/niu_total;
  float P_photo = niu_photo/niu_total;
  


  //Distance of interaction
  float s_int = -log(1-next())*(1/niu_total);  //(mm) 

  float x0_int = x0 + u0 * s_int;
  float y0_int = y0 + v0 * s_int;
  float z0_int = z0 + w0 * s_int;

  
  ImageType::IndexType pixIndex; 
  ImageType::PixelType pixValue;
  pixIndex[0] = (round(x0_int) + (xDim/2))/imageSpacing[0]; 
  pixIndex[1] = (round(y0_int) + (yDim/2))/imageSpacing[1]; 
  pixIndex[2] = (round(z0_int) + (axialDim/2))/imageSpacing[2]; 
  if (pixIndex[0]> (int)size[0] || pixIndex[1]> (int)size[1] || pixIndex[2]> (int)size[2] || pixIndex[0]< 0 || pixIndex[1]< 0 || pixIndex[2]< 0){pixValue = 0;}
  else {pixValue = reader->GetOutput()->GetPixel(pixIndex);}


  while ((E > cutOff) && (pixValue > 0))     
    {                                          //YES, an interaction will take place

      if (E < Min_Ethresh)          //accelerates execution
	{
	  finalCoord.E = 0;
	  return finalCoord;
	}  

      float Prob = next();     
      
      if (Prob <= P_photo) 
	{
	  x0 = x0 + u0 * s_int;
	  y0 = y0 + v0 * s_int;
	  z0 = z0 + w0 * s_int;


	  finalCoord.x = x0;
	  finalCoord.y = y0;
	  finalCoord.z = z0;
	  finalCoord.u = 0;
	  finalCoord.v = 0;
	  finalCoord.w = 0;
	  finalCoord.E = 0;
	  finalCoord.nCompton = nCompton;
	  
	  return finalCoord;
	}
      else
	{
	  nCompton = nCompton + 1;

	  std::vector<float> new_param(2); 
	  new_param = compton(E, generator); //returns an array: [E, THETA]

	  E = new_param[0];              //(keV)
	  float THETA_compt = new_param[1];    //(rad)
	  float PHI = 2 * M_PI * next(); //(rad) 


	  //new coordinates after scatter
	  x0 = x0 + u0 * s_int;
	  y0 = y0 + v0 * s_int;
	  z0 = z0 + w0 * s_int;

	  u0 = u0 * cos(THETA_compt) + sin(THETA_compt)*(w0*cos(PHI)*cos(phi0) - sin(PHI)*sin(phi0));  //x-axis direction cosine of the particle
	  v0 = v0 * cos(THETA_compt) + sin(THETA_compt)*(w0*cos(PHI)*sin(phi0) + sin(PHI)*cos(phi0));  //y-axis direction cosine of the particle
	  w0 = w0 * cos(THETA_compt) - sin(THETA_compt)*sin(theta0)*cos(PHI);                          //z-axis direction cosine of the particle
        

	  sigma = crossSection(E, crossSecInfo); //(cm2/g)
	  sigma_compt = sigma[0];
	  sigma_photo = sigma[1];
	  niu_compt = densH2O * sigma_compt;          //(cm-1)
	  niu_compt = niu_compt * (1/10.);            //(mm-1)
	  niu_photo = densH2O * sigma_photo;          //(cm-1)
	  niu_photo = niu_photo * (1/10.);            //(mm-1)
	  niu_total = niu_compt + niu_photo;          //(mm-1)
	  P_compt = niu_compt/niu_total;
	  P_photo = niu_photo/niu_total;
	  

	  //Distance of interaction
	  s_int = -log(1-next())*(1/niu_total);   //(mm) 


	  //Interaction?? Next position after interaction has to be inside object => pixValue = 1
	  // if not, x y z keep the same value (no interaction) and are projected to the DHs.  
	  x0_int = x0 + u0 * s_int;
	  y0_int = y0 + v0 * s_int;
	  z0_int = z0 + w0 * s_int;

	  pixIndex[0] = (round(x0_int) + (xDim/2))/imageSpacing[0]; 
	  pixIndex[1] = (round(y0_int) + (yDim/2))/imageSpacing[1]; 
	  pixIndex[2] = (round(z0_int) + (axialDim/2))/imageSpacing[2]; 
	  if (pixIndex[0]> (int)size[0] || pixIndex[1]> (int)size[1] || pixIndex[2]> (int)size[2] || pixIndex[0]< 0 || pixIndex[1]< 0 || pixIndex[2]< 0){pixValue = 0;}
	  else {pixValue = reader->GetOutput()->GetPixel(pixIndex);}
	}
    }


  finalCoord.x = x0;
  finalCoord.y = y0;
  finalCoord.z = z0;
  finalCoord.u = u0;
  finalCoord.v = v0;
  finalCoord.w = w0;
  finalCoord.E = E;
  finalCoord.nCompton = nCompton;


  return finalCoord;
    
}

