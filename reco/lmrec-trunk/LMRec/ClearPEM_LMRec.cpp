#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "math.h"
#include <time.h>
#include <assert.h>
#include <string.h>
#include <map>
#include <set>
#include <sys/timeb.h>
#include <getopt.h>
#include <iostream>
#include <boost/lexical_cast.hpp>

#include "../lmFormats/formats.hpp"
#include "../lmFormats/EventStore.hpp"
#include "lmrec_common.hpp"
#include "pem_geometry.hpp"
#include "interfile.hpp"
#include "Metz.hpp"

//Multithread
#include <pthread.h>


int main(int argc, char *argv[])
{
	using namespace std;
	

	float pixelLength = 1;
	float platesDistance = -1;
	char * inputFileName = NULL;
	char * outputPrefix = NULL;
	char * normPrefix = NULL;
	int nSubsets = 4;
	int nIterations = 8;
	int saveIterations = 0;
	char * randFileName = NULL;
	float rand_cut = 0.0;

	bool interfile = true;
	
	float filter_fwhm = 0;
	int filterIterations = 2;
	float filter_power = 0;
	float partialEvents = 1;

	int nThreads = sysconf(_SC_NPROCESSORS_ONLN);
	int nRays = 0;

	static struct option longOptions[] = {
			{ "pixel-length", required_argument, 0, 0 },
			{ "crystal-distance", required_argument, 0, 0 },
			{ "iterations", required_argument, 0, 0 },
			{ "osem", required_argument, 0, 0 },		
			{ "input", required_argument, 0, 0 }, 
			{ "output", required_argument, 0, 0 },
			{ "norm-prefix", required_argument, 0, 0 },
			{ "filter-fwhm", required_argument, 0, 0 },
			{ "filter-iter", required_argument, 0, 0 },
			{ "filter-power", required_argument, 0, 0 },
			{ "no-interfile", no_argument, 0, 0 },
			{ "partial", required_argument, 0, 0 },
			{ "threads", required_argument, 0, 0 },
			{ "rays", required_argument, 0, 0 },
			{ "rand-image", required_argument, 0, 0 },
			{ "rand-cut", required_argument, 0, 0 },
			{ "save-iterations", required_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
		};
	while(1) {
		
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:d:n:", longOptions, &optionIndex);
		
		if (c == -1) {
			break;
		}
		
		if (c == 'i') 
			inputFileName = (char *)optarg;
		else if (c == 'o') 
			outputPrefix = (char *)optarg;
		else if (c == 'd')
			platesDistance = boost::lexical_cast<float>((char *)optarg);
		else if (c == 'n') 
			normPrefix = (char *)optarg;
		else if (c == 0 && optionIndex == 0)
			pixelLength = boost::lexical_cast<float>((char *)optarg);
		else if (c == 0 && optionIndex == 1)
			platesDistance = boost::lexical_cast<float>((char *)optarg);
		else if (c == 0 && optionIndex == 2)
			nIterations = boost::lexical_cast<int>((char *)optarg);
		else if (c == 0 && optionIndex == 3)
			nSubsets = boost::lexical_cast<int>((char *)optarg);
		else if (c == 0 && optionIndex == 4)
			inputFileName = (char *)optarg;
		else if (c == 0 && optionIndex == 5)
			outputPrefix = (char *)optarg;
		else if (c == 0 && optionIndex == 6)
			normPrefix = (char *)optarg;
		else if (c == 0 && optionIndex == 7) {
			filter_fwhm = boost::lexical_cast<float>((char *)optarg);
		}
		else if (c == 0 && optionIndex == 8)
			filterIterations = boost::lexical_cast<int>((char *)optarg);
		else if (c == 0 && optionIndex == 9) {
			filter_power = boost::lexical_cast<float>((char *)optarg);
		}
		else if (c == 0 && optionIndex == 10)
			interfile = false;
		else if (c == 0 && optionIndex == 11) {
			partialEvents = boost::lexical_cast<float>((char *)optarg)/100;
		}
		else if (c == 0 && optionIndex == 12) {
			nThreads = boost::lexical_cast<int>((char *)optarg);
		}
		else if (c == 0 && optionIndex == 13) {
			nRays = boost::lexical_cast<int>((char *)optarg);
		}
		else if (c == 0 && optionIndex == 14) {
			randFileName = (char *)optarg;
		}
		else if (c == 0 && optionIndex == 15) {
			rand_cut = boost::lexical_cast<float>((char *)optarg);
		}
		else if (c == 0 && optionIndex == 16) {
			saveIterations = boost::lexical_cast<int>((char *)optarg);
		}
			 
		else {
			std::cout	<< "Usage: " << argv[0]
					<< "[ -i | --input <input file> ] " 
					<< "[ -o | --output <output prefix> ] "
					<< "[ -d | --crystal-distance <crystal distance in mm> ] "
					<< "[ -n | --norm-prefix <norm files prefix> ] "
					<< "[ --pixel-length <voxel size in mm> ] "
					<< "[ --iterations <number of iterations> ] "
					<< "[ --osem <0|1> ] "
					<< "[ --filter-fwhm <0 to disable or FWHM in mm> ] "
					<< "[ --filter-power < Metz power, 0 for gaussian ] "
					<< "[ --filter-iter <filter each n iterations> ]"
				        << "[ --no-interfile ]"
					<< "[ --partial <event fraction %>"
					<< "[ --rand-image <random file> ]"
					<< "[ --rand-cut <random cut value, 0 to 0.5> ]"
					<< "[ --save-iterations <N> ] "
					<< endl;

			return 1;
		}
	}

	if (platesDistance < 0) {
		std::cerr << "Need to set plates distance!" << endl;
		return 1;
	}

	if (inputFileName == NULL) {
		std::cerr << "Need to set input file name!" << endl;
		return 1;
	}
		

	if (outputPrefix == NULL) {
		std::cerr << "Need to set output prefix!" << endl;
		return 1;
	}	

	if (normPrefix == NULL) {
		std::cerr << "Need to set normalization prefix!" << endl;
		return 1;		
	}

	fprintf(stderr, "Using %d threads\n", nThreads);

	float transaxialSize = platesDistance - 60 < PEM_yLength ? platesDistance - 60 : PEM_yLength;
	struct image_params ip = {
		transaxialSize, 
		PEM_xLength, 
		platesDistance, 
		pixelLength,
		multiset<float>(),
		int(2 + ceil(PEM_xLength / pixelLength / 2)*2), //Ensure dimensions are even
		int(2 + ceil(transaxialSize / pixelLength / 2)*2),
		int(2 + ceil(transaxialSize / pixelLength / 2)*2)	
	};

	fprintf(stderr, "Image size is %d x %d x %d; pixel size is %f\n", ip.xDim, ip.yDim, ip.zDim, ip.pixelLength);	


	if (filterIterations > 0 &&  filter_fwhm > 0) {
		fprintf(stderr, "Filtering with 3D Metz filter of FHWM = %f mm and mp = %f at every %d iterations\n", filter_fwhm, filter_power, filterIterations);
	}
	
	if(saveIterations == 0)
		saveIterations = nIterations;

	FILE *lmFile = NULL;
	if(strcmp(inputFileName, "-") == 0)
		lmFile = stdin;
	else
		lmFile = fopen(inputFileName, "rb");
	assert(lmFile != NULL);

	char fnBuffer[1024];
	//sprintf(fnBuffer, "%s_crystals", normPrefix);
	//FILE *normFile = fopen(fnBuffer, "rb");
	sprintf(fnBuffer, "%s_total", normPrefix);
	FILE *totalImageFile = fopen(fnBuffer, "rb");

	//assert(normFile != NULL);
	assert(totalImageFile != NULL);

	float thr = 1;
 
	//read normalization file   
	//float *normImage1 = new float[PEM_uDim*PEM_vDim];
	//float *normImage2 = new float[PEM_uDim*PEM_vDim];
	//assert(fread(normImage1, sizeof(float), PEM_uDim*PEM_vDim, normFile) == PEM_uDim*PEM_vDim);
	//assert(fread(normImage2, sizeof(float), PEM_uDim*PEM_vDim, normFile) == PEM_uDim*PEM_vDim);
	//fclose(normFile);
    
	//read total image
	float *totalImage = new float[ip.xDim*ip.yDim*ip.zDim];
	char dump;
	assert(fread(totalImage, sizeof(float), ip.xDim*ip.yDim*ip.zDim, totalImageFile) == ip.xDim*ip.yDim*ip.zDim);
	assert(fread(&dump, 1, 1, totalImageFile) <= 0);
	fclose(totalImageFile);

	//read random correction image
	double randTotal = 0;
	float *randImage;
	randImage = new float[ip.xDim*ip.yDim*ip.zDim];
//	float *oldRecImage;
	if (randFileName != NULL)
	{
		FILE *randImageFile = fopen(randFileName, "rb");
		fread(randImage, sizeof(float), ip.xDim*ip.yDim*ip.zDim, randImageFile);
		fclose(randImageFile);
		for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
		{
			randTotal += randImage[i];
		}
/* 		oldRecImage = new float[ip.xDim*ip.yDim*ip.zDim];
		for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
			oldRecImage[i] = 1.0; */
	}
	else
	{
		for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
			randImage[i] = 0.0;
	}

	set<float> angleSet;
	double eventWeightSum = 0;
	long nEvents = 0;
	
	if (nSubsets <= 1) nSubsets = 1;
	vector<EventStore *> eventStore(nSubsets, NULL);
	for(int i = 0; i < nSubsets; i++) {
		eventStore[i] = new EventStore();
	}
	
	int blockSize = 100000*nThreads;
	DKFZFormat *buffer = new DKFZFormat [blockSize];
	while(true) {        
		int n = fread((void*)buffer, sizeof(DKFZFormat), blockSize, lmFile);
		if (n <= 0) break;	
		
		for(int i = 0; i < n; i++) {
			int subset = i % nSubsets;
			eventStore[subset]->addEvents(buffer+i, 1);
			
			ip.angles.insert(buffer[i].angle);
			eventWeightSum += buffer[i].weight;
			nEvents++;
		}
		
	}	
	fclose(lmFile);
	      
	float* recImage = new float[ip.xDim*ip.yDim*ip.zDim];
	for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
		recImage[i] = 1.0;
	if ( (randFileName != NULL)&&(rand_cut > 0) )
	{
		for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
			recImage[i] = randImage[i] * (float)(nEvents) / randTotal;
	}

	float** diffImage = (float**)malloc(sizeof(float*)*nThreads); 
	for (int i=0; i<nThreads; i++)
	{
		diffImage[i] = new float[ip.xDim*ip.yDim*ip.zDim];
		for (int j=0; j<ip.xDim*ip.yDim*ip.zDim; j++)
			diffImage[i][j] = 0;
	}

	pthread_t threads[nThreads];
	struct rec_kernel_parameter parameter[nThreads];
	for(int i = 0; i < nThreads; i++) {
		parameter[i].nRays = nRays;
		parameter[i].generator = new boost::mt19937();
		parameter[i].generator->seed(time_t(i));
	}

	for (int iter=0; iter<nIterations; iter++)
	{
		timeb t0;
		ftime(&t0);

		for(int subset = 0; subset < nSubsets; subset++) {
			eventStore[subset]->rewind();			
			while(true) {
				int nEvents = eventStore[subset]->readEvents(buffer, blockSize);
				if(nEvents <= 0) break;
				
				int eventsPerThread = int(nEvents / nThreads * partialEvents);				
				if (eventsPerThread < 1) eventsPerThread = 1;
			
				for(int i=0; i<nThreads; i++) {
	
					int start = i * eventsPerThread;
					int end = (i+1) * eventsPerThread;
					if (end > nEvents) end = nEvents;
					if(end < start) end = start;
			
					parameter[i].ip = &ip;
					parameter[i].listData = buffer;
					parameter[i].beginNum = start;
					parameter[i].totalNum = (end - start);
					parameter[i].diffImage = diffImage[i];
					parameter[i].recImage = recImage;
					parameter[i].randImage = randImage;
					parameter[i].rand_cut = rand_cut;
					//parameter[i].normImage1 = normImage1;
					//parameter[i].normImage2 = normImage2;
	
					pthread_create(&threads[i], NULL, rec_kernel, (void*)&parameter[i]);	
				}
	
				for (int i=0; i<nThreads; i++)
					pthread_join(threads[i], NULL);
			}

			{
				// Update
				for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
				{
					float diffThread = 0;
					for (int j=0; j<nThreads; j++)
					{
						diffThread += diffImage[j][i];
					} 
					recImage[i] = recImage[i] * diffThread;
					if (totalImage[i] > thr)
						recImage[i] = recImage[i] / totalImage[i]; 
				}

				zero_out(&ip, recImage);

				double total = 0.0;
		                for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
                		        total += recImage[i];
		                if(total == 0) total == 1.0; // avoid division by zero
		                for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
		                        recImage[i] = recImage[i] * (eventWeightSum-randTotal) / total;

				
                                for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++) {
                                    for (int j=0; j<nThreads; j++) {
                                        diffImage[j][i] = 0.0;
                                    }
                                }

			}

			//fprintf(stderr, "Finished angle %4.1f\r", theta * 180 / M_PI); fflush(stdout);
		}
		
		//smooth
		if ((filterIterations > 0) && (iter % filterIterations == filterIterations - 1)) {

			if(filter_fwhm > 0 ) {
				metz_filter_3d(recImage, filter_fwhm / 2.35482, filter_power, &ip);
			}
		}			
		
		// Final normalization
		double total = 0.0;
		for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
			total += recImage[i];
		if(total == 0) total == 1.0; // avoid division by zero
		for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
			recImage[i] = recImage[i] * (eventWeightSum-randTotal) / total;
		double total2 = 0.0;
		for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
			total2 += recImage[i];
		fprintf(stderr, "eventSum = %lf, randTotal = %lf, total = %lf, total2 = %lf\n", eventWeightSum, randTotal, total, total2);
		if(iter % saveIterations == saveIterations - 1) {
			char filename[4096];
			sprintf(filename, "%s_%03d", outputPrefix, iter);

			if (!interfile) {
				FILE *outputFile = fopen(filename, "wb");
				fwrite(recImage, sizeof(float), ip.xDim*ip.yDim*ip.zDim, outputFile);
				fclose(outputFile);
			}
			else {
				writeInterFile(filename, &ip, recImage);
			}
		}
		
		timeb t1;
		ftime(&t1);
		float dt = t1.time - t0.time + 0.001*(t1.millitm - t0.millitm);
		fprintf(stderr, "Finished iteration %02d: %ld events in %f seconds; %f events/second\n", 
		       iter, 
		       nEvents, dt,
		       nEvents/dt);
	}
	
//all done    
	return 0;
}

