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

	int nThreads = sysconf(_SC_NPROCESSORS_ONLN);
	int nRays = 0;
	static struct option longOptions[] = {
		{ "threads", required_argument, 0, 0 },
		{ "rays", required_argument, 0, 0 },
		{ NULL, 0, 0, 0 }
	};
	
	while(1) {
		
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:d:n:", longOptions, &optionIndex);
		
		if (c == -1) {
			break;
		}

		if (c == 0 && optionIndex == 0) {
			nThreads = boost::lexical_cast<int>((char *)optarg);
		}
		else if (c == 0 && optionIndex == 1) {
			nRays = boost::lexical_cast<int>((char *)optarg);
		}
		else {
			fprintf(stderr, "Usage: %s [--threads nThreads ][ --rays nRays] <plates distance> <pixel length> <norm prefix> <input file> <output file>\n", argv[0]);
			return 1;
		}
		
	}

	if ((argc - optind) < 5) {
			fprintf(stderr, "Usage: %s [--threads nThreads ][ --rays nRays] <plates distance> <pixel length> <norm prefix> <input file> <output file>\n", argv[0]);
			
		return 1;
	}

	fprintf(stderr, "Using %d threads\n", nThreads);

	float platesDistance = boost::lexical_cast<float>(argv[optind+0]);
	float pixelLength = boost::lexical_cast<float>(argv[optind+1]);
	char * normPrefix = argv[optind+2];
	char * inputFileName = argv[optind+3];
	char * outputPrefix = argv[optind+4];

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

	float *randImage = new float[ip.xDim*ip.yDim*ip.zDim];
	for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
		randImage[i] = 0.0;
	float rand_cut = 0.0;
	
	set<float> angleSet;
	double eventWeightSum = 0;
	long nEvents = 0;
	
	int nSubsets = 1;
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
			
			if (buffer[i].weight > 0)
				continue;
				
			buffer[i].weight *= -1;
		
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

	timeb t0;
	ftime(&t0);

	for(int subset = 0; subset < nSubsets; subset++) {
		eventStore[subset]->rewind();			
		while(true) {
			int nEvents = eventStore[subset]->readEvents(buffer, blockSize);
			if(nEvents <= 0) break;
			
			int eventsPerThread = int(nEvents / nThreads);				
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

			//normalize
			double total = 0.0;
			for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
				total += recImage[i];
			if(total == 0) total = 1.0; // Avoid division by zero
			for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
				recImage[i] = recImage[i] * eventWeightSum / total;
			
			for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
				for (int j=0; j<nThreads; j++)
				{
					diffImage[j][i] = 0.0;
				}

		}
	}
		//fprintf(stderr, "Finished angle %4.1f\r", theta * 180 / M_PI); fflush(stdout);
	
	//filter
	gaussian_filter_2d(recImage, ip.pixelLength, &ip);
	
	// Final normalization
	double total = 0.0;
	for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
		total += recImage[i];
	if (total == 0) total = 1.0; // AVoid division by zero
	for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
		recImage[i] = recImage[i] * eventWeightSum / total;

	char filename[4096];
	sprintf(filename, "%s", outputPrefix);

	FILE *outputFile = fopen(filename, "wb");
	fwrite(recImage, sizeof(float), ip.xDim*ip.yDim*ip.zDim, outputFile);
	fclose(outputFile);

	timeb t1;
	ftime(&t1);
	float dt = t1.time - t0.time + 0.001*(t1.millitm - t0.millitm);
	fprintf(stderr, "Finished: %ld events in %f seconds; %f events/second\n", 
		   nEvents, dt,
		   nEvents/dt);

	
//all done    
	return 0;
}

