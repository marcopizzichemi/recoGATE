#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "math.h"
#include <boost/lexical_cast.hpp>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <map>
#include <sys/timeb.h>
#include <getopt.h>

#include "lmrec_common.hpp"
#include "pem_geometry.hpp"



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
			fprintf(stderr, "Usage: %s [--threads nThreads ][ --rays nRays] <plates distance> <pixel length> <angles file> <input file> <output prefix>\n", argv[0]);
			return 1;
		}
		
	}

	if ((argc - optind) < 5) {
			fprintf(stderr, "Usage: %s [--threads nThreads ][ --rays nRays] <plates distance> <pixel length> <angles file> <input file> <output prefix>\n", argv[0]);
			
		return 1;
	}
	
	fprintf(stderr, "Using %d threads\n", nThreads);

	float platesDistance = boost::lexical_cast<float>(argv[optind+0]);
	float pixelLength = boost::lexical_cast<float>(argv[optind+1]);
	char * angleFileName = argv[optind+2];
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

	FILE *angleFile = fopen(angleFileName, "r");
	assert(angleFile != NULL);
	while(1) {
		float angle;
		int r = fscanf(angleFile, "%f\n", &angle);
		if(r != 1) break;
		ip.angles.insert(angle*M_PI/180);
	}
	fclose(angleFile);
	
	std::set<float>	 uniqueAngles;
	for(std::multiset<float>::iterator i = ip.angles.begin(); i != ip.angles.end(); i++)
		uniqueAngles.insert(*i);

	fprintf(stderr, "Image size is %d x %d x %d; pixel size is %f\n", ip.xDim, ip.yDim, ip.zDim, ip.pixelLength);
	fprintf(stderr, "Angles: ");
	for(std::set<float>::iterator i = uniqueAngles.begin(); i != uniqueAngles.end(); i++) {
		fprintf(stderr, "%1.1fx%u ", (*i) * 180/M_PI, unsigned(ip.angles.count(*i)));
	}
	fprintf(stderr, "\n");


	char fnBuffer[1024];
    
	//begin to calculate totalImage

	float PIXELLENGTH = ip.pixelLength;
	int XDim = ip.xDim;
	int YDim = ip.yDim;
	int ZDim = ip.zDim;

	float* recImage = new float[XDim*YDim*ZDim];
	for (int i=0; i<XDim*YDim*ZDim; i++)
		recImage[i] = 1.0;

	float** diffImage = (float**)malloc(sizeof(float*)*nThreads);    
	for (int i=0; i<nThreads; i++)
	{
		diffImage[i] = new float[XDim*YDim*ZDim];
		for (int j=0; j<XDim*YDim*ZDim; j++)
			diffImage[i][j] = 0;
		
	}
    
	float *randImage = new float[ip.xDim*ip.yDim*ip.zDim];
	for (int i=0; i<ip.xDim*ip.yDim*ip.zDim; i++)
		randImage[i] = 0.0;
	float rand_cut = 0.0;
		
	double eventSum = 0;
	FILE *lmFile = NULL;
	if(strcmp(inputFileName, "-") == 0)
		lmFile = stdin;
	else 
		lmFile = fopen(inputFileName, "rb");

	pthread_t threads[nThreads];  
	struct rec_kernel_parameter parameter[nThreads];
	for(int i = 0; i < nThreads; i++) {
		parameter[i].nRays = nRays;
		parameter[i].generator = new boost::mt19937();
		parameter[i].generator->seed(time_t(i));
	}

	int blockSize = 100000*nThreads;
	DKFZFormat *buffer1 = new DKFZFormat [blockSize];
	DKFZFormat *buffer2 = new DKFZFormat [blockSize];
	while(true) {        
		int n = fread((void*)buffer1, sizeof(DKFZFormat), blockSize, lmFile);
		if (n <= 0) break;

		for (int i = 0; i < n; i++)
			eventSum += buffer1[i].weight;

		memcpy(buffer2, buffer1, n * sizeof(DKFZFormat));
		for(std::set<float>::iterator i = uniqueAngles.begin(); i != uniqueAngles.end(); i++) {
			float theta = *i;
			int count = ip.angles.count(theta);
			
			

			for(int j = 0; j < n; j++) {
				buffer2[j].angle = theta;
				buffer2[j].weight *= buffer1[j].weight * count;
			}

			int eventsPerThread = n / nThreads;
			if (eventsPerThread < 1) eventsPerThread = 1;

			for(int i=0; i<nThreads; i++) {
				int start = i * eventsPerThread;
				int end = (i+1) * eventsPerThread;
				if (end > n) end = n;
				if(end < start) end = start;
	
				parameter[i].ip = &ip;
				parameter[i].listData = buffer2;
				parameter[i].beginNum = start;
				parameter[i].totalNum = (end - start);
				parameter[i].diffImage = diffImage[i];
				parameter[i].recImage = recImage;
				parameter[i].randImage = randImage;
				//parameter[i].normImage1 = normImage1;
				//parameter[i].normImage2 = normImage2;

				pthread_create(&threads[i], NULL, rec_kernel, (void*)&parameter[i]);	
			}

			for (int i=0; i<nThreads; i++)
				pthread_join(threads[i], NULL);	
		}
	}
	delete [] buffer2;
	delete [] buffer1;
	fclose(lmFile);



	// Update
	for (int i=0; i<XDim*YDim*ZDim; i++)
	{
		float diffThread = 0;
		for (int j=0; j<nThreads; j++)
		{
			diffThread += diffImage[j][i];
		}
		// recImage[i] = diffThread (sum of diffImage for each thread)
		recImage[i] = recImage[i] * diffThread;
	}
	

	if(eventSum == 0) {
		// If we had not events, we'll just force the output image to be full of ones
		for (int i=0; i<XDim*YDim*ZDim; i++)
                recImage[i] = 1.0;
	} 
	else {
		// Normalize
		double total = 0.0;
		for (int i=0; i<XDim*YDim*ZDim; i++)
			total += recImage[i];
		if (total == 0) total = 1.0; // Avoid division by zero
		for (int i=0; i<XDim*YDim*ZDim; i++)
			recImage[i] = recImage[i] * eventSum / total;
		
	}
		
	sprintf(fnBuffer, "%s_total", outputPrefix);
	FILE *outputFile = fopen(fnBuffer, "wb");
	fwrite(recImage, sizeof(float), XDim*YDim*ZDim, outputFile);
	fclose(outputFile);
	return 0;
}

