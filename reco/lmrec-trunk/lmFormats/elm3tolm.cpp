#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <set>

#include "formats.hpp"

static const float eMin = 400;
static const float eMax = 650;
static const float dMax = 4E-9;
static const int   maxHits = 4;

static const rejectionWindow = 20E-9;

int main(int argc, char * argv[]) 
{
	char * outputFileName = NULL;
	char * randomFileName = NULL;

	bool noDOI = false;
	bool lmSTIR = false;
	
	static struct option longOptions[] = {
			{ "energy-low", required_argument, 0, 0 }, 
			{ "energy-high", required_argument, 0, 0 }, 
			{ "time-window", required_argument, 0, 0 },			
			{ "random-file", required_argument, 0, 0 },
			{ "no-doi", no_argument, 0, 0 },
			{ "stir", no_argument, 0, 0 },
			{ "dkfz", no_argumentm, 0, 0 },
			//{ "dw-correction", no_argument, 0, 0},
			{ NULL, 0, 0, 0 }
		};


	while(1) {

		int optionIndex = 0;
		int c = getopt_long(argc, argv, "o:", longOptions, &optionIndex);
		
		if (c == -1) {
			break;
		}
		
	
		if (c == 'o') {
			outputFileName = optarg;
		}
		else if(c == 0 && optionIndex == 0) {
			eMin = boost::lexical_cast<float>(optarg);
		}
		else if (c == 0 && optionIndex == 1) {
			eMax = boost::lexical_cast<float>(optarg);
		}
		else if (c == 0 && optionIndex == 2) {
			dMax = boost::lexical_cast<float>(optarg) * 1E-9;
		}
		else if (c == 0 && optionIndex == 3) {
			randomFileName = optarg;
		}
		else if (c == 0 && optionIndex == 4) {
			noDOI = true;
		}
		else if(c == 0 && optionIndex = 5) {
			lmSTIR = true;
		}
		else if(c == 0 && optionIndex = 6) {
			lmStir = false;
		}
		else {
			std::cout << "Usage: " << argv[0]
				  << "[ -o <outfile.lm> "
				  << "[ --energy-low <energy in keV> ] " 
				  << "[ --energy-high <energy in keV> ] "
				  << "[ --time-window <time window in ns> ] "
				  << "[ --randoms-file <randoms.lm> ] "
				  << "[ --stir || --dkfz ] "
		//		  << "[ --dw-correction ]" 
				  << std::endl;
			return 1;
		}		
	}

	if ((argc - optind) < 1) {
		fprintf(stderr, "Usage: %s [ options ] <file1> <file2> ... <fileN>\n", argv[0]);
		return 1;
	}



	double t0 = -1;
	FILE *outputFile = NULL;
	if (outputFileName != NULL) {
		if (strcmp(outputFileName, "-") == 0)
			outputFile = stdout;
		else
			outputFile = fopen(outputFileName, "wb");
					
		if (outputFile == NULL) {
			fprintf(stderr, "Could not open %s for writing\n", outputFileName);
			return 1;
		}
	}

	FILE * randomFile = NULL;
	if (randomFileName != NULL) {
		if (strcmp(randomFileName, "-") == 0)
			randomFile = stdout;
		else
			randomFile = fopen(randomFileName, "wb");
			
		if (randomFile == NULL) {
			fprintf(stderr, "Could not open %s for writing\n", randomFileName);
			return 1;
		}
	}

	long long nRead = 0;
	long long nPrompts = 0;
	long long nRandoms = 0;
	long long nBad = 0;

	std::set<float> angleSet;

	for(int i = optind; i < argc; i++) {
		FILE * fIn = NULL;
		if(strcmp(argv[i], "-") == 0)
			fIn = stdin;
		else
			fIn = fopen(argv[i], "rb");
				
		if (fIn == NULL) {
			fprintf(stderr, "File %s does not exist\n", argv[i]);
			continue;
		}
		
		u_int8_t h;
		while(fread((void*)&h, sizeof(h), 1, fIn) == 1) {
			ELM3::event_t event;
			ELM3::marker_t marker;
			char buffer[128];
		
			switch(h) {
				case(0)
			}
		
		
			if(nRead == 0) 
				fprintf(stderr, "Distance is %f\n", fe.d);

			nRead++;
			if(t0 < 0) t0 = fe.ts;

			// Place events at plates surface when DOI is not to be used
			if(noDOI) {
				fe.z1 = fe.z1 < 0 ? -fe.d/2.0 : fe.d/2.0;
				fe.z2 = fe.z2 < 0 ? -fe.d/2.0 : fe.d/2.0;
			}

 			float u = fe.x1 - fe.x2;
			float v = fe.y1 - fe.y2;
			float w = fe.z1 - fe.z2;
		

			if(fabs(w) < fe.d/2.0) {
				nBad++;
				continue;
			}

			angleSet.insert(fe.yozRot);

 			if(fe.e1 < eMin || fe.e1 > eMax ||
 			   fe.e2 < eMin || fe.e2 > eMax ||
 			   fabs(fe.dt) > dMax)
			   continue;

			if(fe.n1 > maxHits || fe.n2 > maxHits)
				continue;
			
			if (fe.z1 > 0)
				fe.dt *= -1;

			DKFZFormat lmEntry = { fe.x1, fe.y1, fe.z1, fe.x2, fe.y2, fe.z2 , fe.yozRot, 1 };
			
			if(fe.random != 0)
				lmEntry.weight *= -1;
			
			if(fe.random == 0) 
				if(outputFile != NULL)
					fwrite((void*) &lmEntry, sizeof(lmEntry), 1, outputFile);			
			
			if ((fe.random != 0) && (randomFile != NULL))
					fwrite((void*) &lmEntry, sizeof(lmEntry), 1, randomFile);
			

			if(fe.random == 0)
				nPrompts++;
			else
				nRandoms++;
		}

		fclose(fIn);
	}


	if (outputFile != NULL)
		fclose(outputFile);
	if (randomFile != NULL)
		fclose(randomFile);

	if(nBad > 0) fprintf(stderr, "Found %lld insane events\n", nBad);
	fprintf(stderr, "%lld prompts, %lld randoms\n", nPrompts, nRandoms);

	fprintf(stderr, "Angles: ");
	for(std::set<float>::iterator i = angleSet.begin(); i != angleSet.end(); i++) {
		fprintf(stderr, "%f ", (*i) * 180/M_PI);
	}
	fprintf(stderr, "\n");

	return 0;
	
}
