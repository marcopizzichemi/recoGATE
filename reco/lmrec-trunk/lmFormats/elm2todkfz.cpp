#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <set>

#include "formats.hpp"

static float eMin = 350;
static float eMax = 700;
static float dMax = 20E-9;
static int   maxHits = 2;
static float halfLifeTime = 0;
static double t0 = -1;


int main(int argc, char * argv[])
{
	char * outputFileName = NULL;
	char * randomFileName = NULL;
	char * scatterFileName = NULL;

	bool noDOI = false;
        bool directCorrection = false;

	static struct option longOptions[] = {
			{ "energy-low", required_argument, 0, 0 },
			{ "energy-high", required_argument, 0, 0 },
			{ "time-window", required_argument, 0, 0 },
			{ "randoms-file", required_argument, 0, 0 },
			{ "no-doi", no_argument, 0, 0 },
			{ "half-life", required_argument, 0, 0 },
			{ "initial-time", required_argument, 0, 0 },
			{ "scatter-total-events", required_argument, 0, 0 },
                        { "use-direct-correction", no_argument, 0, 0 },
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
		else if (c == 0 && optionIndex == 5) {
			halfLifeTime = boost::lexical_cast<float>(optarg) * 60;
		}
		else if (c == 0 && optionIndex == 6) {
			t0 = boost::lexical_cast<double>(optarg);
		}
                else if (c == 0 && optionIndex == 7) {
                        scatterFileName = optarg;
                }
                else if (c == 0 && optionIndex == 8) {
                        directCorrection = true;
                }
		else {
			std::cout << "Usage: " << argv[0]
				  << "[ -o <outfile.lm> ]"
				  << "[ --scatter-total-events <count.txt> ]"
				  << "[ --energy-low <energy in keV> ] "
				  << "[ --energy-high <energy in keV> ] "
				  << "[ --time-window <time window in ns> ] "
				  << "[ --randoms-file <randoms.lm> ] "
				  << "[ --half-life <isotope half life in minutes> ]"
				  << "[ --initial-time <intial time for decay correction (unix epoch)> ]"
                                  << "[ --use-direct-correction ]"
				  << std::endl;
			return 1;
		}
	}

	if ((argc - optind) < 1) {
		fprintf(stderr, "Usage: %s [ options ] <file1> <file2> ... <fileN>\n", argv[0]);
		return 1;
	}

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


	FILE * scatterFile = NULL;
	if (scatterFileName != NULL) {
		if (strcmp(scatterFileName, "-") == 0)
			scatterFile = stdout;
		else
			scatterFile = fopen(scatterFileName, "wb");

		if (scatterFile == NULL) {
			fprintf(stderr, "Could not open %s for writing\n", scatterFileName);
			return 1;
		}
	}

	long long nRead = 0;
	double nPrompts = 0;
	double nRandoms = 0;
	long long nBad = 0;

	double tMin = INFINITY;
	double tMax = -INFINITY;

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

		ELM2Format fe;
		while(fread((void*)&fe, sizeof(fe), 1, fIn) == 1) {
			if(nRead == 0)
				fprintf(stderr, "Diameter is %f\n", fe.d);

			//fprintf(stderr, "%f %f %f -- %f %f %f --- %f\n", fe.x1, fe.y1, fe.z1, fe.x2, fe.y2, fe.z2, fe.dt);

			nRead++;
			if(t0 == -1) t0 = fe.ts;
			tMin = tMin < fe.ts ? tMin : fe.ts;
			tMax = tMax > fe.ts ? tMax : fe.ts;

			// Place events at plates surface when DOI is not to be used
      // FIXME to be modified if we want to use this flag in cylindrical conf
			if(noDOI) {
				fe.z1 = fe.z1 < 0 ? -fe.d/2.0 : fe.d/2.0;
				fe.z2 = fe.z2 < 0 ? -fe.d/2.0 : fe.d/2.0;
			}

 			float u = fe.x1 - fe.x2;
			float v = fe.y1 - fe.y2;
			float w = fe.z1 - fe.z2;

      //the "insanity condition" now is just "one of the two hits is inside the FOV"
			if( ( fabs(sqrt(pow(fe.y1,2)+pow(fe.z1,2))) < fe.d/2.0 ) | ( fabs(sqrt(pow(fe.y2,2)+pow(fe.z2,2))) < fe.d/2.0 ) ) {
				nBad++;
				continue;
			}

			angleSet.insert(fe.yozRot);

 			if(fe.e1 < eMin || fe.e1 > eMax ||
 			   fe.e2 < eMin || fe.e2 > eMax)
			   continue;
			fe.dt = fabs(fe.dt);
			float weight = 1;
			if(fe.dt > 90E-9)
				continue;
			else if(fe.dt > 20E-9) {
				fe.random = 1;
				weight = dMax/(90E-9 - 20E-9);
			}
			else if(fe.dt > dMax)
				continue;

			if(fe.n1 > maxHits || fe.n2 > maxHits)
				continue;

			// if (fe.z1 > 0)
				// fe.dt *= -1;

			float decay = (halfLifeTime == 0) ? 1.0 : pow(2, -(fe.ts - t0)/halfLifeTime);
			weight *= 1/decay;

			DKFZFormat lmEntry = { fe.x1, fe.y1, fe.z1, fe.x2, fe.y2, fe.z2 , fe.yozRot, weight };
			if(fe.random != 0)
				lmEntry.weight *= -1;

			if (fe.random == 0 || directCorrection)
				if(outputFile != NULL)
					fwrite((void*) &lmEntry, sizeof(lmEntry), 1, outputFile);

			if ((fe.random != 0) && (randomFile != NULL))
					fwrite((void*) &lmEntry, sizeof(lmEntry), 1, randomFile);

			//printf("%d %f %f\n", fe.random, weight, fe.dt * 1E9);

			if(fe.random == 0)
				nPrompts += weight;
			else
				nRandoms += weight;
		}

		fclose(fIn);
	}

	fprintf(scatterFile, "%10.0lf", nPrompts-nRandoms);

	if (outputFile != NULL)
		fclose(outputFile);
	if (randomFile != NULL)
		fclose(randomFile);

	if(nBad > 0) fprintf(stderr, "Found %lld insane events\n", nBad);
	fprintf(stderr, "%10.0lf prompts, %10.0lf randoms\n", nPrompts, nRandoms);

	fprintf(stderr, "Angles: ");
	for(std::set<float>::iterator i = angleSet.begin(); i != angleSet.end(); i++) {
		fprintf(stderr, "%f ", (*i) * 180/M_PI);
	}
	fprintf(stderr, "\n");

	printf("tMin = %f, tMax = %f, acqTime = %f\n", tMin, tMax, tMax - tMin);

	if(halfLifeTime != 0 && tMin < t0) {
		fprintf(stderr, "WARNING: t0 for decay correction is %lf, but minium time found in data was %lf seconds earlier\n", t0, t0-tMin);
		fprintf(stderr, "WARNING: Maybe this should be re-run with --initial-time %lf\n", tMin);

	}
	return 0;

}
