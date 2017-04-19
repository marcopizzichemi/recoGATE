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

// #include "../lmFormats/formats.hpp"
// #include "../lmFormats/EventStore.hpp"
// #include "lmrec_common.hpp"
// #include "pem_geometry.hpp"
// #include "interfile.hpp"
// #include "Metz.hpp"

//Multithread
#include <pthread.h>


int main(int argc, char *argv[])
{
	using namespace std;

  int XDim = atoi(argv[1]);
  int YDim = atoi(argv[2]);
  int ZDim = atoi(argv[3]);

  float* recImage = new float[XDim*YDim*ZDim];
  for (int i=0; i<XDim*YDim*ZDim; i++)
		recImage[i] = 1.0;
  
	FILE *outputFile = fopen("unit.lm", "wb");
	fwrite(recImage, sizeof(float), XDim*YDim*ZDim, outputFile);
	fclose(outputFile);
  return 0;
}
