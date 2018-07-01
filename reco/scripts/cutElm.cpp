// compile with
// g++ -o cutElm cutElm.cpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <getopt.h>
#include <iostream>
#include <set>
#include <assert.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <vector>



struct EventFormat {
  double ts;                     // time of the event, in seconds. if i'm not mistaken is the absolute machine time in seconds
  u_int8_t random;               // random flag, seti it to 0
  float d;                       // distance between the heads, fixed to the head distance value...
  float yozRot;                  // rotation of the heads
  float x1;                      // x coordinate of the event, for detector 1
  float y1;                      // y coordinate of the event, for detector 1
  float z1;                      // z coordinate of the event, for detector 1
  float e1;                      // energy deposited in detector 1
  u_int8_t n1;                   // multiplicity of event, set it to 1
  float x2;                      // x coordinate of the event, for detector 2
  float y2;                      // y coordinate of the event, for detector 2
  float z2;                      // z coordinate of the event, for detector 2
  float e2;                      // energy deposited in detector 2
  u_int8_t n2;                   // multiplicity of event, set it to 1
  float dt;                      // delta time between event in detector 1 and event in detector 2
} __attribute__((__packed__));


int main(int argc, char * argv[])
{
  long long nEvents = 0;

  std::ofstream* CurrentOutput = new std::ofstream("cut_output.elm2", std::ios::binary);

  FILE * fIn = NULL;
  if(strcmp(argv[1], "-") == 0)
  {
    fIn = stdin;
  }

  else
  {
    fIn = fopen(argv[1], "rb");
  }


  if (fIn == NULL) {
    fprintf(stderr, "File %s does not exist\n", argv[1]);
    return 1;
  }

  long long int events = atoi(argv[2]);
  printf("Reading %s up to event %lld\n", argv[1],events);
  EventFormat fe;
  //read current input file
  while(fread((void*)&fe, sizeof(fe), 1, fIn) == 1 && nEvents < (events *1))
  {
    if(false) // no output
    {
      std::cout << fe.ts     << " ";
      std::cout << (int) fe.random << " ";
      std::cout << fe.d      << " ";
      std::cout << fe.yozRot << " ";
      std::cout << fe.x1     << " ";
      std::cout << fe.y1     << " ";
      std::cout << fe.z1     << " ";
      std::cout << fe.e1     << " ";
      std::cout << (int) fe.n1     << " ";
      std::cout << fe.x2     << " ";
      std::cout << fe.y2     << " ";
      std::cout << fe.z2     << " ";
      std::cout << fe.e2     << " ";
      std::cout << (int) fe.n2     << " ";
      std::cout << fe.dt         ;
      std::cout << std::endl;
    }
    CurrentOutput->write((char*)&fe,sizeof(fe));
    nEvents++;
  }
  CurrentOutput->close();
  std::cout << "done!" << std::endl;


  return 0;
}
