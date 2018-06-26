#include <vector>

//struct of energy depositions events.
// where (crystal id and position in the crystal) was the energy deposited, how much and when
// local coordinates: remember that in GATE simulation,
struct enDep
{
  Float_t  localPosX;
  Float_t  localPosY;
  Float_t  localPosZ;
  Int_t    parentID;
  Int_t    trackID;
  Int_t    globalCryID;
  Int_t    eventID;
  Int_t    gantryID;
  Int_t    rsectorID;
  Int_t    moduleID;
  Int_t    submoduleID;
  Int_t    crystalID;
  Int_t    layerID;
  Int_t    sourceID;
  Double_t time;
  Float_t  edep;
};

struct avgCryEnergyDep
{
  Int_t eventID;
  Int_t parentID;
  Int_t trackID;
  int   globalCryID;
  Int_t    gantryID;
  Int_t    rsectorID;
  Int_t    moduleID;
  Int_t    submoduleID;
  Int_t    crystalID;
  Int_t    layerID;
  Int_t    sourceID;
  float x;
  float y;
  float z;
  float sx;
  float sy;
  float sz;
  float energy;
  float time;
};

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
