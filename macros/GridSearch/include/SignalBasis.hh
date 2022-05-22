#ifndef SIGNALBASIS_HH
#define SIGNALBASIS_HH

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <string.h>

#include "Global.hh"

using namespace std;

class SignalBasis{
public:
  SignalBasis();
  ~SignalBasis();

  void Read(string basisfile);
  void MakeCoarseFineList();
  void MakeCoarseFineList(int seg);
  void MakeSegmentMap(int neighbours);

  vector<pointPsa> *GetPts(int seg){ return segPts[seg];}
  vector<vector<int>> *GetPtlist(int seg){ return Ptlist[seg];}
  char *GetMask(int seg){ return hmask[seg];}

  void GetPtPos(int netChSeg, int id, float *pos){
    if(id<0){
      pos[0] = averPt[netChSeg].x;
      pos[1] = averPt[netChSeg].y;
      pos[2] = averPt[netChSeg].z;
    }else{
      pos[0] = segPts[netChSeg]->at(id).x;
      pos[1] = segPts[netChSeg]->at(id).y;
      pos[2] = segPts[netChSeg]->at(id).z;
    }
  }

  int NearestPoint(float px, float py, float pz, int iseg, float& dist);
  int NearestCoarse(float px, float py, float pz, int iseg, float& dist);
  bool AssignIt(float px, float py, float pz, float cx, float cy, float cz);
  bool IsCoarse(int jj, int seg);
  
private:

  const int fstep = 2; // mm fine grid
  const int cstep = 6; // mm coarse grid
  const int cstepx = cstep; // mm coarse grid
  const int cstepy = cstep; // mm coarse grid
  const int cstepz = cstep; // mm coarse grid

  vector<pointPsa> *segPts[NSEGS];
  pointPsa averPt[NSEGS];
  vector<vector<int>> *Ptlist[NSEGS]; // coarse-fine list, 0 coarse grid
  char hmask[NSEGS][NCHAN+1];
};

#endif // #ifndef SIGNALBASIS_HH
