#ifndef PSAFILTERGRIDSEARCH_HH
#define PSAFILTERGRIDSEARCH_HH

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <string.h>

#include "Global.hh"
#include "SignalBasis.hh"

using namespace std;

class PSAFilterGridSearch{
public:
  PSAFilterGridSearch(string basisfile);
  ~PSAFilterGridSearch();

  void ProcessOneEvent(pointExp &PE);
  void PrepareEvent(pointExp &PE, bool doReset);
  void InitialSolution(pointExp &PE, bool keepPrevious);
  int ProcessTheEvent(pointExp &PE);
  int SearchFullGrid(pointExp &PE, int netChSeg, float netChEner);
  int SearchAdaptive(pointExp &PE, int netChSeg, bool bCoarseOnly);

  int MakeSearchWave(pointExp &PE);
  int MakeLocalMask(pointExp &PE);

  
  void GetPtPos(int netChSeg, int id, float *pos){
    sbasis->GetPtPos(netChSeg, id, pos);
  }
  
private:

  SignalBasis *sbasis;

  float Chi2InnerLoop(const float *pReal, const float *pBase, float bScale);

};

#endif // #ifndef PSAFILTERGRIDSEARCH_HH
