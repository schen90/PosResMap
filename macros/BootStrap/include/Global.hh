#ifndef GLOBAL_HH
#define GLOBAL_HH

#include <stdlib.h>
#include <iostream>
#include <x86intrin.h>

#define NOISE 1000000

#define FIXED_METRIC_NONE   0
#define FIXED_METRIC_ABSVAL 1
#define FIXED_METRIC_SQUARE 2
#define FIXED_METRIC_1SQRT  3
#define FIXED_METRIC_2SQRT  4

#define FIXED_METRIC FIXED_METRIC_ABSVAL

const bool  gUseAdaptive = true;  // default
const bool  gCoarseOnly = false;

const int   NSLIC = 6;            // 6 slices
const int   NSECT = 6;            // 6 sectors
const int   NSEGS = NSLIC*NSECT;  // 36 segments
const int   NCHAN = NSEGS + 1;    // 37 the 36 segments and CC

const int   BSIZE = 56;
const int   TSTEP = 10;           // time step of the experimental data and of the internal basis signals (ns)
const int   basis_tzero =   10;   // channel of T0

const int   LOOP_SAMP = 40;       // 400 ns
const int   LOOP_4SSE = LOOP_SAMP/4;

const float  minSegEnergy = 10.f; // energy threshold for a segment to be part of the decomposition
const float lowSegEnergy = 50.f;  // segment-energies smaller than this are decomposed in a simplified way

struct pointPsa{
  float   x, y, z;
  int     netChSeg;
  float   Amp[NCHAN][BSIZE];
};

struct pointExp{
  int     numNetCharges;
  int     netChargeSegnum[NSEGS]; // normally given in the order of energy release
  float   netChargeEnergy[NSEGS];
  float   tAmp[NCHAN][BSIZE]; // the total original experimental data
  float   wAmp[NCHAN][BSIZE]; // the experimental data to be used and modified by the grid search ==> residuals
  float   sAmp[NCHAN][BSIZE]; // experimental trace manipulated for being used by the "chi2" loops
  float   rAmp[NCHAN][BSIZE]; // the "fitted" trace
  
  int     nhits[NSEGS];
  float   x[NSEGS];
  float   y[NSEGS];
  float   z[NSEGS];
  int     resPt[NSEGS];  // the best solution point for the numNetCharges segments
  float   resFac[NSEGS]; // relative amplitude of the point
  
  bool    isValid;
  bool    isInitialized;
  float   baseScale;          // temporarily used to manage the chi2 mapping
  int     netChSeg;
  char    localMask[NCHAN];
  int     bestPt;
  float   chi2min;                    // not very meaningful outside the grid-search

  // Produces sAmplitude from wAmplitude.
  void MakeSearchWave(char *mask){
    memset(sAmp, 0, sizeof(sAmp));
    for(int iSegm = 0; iSegm < NCHAN; iSegm++){
      if(mask[iSegm] != '0'){
	float *sA = sAmp[iSegm]; // the wave to be composed
	float *wA = wAmp[iSegm]; // taking the data from this
	for(int kk=0; kk<BSIZE; kk++){
	  sA[kk] = wA[kk];
	}
      }
    }
  }

  // Modify wAmplitude (float) adding the base point passed in pPsa
  void AddBaseTrace(pointPsa &pPsa, const float fact) {
    int ptNetCharge = pPsa.netChSeg;
    for(int iSegm = 0; iSegm < NCHAN; iSegm++){ // segments and core
      float *realTrace = wAmp[iSegm]; // from the experimental data
      float *baseTrace;               // signal basis samples
      baseTrace = pPsa.Amp[iSegm];   // take them as they are

      for(int kk=0; kk<BSIZE; kk++) {
        (*realTrace++) += (*baseTrace++)*fact;
      }
    }
    
  }
  
};

#endif 
