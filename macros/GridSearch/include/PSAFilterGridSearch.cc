#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <x86intrin.h>

#include "PSAFilterGridSearch.hh"

using namespace std;

PSAFilterGridSearch::PSAFilterGridSearch(string basisfile){

  // SignalBasis
  sbasis = new SignalBasis();
  sbasis->Read(basisfile);
  sbasis->MakeCoarseFineList();

}


void PSAFilterGridSearch::ProcessOneEvent(pointExp &PE){

  for(int rpt=0; rpt<4; rpt++){
    //////////////////////////////////////
    // 1  prepare event for this iteration
    //////////////////////////////////////

    PrepareEvent(PE, (rpt==0));

    if(!PE.isValid){
      return;
    }

    InitialSolution(PE, (rpt!=0) );  // first time: center of the fired segment; otherwise: position of previous iteration

    /////////////////////////////////////
    // 2  call the 1-hit grid-search algo
    /////////////////////////////////////

    int rv = ProcessTheEvent(PE);

  }

}


void PSAFilterGridSearch::PrepareEvent(pointExp &PE, bool doReset){

  if(doReset){
    // find hit segments
    int    numsegs  = PE.numNetCharges;
    int   *segOrder = PE.netChargeSegnum;
    float *eneOrder = PE.netChargeEnergy;
    float *posx     = PE.x;
    float *posy     = PE.y;
    float *posz     = PE.z;
    int   *nhits    = PE.nhits;
    
    if(numsegs<1 || numsegs>12){
      PE.isValid = false;
      return;
    }

    PE.netChSeg = -1;
    PE.isValid = true;
    if(numsegs>1){
      // order the net-charge segments according to the released energy (not done for GridSearch==SegCenter)
      for(int n1=numsegs-1; n1>=0; n1--){
	for(int n2=0; n2<n1; n2++){
	  if(eneOrder[n2]<eneOrder[n2+1]){
	    swap(segOrder[n2], segOrder[n2+1]);
	    swap(eneOrder[n2], eneOrder[n2+1]);
	    swap(posx[n2], posx[n2+1]);
	    swap(posy[n2], posy[n2+1]);
	    swap(posz[n2], posz[n2+1]);
	    swap(nhits[n2], nhits[n2+1]);
	  }
	}
      }
    }
    
  }else{
    PE.netChSeg = -1;
    PE.isValid = true;
  }
  
  // save tAmplitude into wAmplitude, where is can be dynamically modified by the search algorithms
  memcpy(PE.wAmp, PE.tAmp, sizeof(PE.wAmp));
}


// place everything at the center or keep what already found in the previous minimization
void PSAFilterGridSearch::InitialSolution(pointExp &PE, bool keepPrevious){
  for(int snum = 0; snum < PE.numNetCharges; snum++){
    int   netChSeg  = PE.netChargeSegnum[snum];
    float netChEner = PE.netChargeEnergy[snum];
    float scaleFact = netChEner;

    int bestPt = (keepPrevious && PE.isInitialized) ? PE.resPt[snum] : sbasis->GetPtlist(netChSeg)->at(0)[0];

    PE.resPt[snum]  = bestPt;
    PE.resFac[snum] = 1.f;

    PE.netChSeg = netChSeg;   // now working on this

    // subtract the result of this point from PF.wAmplitude
    pointPsa &bestPoint = sbasis->GetPts(netChSeg)->at(bestPt);
    PE.AddBaseTrace(bestPoint, -scaleFact);
  }
  
  PE.isInitialized = true;
}


int PSAFilterGridSearch::ProcessTheEvent(pointExp &PE){

  int chiDone = 0;
  
  if(!PE.isValid) return chiDone;

  for(int snum = 0; snum < PE.numNetCharges; snum++){
    int   netChSeg  = PE.netChargeSegnum[snum];
    float netChEner = PE.netChargeEnergy[snum];

    float scaleFact = netChEner;
    PE.baseScale   = scaleFact;
    PE.netChSeg    = netChSeg; // now working on this

    // add-back the initial/previous solution to PE.wAmplitude
    int bestPt = PE.resPt[snum];
    pointPsa &bestPoint = sbasis->GetPts(netChSeg)->at(bestPt);
    PE.AddBaseTrace(bestPoint, scaleFact);

    // prepare sAmplitude from the subset of active segments
    int nActive = MakeSearchWave(PE);

    PE.chi2min = float(1.e30);
    PE.bestPt = -1;
    
    bool onlyCoarse = gCoarseOnly;
    if(snum>1 && netChEner<lowSegEnergy) { // if mult>1 and the energy of this point is small, make only the coarse search
      onlyCoarse = true;
    }

    if(gUseAdaptive || onlyCoarse){
      int rv = SearchAdaptive(PE, netChSeg, onlyCoarse);
      chiDone += abs(rv);
    }else{
      int rv = SearchFullGrid(PE, netChSeg, netChEner);
      chiDone += abs(rv);
    }

    
    pointPsa &bestPoint1 = sbasis->GetPts(netChSeg)->at(PE.bestPt);

    // Subtract the result of this point from PF.wAmplitude, which is the the spectrum of residuals.
    PE.AddBaseTrace(bestPoint1, -scaleFact);
    PE.resPt[snum] = PE.bestPt;
    PE.resFac[snum] = 1.f;

  }  // loop over the net-charge segments

  return chiDone;
}


int PSAFilterGridSearch::SearchFullGrid(pointExp &PE, int netChSeg, float netChEner){
  char *lMask = PE.localMask;

  int   bestPt  = 0;
  float chi2min = PE.chi2min;
  int   chiDone = 0;

  // the base points of this segment
  vector<pointPsa> *pPtSeg = sbasis->GetPts(netChSeg);

  int iPts = 0;
  int nPts = pPtSeg->size();

  for(; iPts < nPts; iPts++){
    chiDone++;
    float baseScale = PE.baseScale;
    float chi2 = 0;
    for(int iseg=0; iseg<NCHAN; iseg++){
      if(lMask[iseg] != '0'){
	float realTrace[BSIZE];
	memcpy(realTrace, PE.sAmp[iseg], sizeof(realTrace));
	float *baseTrace = pPtSeg->at(iPts).Amp[iseg];
	chi2 += Chi2InnerLoop(realTrace, baseTrace, baseScale);
	if(chi2 > chi2min) break;
      }
    } // end loop over the segments
    if(chi2 < chi2min){
      bestPt  = iPts;
      chi2min = chi2;
    }    
  }// end loop over the base points iPts

  PE.bestPt  = bestPt;
  PE.chi2min = chi2min;

  return chiDone;
}


int PSAFilterGridSearch::SearchAdaptive(pointExp &PE, int netChSeg, bool bCoarseOnly){
  float baseScale = PE.baseScale; // scaling signals to data
  char *lMask     = PE.localMask; // selection of segments

  // the basis points of this segment
  vector<pointPsa> *pPtSeg = sbasis->GetPts(netChSeg);

  // the coarse-fine structure of this segment
  vector<vector<int>> *pPtlist = sbasis->GetPtlist(netChSeg);

  int bestPt = -1;
  float chi2min = PE.chi2min;

  /////////////////////////////
  // loop on the coarse grid //
  /////////////////////////////
  int   kbest   = -1;
  int   knpts   = pPtlist->size();
  float chi2    = 0;
  int   chiDone = 0;

  // prepare the two pointers for each valid segment and call Chi2InnerLoop
  for(int kloop = 0; kloop < knpts; kloop++){
    int kPts = pPtlist->at(kloop)[0];
    chiDone++;
    chi2 = 0;
    for(int iseg=0; iseg<NCHAN; iseg++){
      if(lMask[iseg] != '0'){
	//float *realTrace = PE.sAmp[iseg];
	float realTrace[BSIZE];
	memcpy(realTrace, PE.sAmp[iseg], sizeof(realTrace));
	float *baseTrace = pPtSeg->at(kPts).Amp[iseg];
	chi2 += Chi2InnerLoop(realTrace, baseTrace, baseScale);
	if(chi2 > chi2min) break;
      }
    } // end loop over the segments
    if(chi2 < chi2min){
      bestPt  = kPts;
      chi2min = chi2;
      kbest   = kloop;
    }
  } // end loop over the base points iPts

  if(bestPt < 0){
    return -chiDone;
  }

  PE.bestPt  = bestPt;
  PE.chi2min = chi2min;

  if(bCoarseOnly){
    return chiDone;
  }

  
  /////////////////////////////
  // loop on the fine   grid //
  /////////////////////////////
  bestPt = -1;
  int   jnpts   = pPtlist->at(kbest).size();;
  for(int jloop=1; jloop < jnpts; jloop++){
    int jPts = pPtlist->at(kbest)[jloop];
    chiDone++;
    chi2 = 0;
    for(int iseg=0; iseg<NCHAN; iseg++){
      if(lMask[iseg] != '0'){
	//float *realTrace = PE.sAmp[iseg];
	float realTrace[BSIZE];
	memcpy(realTrace, PE.sAmp[iseg], sizeof(realTrace));
	float *baseTrace = pPtSeg->at(jPts).Amp[iseg];
	chi2 += Chi2InnerLoop(realTrace, baseTrace, baseScale);
	if(chi2 > chi2min) break;
      }
    } // end loop over the segments
    if(chi2 < chi2min){
      bestPt  = jPts;
      chi2min = chi2;
    }
  } // end loop over the base points jPts   

  if(bestPt < 0){
    return -chiDone;
  }

  PE.bestPt  = bestPt;
  PE.chi2min = chi2min;

  return chiDone;

}


inline float PSAFilterGridSearch::Chi2InnerLoop(const float *pReal, const float *pBase, float bScale)
{
  const __m128 masks = _mm_set1_ps(-0.0f);
  const __m128 zeros = _mm_setzero_ps();

  __m128 loopFac    = _mm_set_ps(bScale, bScale, bScale, bScale);
  
  __m128* realTrace = (__m128*)pReal;
  __m128* baseTrace = (__m128*)pBase;

  __m128  wave = _mm_setzero_ps();
  __m128  diff = _mm_setzero_ps();
  __m128  chis = _mm_setzero_ps();

  for(int nn = 0; nn < LOOP_4SSE; nn++) {
    wave = _mm_mul_ps(loopFac, baseTrace[nn]);
    diff = _mm_sub_ps(realTrace[nn], wave);
#if   FIXED_METRIC == FIXED_METRIC_ABSVAL                               // pow(|d|,1  )
    chis = _mm_add_ps(chis, _mm_andnot_ps(masks, diff) );
#elif FIXED_METRIC == FIXED_METRIC_SQUARE                               // pow( d ,2  )
    chis = _mm_add_ps(chis, _mm_mul_ps(diff, diff) );
#elif FIXED_METRIC == FIXED_METRIC_1SQRT                                // pow(|d|,1/2)
    chis = _mm_add_ps(chis, _mm_sqrt_ps(_mm_andnot_ps(masks, diff)) );  
#elif FIXED_METRIC == FIXED_METRIC_2SQRT                                // pow(|d|,1/4)
    chis = _mm_add_ps(chis, _mm_sqrt_ps(_mm_sqrt_ps(_mm_andnot_ps(masks, diff )) ) );
#else
# error Inconsistency in the definition of the distance metric when using the SSE versions
#endif
  }

  float chi2 = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(chis, zeros), zeros)); // sum the 4 values and save in chi2
  return chi2;
}


int PSAFilterGridSearch::MakeSearchWave(pointExp &PE){

  // select the active segments
  int nActive = MakeLocalMask(PE);

  // calculate sAmplitude
  PE.MakeSearchWave(PE.localMask);

  return nActive;
}

int PSAFilterGridSearch::MakeLocalMask(pointExp &PE){
  // set mask of segments for the search loop
  char *localMask = PE.localMask;
  strcpy(localMask,sbasis->GetMask(PE.netChSeg));
  int sMult = PE.numNetCharges;
  // removing common segments with the other net charges
  for(int ii=0; ii<sMult; ii++){
    int ncseg = PE.netChargeSegnum[ii];
    if(ncseg!=PE.netChSeg){
      // net-charge segments are always removed
      localMask[ncseg] = '0';
    }
  }


  int ns = 0;
  for(int ii=0; ii<NCHAN; ii++){
    if(localMask[ii]!='0')
      ns++;
  }
  return ns;
}


