#include <stdlib.h>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TRandom.h"
#include "TVirtualFFT.h"
#include "time.h"

using namespace std;

const int NSLIC = 6;
const int NSECT = 6;
const int NSEGS = NSLIC*NSECT;
const int NCHAN = NSEGS + 1;
const int NSIGS = 121; // npoints in 5ns step
const int NZERO = 10;
const int NSIGS2 = 80; // npoints in 5ns step

const int fstep = 2;
const int cstep = 6;
const int cstepx = cstep;
const int cstepy = cstep;
const int cstepz = cstep;

struct pointPsa{
  int   index;
  float pos[3];
  int   netChSeg;
  float Amp[NCHAN][NSIGS];
};

vector<pointPsa> *segPts[NSEGS];
vector<pointPsa> *segDev[NSEGS];
pointPsa averPt[NSEGS];
vector<int> cPtlist[NSEGS]; // coarse list
char hmask[NSEGS][NCHAN+1];
int nfinish = 0;

void FindPsa(int index, int &seg, int &ii){
  for(seg=0; seg<NSEGS; seg++){
    for(ii=0; ii<segPts[seg]->size(); ii++){
      if(segPts[seg]->at(ii).index==index) return;
      if(segPts[seg]->at(ii).index>index) break;
    }
  }

  seg = -1;
  ii = -1;
  return;
}

void MakeSegmentMap(int neighbours){
  cout<<"Make Segment Map...";
  memset(hmask, '0', sizeof(hmask));
  // '1' self, '2' neighbour, '9' CC
  for(int iseg=0; iseg<NSEGS; iseg++){
    int isec = iseg/NSLIC;
    int isli = iseg%NSLIC;

    for(int jseg=0; jseg<NSEGS; jseg++){
      int jsec = jseg/NSLIC;
      int jsli = jseg%NSLIC;

      int distV = abs(jsli - isli);
      int distH = abs(jsec - isec);
      distH = min(distH, abs(jsec - isec + NSECT));
      distH = min(distH, abs(jsec - isec - NSECT));
      int mdist = distV + distH;
      if(mdist<=neighbours && distH<neighbours && distV<neighbours){
        hmask[iseg][jseg] = (iseg==jseg) ? '1' : '2';
      }
    }

    hmask[iseg][NCHAN-1] = '9';
    hmask[iseg][NCHAN] = 0; // to close each line as a string
  }
  cout<<endl;
  return;
}


int Read(string basisfile){
  // read basis
  TFile *fbasis = new TFile(basisfile.c_str());
  if(!fbasis->IsOpen()){
    cerr<<"cannot find basisfile "<<basisfile<<endl;
    return 0;
  }
  cout<<"read pulse signal basis from "<<basisfile<<endl;

  Int_t seg;
  Double_t pos[3];
  Double_t core[121];
  Double_t spulse[4356];
  TTree *basistree = (TTree *)fbasis->Get("tree");
  basistree->SetBranchAddress("seg",&seg);
  basistree->SetBranchAddress("pos",pos);
  basistree->SetBranchAddress("core",core);
  basistree->SetBranchAddress("spulse",spulse);
  int npoint = basistree->GetEntriesFast();

  for(int iseg=0; iseg<NSEGS; iseg++)
    averPt[iseg].pos[0] = averPt[iseg].pos[1] = averPt[iseg].pos[2] = 0;

  for(int ipoint=0; ipoint<npoint; ipoint++){
    if(ipoint%1000==0) cout<<"\r ipoint = "<<ipoint<<flush;
    basistree->GetEntry(ipoint);

    seg = seg-1;

    pointPsa Pt;
    Pt.index = ipoint;
    Pt.netChSeg = seg;
    Pt.pos[0] = pos[0];  Pt.pos[1] = pos[1];  Pt.pos[2] = pos[2];

    averPt[seg].pos[0] += pos[0];
    averPt[seg].pos[1] += pos[1];
    averPt[seg].pos[2] += pos[2];

    for(int iseg=0; iseg<NSEGS; iseg++){
      for(int isig=0; isig<NSIGS; isig++)
        Pt.Amp[iseg][isig] = spulse[iseg*NSIGS+isig];
    }
    for(int isig=0; isig<NSIGS; isig++)
      Pt.Amp[NCHAN-1][isig] = core[isig];

    segPts[seg]->push_back(Pt);

    pointPsa PtDev;
    PtDev.index = ipoint;
    PtDev.netChSeg = seg;
    PtDev.pos[0] = pos[0];  PtDev.pos[1] = pos[1];  PtDev.pos[2] = pos[2];
    for(int iseg=0; iseg<NCHAN; iseg++){
      for(int isig=0; isig<NSIGS; isig++)
        PtDev.Amp[iseg][isig] = 0;
    }
    segDev[seg]->push_back(PtDev);
  }

  for(int iseg=0; iseg<NSEGS; iseg++){
    averPt[iseg].pos[0] = averPt[iseg].pos[0] / segPts[iseg]->size();
    averPt[iseg].pos[1] = averPt[iseg].pos[1] / segPts[iseg]->size();
    averPt[iseg].pos[2] = averPt[iseg].pos[2] / segPts[iseg]->size();
  }

  cout<<"\r load "<<npoint<<" points from basis"<<endl;
  fbasis->Close();
  return npoint;  
}


// input value of dist used to stop search when a point closer than this
inline int NearestPoint(float px, float py, float pz, int iseg, float& dist){
  vector<pointPsa> *ppts = segPts[iseg];
  int               npts = ppts->size();

  int cnear = -1;
  float cprec = dist*dist;

  float cdist = (float)1.e20;
  for(int jj = 0; jj < npts; jj++){
    float dx = ppts->at(jj).pos[0] - px;
    float dy = ppts->at(jj).pos[1] - py;
    float dz = ppts->at(jj).pos[2] - pz;
    float d2 = dx*dx + dy*dy + dz*dz;
    if(d2 < cdist){
      cdist = d2;
      cnear = jj;
      if(cdist < cprec)
        break;
    }
  }
  dist = (float)sqrt(cdist);
  return cnear;
}


// input value of dist used to stop search when a point closer than this
inline int NearestCoarse(float px, float py, float pz, int iseg, float& dist){
  vector<pointPsa> *ppts = segPts[iseg];
  int               numC = cPtlist[iseg].size();

  int cnear = -1;
  float cprec = dist*dist;

  float cdist = (float)1.e20;
  for(int kk = 0; kk < numC; kk++){
    int jj = cPtlist[iseg][kk];
    float dx = ppts->at(jj).pos[0] - px;
    float dy = ppts->at(jj).pos[1] - py;
    float dz = ppts->at(jj).pos[2] - pz;
    float d2 = dx*dx + dy*dy + dz*dz;
    if(d2 < cdist){
      cdist = d2;
      cnear = kk;
      if(cdist < cprec)
        break;
    }
  }
  dist = (float)sqrt(cdist);
  return cnear;
}


inline bool AssignIt(float px, float py, float pz, float cx, float cy, float cz){
  const float surface = 0.1f;  // to include the surface
  
  // all points inside a "cube" centered on this coarse grid point
  float dx = px - cx; if(dx < 0) dx = -dx;
  float dy = py - cy; if(dy < 0) dy = -dy;
  float dz = pz - cz; if(dz < 0) dz = -dz;
  bool inside = ( dx<(cstepx+surface) && dy<(cstepy+surface) && dz<(cstepz+surface));

  return inside;
}


inline bool IsCoarse(int jj, int seg){
  for(int cc=0; cc<cPtlist[seg].size(); cc++){
    if(cPtlist[seg][cc] == jj)
      return true;
  }
  return false;
}



void MakeCoarseFineList(int seg){
  vector<pointPsa> *ppts = segPts[seg];
  int               npts = ppts->size();

  // the (precalculated) center of the segment
  float xaver = averPt[seg].pos[0];
  float yaver = averPt[seg].pos[1];
  float zaver = averPt[seg].pos[2];

  // precision of the search for the nearest grid point
  float dprec = fstep*0.49f; // slightly less than half fine step

  float cdist = dprec;
  int   cnear = NearestPoint(xaver, yaver, zaver, seg, cdist);
  float ccx = ppts->at(cnear).pos[0];
  float ccy = ppts->at(cnear).pos[1];
  float ccz = ppts->at(cnear).pos[2];

  // find the bounding box of the segment
  float xmin = ccx; float xmax = ccx;
  float ymin = ccy; float ymax = ccy;
  float zmin = ccz; float zmax = ccz;
  for(int jj = 0; jj < npts; jj++){
    xmin = min(xmin, ppts->at(jj).pos[0]);  xmax = max(xmax, ppts->at(jj).pos[0]);
    ymin = min(ymin, ppts->at(jj).pos[1]);  ymax = max(ymax, ppts->at(jj).pos[1]);
    zmin = min(zmin, ppts->at(jj).pos[2]);  zmax = max(zmax, ppts->at(jj).pos[2]);
  }
  // expand the box to integer multiples of cstep with respect to the center of segment
  int   nn; float vn;
  const float prec = 0.1f;

  nn = (int)((ccx - xmin)/cstepx);
  vn = ccx-cstepx*nn;
  xmin = abs(vn-xmin)<prec ? vn : vn - cstepx;
  nn = (int)((xmax - ccx)/cstepx);
  vn = ccx+cstepx*nn;
  xmax = abs(vn-xmax)<prec ? vn : vn + cstepx;

  nn = (int)((ccy - ymin)/cstepy);
  vn = ccy-cstepy*nn;
  ymin = abs(vn-ymin)<prec ? vn : vn - cstepy;
  nn = (int)((ymax - ccy)/cstepy);
  vn = ccy+cstepy*nn;
  ymax = abs(vn-ymax)<prec ? vn : vn + cstepy;

  nn = (int)((ccz - zmin)/cstepz);
  vn = ccz-cstepz*nn;
  zmin = abs(vn-zmin)<prec ? vn : vn - cstepz;
  nn = (int)((zmax - ccz)/cstepz);
  vn = ccz+cstepz*nn;
  zmax = abs(vn-zmax)<prec ? vn : vn + cstepz;

  int npx = (int)((xmax-xmin + 2*prec)/cstepx)+1;
  int npy = (int)((ymax-ymin + 2*prec)/cstepy)+1;
  int npz = (int)((zmax-zmin + 2*prec)/cstepz)+1;
  int npc = npx*npy*npz;

  // Generate a temporary coarse grid with the same size as the fine-grid
  // This is much too big if cstep is larger than fstep and too small if cstep==fstep
  // In principle we should generate just npc elements but npts works always

  // to count assigment of points
  int *nAss = new int[npts]; memset(nAss, 0, sizeof(int  )*npts); // to define the coarse grid
  int *mAss = new int[npts]; memset(mAss, 0, sizeof(int  )*npts); // couts multiple assignments of a fgrid point

  int   cdx = 0, cdy = 0, cdz = 0;
  float ssx = 0, ssy = 0, ssz = 0;

  ssz = ccz;
  cdz = cstepz;
  while(true) { // z

    ssy = ccy;
    cdy = cstepy;
    while(true) { // y

      ssx = ccx;
      cdx = cstepx;
      while(true) { // x

	if(ssx >= xmin && ssx <= xmax && ssy >= ymin && ssy <= ymax && ssz >= zmin && ssz <= zmax) {
	  // Find the closest of the npts fine points, but only if the coarse point is inside the bounding box
	  float jdist = dprec; // stop looping when a point is closer than this
	  int jnear = NearestPoint(ssx, ssy, ssz, seg, jdist); // find point-id of (ssx,ssy,ssz)
	  if(jdist < dprec){ // 1 mm precision
	    if(!nAss[jnear]){ // already assigned
	      nAss[jnear]++;
              cPtlist[seg].push_back(jnear);
	    }
	  }
	}

	// update x and check boundaries
	float ssxn = ccx + cdx;
        if(cdx > 0){
          if(ssxn>xmax && ssx<xmin) break;
          cdx = -cdx;
        }else{
          if(ssxn<xmin && ssx>xmax) break;
          cdx = -cdx + cstepx;
        }
        ssx = ssxn;
	
      } // end of loop x

      // update y and check boundaries
      float ssyn = ccy + cdy;
      if(cdy > 0) {
        if(ssyn>ymax && ssy<ymin) break;
        cdy = -cdy;
      }
      else {
        if(ssyn<ymin && ssy>ymax) break;
        cdy = -cdy + cstepy;
      }
      ssy = ssyn;
      
    } // end of loop y

    // update z and check boundaries
    float sszn = ccz + cdz;
    if(cdz > 0) {
      if(sszn>zmax && ssz<zmin) break;
      cdz = -cdz;
    }
    else {
      if(sszn<zmin && ssz>zmax) break;
      cdz = -cdz + cstepz;
    }
    ssz = sszn;
    
  } // end of loop z


  // find the fine-grid points assigned to each point of the coarse grid
  // if needed, the coarse grid is extended to avoid unassigned points

  bool repeat = true;
  while(repeat) {

    // Count number of fine-grid points assigned to each of the coarse points
    // At the same time mark how many times the fine-grid points are assigned
    repeat = false;
    memset(mAss, 0, sizeof(int  )*npts);

    for(int kk = 0; kk<cPtlist[seg].size(); kk++){
      // count fine grid points, assigning to this coarse point a blob of surrounding points
      int jcc = cPtlist[seg][kk];
      float cpx = ppts->at(jcc).pos[0];
      float cpy = ppts->at(jcc).pos[1];
      float cpz = ppts->at(jcc).pos[2];
      for(int jj = 0; jj < npts; jj++){
	if( AssignIt(ppts->at(jj).pos[0], ppts->at(jj).pos[1], ppts->at(jj).pos[2], cpx, cpy, cpz) ){
	  // check that it is not a grid point
	  if( !IsCoarse(jj,seg) ){
            mAss[jj]++;
          }
	}
      }
      // add the coarse point itself
      mAss[jcc]++;
    }

    // Check if all segment points have been assigned.
    int notassigned = 0;
    float ddnew = 0;
    int   jjnew = -1;
    int   kknew = -1;
    for(int jj = 0; jj < npts; jj++){
      if(mAss[jj] == 0){
        // unassigned point
        float bx = ppts->at(jj).pos[0];
	float by = ppts->at(jj).pos[1];
	float bz = ppts->at(jj).pos[2];
	float kdist = dprec;
	int   knear = NearestCoarse(bx, by, bz, seg, kdist);
	if(jjnew == -1 || kdist > ddnew){
          jjnew = jj;
          kknew = knear;
          ddnew = kdist;
        }
        notassigned++;
      }
    }

    // If there are unassigned points, extend the coarse grid with the most detached of them.
    // Repeat the checks until all fine-grid points have been assigned
    if(notassigned){
      nAss[jjnew]++;
      cPtlist[seg].push_back(jjnew);
      repeat = true;
    }
    
  } // while(repeat)

  delete [] mAss;
  delete [] nAss;
  return;
}


void MakeCoarseDev(int seg){
  vector<pointPsa> *ppts = segDev[seg];
  int               npts = ppts->size();
  char            *lMask = hmask[seg];
  
  for(int jj=0; jj<npts; jj++){
    if( !IsCoarse(jj,seg) ) continue;

    for(int iseg=0; iseg<NCHAN; iseg++){

      if(lMask[iseg]=='0') continue;
      
      // FFT
      Double_t re_full[NSIGS];
      Double_t im_full[NSIGS];
      double k = -0.5;
      for(int isig=0; isig<NSIGS; isig++){
	if(isig==0 || isig>3){
	  re_full[isig] = im_full[isig] = 0;
	  continue;
	}
	re_full[isig] = gRandom->Uniform(-1,1)*exp(k*isig)/40;
	im_full[isig] = gRandom->Uniform(-1,1)*exp(k*isig)/40;

	if(isig==1){
	  re_full[isig] = re_full[isig]/2;
	  im_full[isig] = im_full[isig]/2;
	}else if(isig==2 || isig==3){
	  re_full[isig] = re_full[isig]*2;
	  im_full[isig] = im_full[isig]*2;
	}

	//re_full[0] += -2*re_full[isig];
	//im_full[0] += -2*im_full[isig];
      }

      // FFT backward transform
      TH1 *hb = 0;
      int nsig = NSIGS;
      TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &nsig, "C2R M K");
      fft_back->SetPointsComplex(re_full,im_full);
      fft_back->Transform();
      hb = TH1::TransformHisto(fft_back,hb,"Re");

      for(int isig=0; isig<NSIGS; isig++){
	double bincontent = hb->GetBinContent(isig);
	if(isig<25) bincontent = bincontent*exp(0.1*(isig-25.));
	if(isig>90) bincontent = bincontent*exp(0.1*(90.-isig));
        ppts->at(jj).Amp[iseg][isig] = bincontent;
      }
      delete hb;
    }

    nfinish++;
  }

  return;
}


void MakeCoarseDev2(int seg){
  vector<pointPsa> *ppts = segPts[seg];
  vector<pointPsa> *pdev = segDev[seg];
  int               npts = ppts->size();
  char            *lMask = hmask[seg];
  
  for(int jj=0; jj<npts; jj++){
    //if( !IsCoarse(jj,seg) ) continue;
    if(nfinish%100==0) cout<<"\r finish "<<nfinish<<" grid points..."<<flush;
    
    for(int iseg=0; iseg<NCHAN; iseg++){

      if(lMask[iseg]=='0') continue;

      // FFT
      int nsig = NSIGS2;
      Double_t amp[NSIGS2];
      for(int isig=0; isig<NSIGS2; isig++)
	amp[isig] = ppts->at(jj).Amp[iseg][NZERO+isig];
      
      TVirtualFFT *fft = TVirtualFFT::FFT(1, &nsig, "R2C M K");
      fft->SetPoints(amp);
      fft->Transform();
      Double_t re_full[NSIGS2];
      Double_t im_full[NSIGS2];
      fft->GetPointsComplex(re_full,im_full);

      for(int i=0; i<NSIGS2; i++){
	re_full[i] = re_full[i]/NSIGS2;
	im_full[i] = im_full[i]/NSIGS2;

	if(lMask[iseg]=='2'){ //neighbour
	  if(i==0){
	    re_full[i] = im_full[i] = 0;
	  }else if(i>0 && i<10){
	    re_full[i] = re_full[i] * gRandom->Gaus(1,0.5);
	    im_full[i] = im_full[i] * gRandom->Gaus(1,0.5);
	  }
	  re_full[0] += -2*re_full[i];
	  im_full[0] += -2*im_full[i];

	}else{ // self and core
	  if(i==0){
	    re_full[i] = im_full[i] = 0;
	  }else if(i>0 && i<4){
	    re_full[i] = re_full[i] * gRandom->Gaus(1,0.2);
	    im_full[i] = im_full[i] * gRandom->Gaus(1,0.2);
	  }
	  re_full[0] += -2*re_full[i];
	  im_full[0] += -2*im_full[i];
	}

      }

      // FFT backward transform
      TH1 *hb = 0;
      TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &nsig, "C2R M K");
      fft_back->SetPointsComplex(re_full,im_full);
      fft_back->Transform();
      hb = TH1::TransformHisto(fft_back,hb,"Re");

      for(int isig=0; isig<NSIGS; isig++){
	double bincontent = 0;
	if(isig>8 && isig<NSIGS2+1)
	  bincontent = hb->GetBinContent(isig-9) - ppts->at(jj).Amp[iseg][isig];

	if(isig>NSIGS2 && isig<NSIGS2+10){
	  bincontent =
	    (isig-NSIGS2)*ppts->at(jj).Amp[iseg][isig] +
	    (NSIGS2+10-isig)*hb->GetBinContent(isig-9);
	  bincontent = bincontent / 10.;
	  bincontent = bincontent - ppts->at(jj).Amp[iseg][isig];
	}
	
        pdev->at(jj).Amp[iseg][isig] = bincontent;
      }
      delete hb;
      delete fft;
      delete fft_back;
    }

    nfinish++;
  }

  return;
}


void MakeFineDev(int seg){
  vector<pointPsa> *ppts = segDev[seg];
  int               npts = ppts->size();
  char            *lMask = hmask[seg];
  
  int *nAss = new int[npts]; memset(nAss, 0, sizeof(int  )*npts); // to define the assigned grid
  for(int jj=0; jj<npts; jj++)
    if( IsCoarse(jj,seg) )
      nAss[jj] = 1;

  float maxcdist = cstep;
  bool repeat = true;
  int NAss = 0;
  while(repeat){
    
    for(int jj=0; jj<npts; jj++){
      if( nAss[jj]==1 ) continue;

      // if a gird still not assign Dev, find assigned grid in the same line
      vector <int>   clist;
      vector <float> cdist;
      bool kfound = false;
      int ix;
      for(ix=0; ix<3; ix++){
	clist.clear();
	cdist.clear();
	for(int ajj=0; ajj<npts; ajj++){
	  if( nAss[ajj]==0 ) continue;
	  int ix1 = (ix+1)%3;
	  int ix2 = (ix+2)%3;
	  if(fabs(ppts->at(jj).pos[ix]  - ppts->at(ajj).pos[ix]) > maxcdist) continue;
	  if(fabs(ppts->at(jj).pos[ix1] - ppts->at(ajj).pos[ix1]) > 0.1f) continue;
	  if(fabs(ppts->at(jj).pos[ix2] - ppts->at(ajj).pos[ix2]) > 0.1f) continue;
	  clist.push_back(ajj);
	  cdist.push_back(fabs(ppts->at(jj).pos[ix]  - ppts->at(ajj).pos[ix]));
	}
	if(clist.size()>1) kfound = true;
	if(kfound) break;
      }

      if(!kfound) continue; // cannot find assigned grid in the same line

      // sort according to dist
      for(int i=0; i<cdist.size(); i++){
	for(int ii=i+1; ii<cdist.size(); ii++){
	  if(cdist[i]>cdist[ii]){
	    int itmp = clist[i];
	    clist[i] = clist[ii];
	    clist[ii] = itmp;
	    float ftmp = cdist[i];
	    cdist[i] = cdist[ii];
	    cdist[ii] = ftmp;
	  }
	}
      }


      // interpolation or extrapolation
      float fpos, cpos[2];
      fpos = ppts->at(jj).pos[ix];
      cpos[0] = ppts->at(clist[0]).pos[ix];
      cpos[1] = ppts->at(clist[1]).pos[ix];
      float factor[2];
      factor[0] = (fpos - cpos[1]) / (cpos[0] - cpos[1]);
      factor[1] = (fpos - cpos[0]) / (cpos[1] - cpos[0]);
      for(int iseg=0; iseg<NCHAN; iseg++){
	if(lMask[iseg]=='0') continue;
	for(int isig=0; isig<NSIGS; isig++){
	  ppts->at(jj).Amp[iseg][isig] =
	    factor[0]*ppts->at(clist[0]).Amp[iseg][isig] +
	    factor[1]*ppts->at(clist[1]).Amp[iseg][isig];
	}
      }
      nAss[jj] = 1;
      
    }

    // check if all grid assigned
    repeat = false;
    int tmpNAss = 0;
    for(int jj=0; jj<npts; jj++){
      if(nAss[jj]==0) repeat = true;    
      else tmpNAss++;
    }

    // check if need extend maxcdist
    if(tmpNAss==NAss) maxcdist+=1;
    
    NAss = tmpNAss;
  } // while(repeat)
  
  return;
}


void WriteBasisFile(string basisfile, int npoint){
  TFile *fout = new TFile(basisfile.c_str(),"RECREATE");
  TTree *tree = new TTree("tree","tree with deviation");
  Int_t seg;
  Double_t pos[3];
  Double_t core[121];
  Double_t spulse[4356];
  tree->Branch("seg",&seg,"seg/I");
  tree->Branch("pos",pos,"pos[3]/D");
  tree->Branch("core",core,"core[121]/D");
  tree->Branch("spulse",spulse,"spulse[4356]/D");

  for(int i=0; i<npoint; i++){
    int ii;
    FindPsa(i,seg,ii);
    if(seg<0){ cerr<<"cannot find index "<<i<<endl; continue;}

    for(int ix=0; ix<3; ix++) pos[ix] = segPts[seg]->at(ii).pos[ix];

    for(int iseg=0; iseg<NSEGS; iseg++){
      for(int isig=0; isig<NSIGS; isig++)
	spulse[iseg*NSIGS+isig] =
	  segPts[seg]->at(ii).Amp[iseg][isig] +
	  segDev[seg]->at(ii).Amp[iseg][isig];
    }
    for(int isig=0; isig<NSIGS; isig++)
      core[isig] =
	segPts[seg]->at(ii).Amp[NCHAN-1][isig] +
	segDev[seg]->at(ii).Amp[NCHAN-1][isig];

    seg = seg+1;
    tree->Fill();
  }
  fout->cd();
  tree->Write();
  fout->Close();

  return;
}


void AddDev(string basisfile){
  MakeSegmentMap(2);
  for(int iseg=0; iseg<NSEGS; iseg++){
    segPts[iseg] = new vector<pointPsa>();
    segDev[iseg] = new vector<pointPsa>();
  }

  int npoint = Read(basisfile);

  cout<<"Generating coarse point list:";
  for(int iseg=0; iseg<NSEGS; iseg++)
    MakeCoarseFineList(iseg);
  cout<<endl;

  int ccount = 0;
  for(int iseg=0; iseg<NSEGS; iseg++)
    ccount += cPtlist[iseg].size();
  cout<<"Generate "<<ccount<<" coarse points in total"<<endl;

  time_t t;
  gRandom->SetSeed(time(&t));

  cout<<"Make Random Deviation for coarse grid:"<<endl;
  for(int iseg=0; iseg<NSEGS; iseg++){
    MakeCoarseDev2(iseg);  
  }
  cout<<"\r finish "<<nfinish<<" grid points..."<<endl;

  /*
  cout<<"Interpolation Deviation for fine grid:";
  for(int iseg=0; iseg<NSEGS; iseg++){
    MakeFineDev(iseg);  
  }
  cout<<endl;
  */

  cout<<"Write Basis file";
  string outputfile = "Dev/"+basisfile;
  WriteBasisFile(outputfile, npoint);
  cout<<endl;
  
  return;
}


#ifndef __CINT__
int main(int argc, char *argv[]){

  if(argc>1){
    AddDev( string(argv[1]) );
  }
  
  return 0;
}
#endif
