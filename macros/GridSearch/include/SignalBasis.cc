#include "TFile.h"
#include "TTree.h"
#include "SignalBasis.hh"

SignalBasis::SignalBasis(){
  MakeSegmentMap(2);
  for(int iseg=0; iseg<NSEGS; iseg++){
    segPts[iseg] = new vector<pointPsa>();
    Ptlist[iseg] = new vector<vector<int>>();
  }
}

void SignalBasis::Read(string basisfile){
  // read basis
  TFile *fbasis = new TFile(basisfile.c_str());
  if(!fbasis->IsOpen()){
    cerr<<"cannot find basisfile "<<basisfile<<endl;
    return;
  }
  cout<<"read pulse signal basis from "<<basisfile<<endl;

  const int nsig = 121; // npoints in 5ns step
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
    averPt[iseg].x = averPt[iseg].y = averPt[iseg].z = 0;
  
  for(int ipoint=0; ipoint<npoint; ipoint++){
    if(ipoint%1000==0) cout<<"\r ipoint = "<<ipoint<<flush;
    basistree->GetEntry(ipoint);

    seg = seg-1;

    pointPsa Pt;
    Pt.netChSeg = seg;
    Pt.x = pos[0];  Pt.y = pos[1];  Pt.z = pos[2];

    averPt[seg].x += pos[0];
    averPt[seg].y += pos[1];
    averPt[seg].z += pos[2];
    
    for(int iseg=0; iseg<NSEGS; iseg++){
      for(int isig=0; isig<BSIZE; isig++)
	Pt.Amp[iseg][isig] = (spulse[iseg*nsig+basis_tzero+isig*2] + spulse[iseg*nsig+basis_tzero+isig*2+1])/2;
    }
    for(int isig=0; isig<BSIZE; isig++)
      Pt.Amp[NCHAN-1][isig] = (core[basis_tzero+isig*2] + core[basis_tzero+isig*2+1])/2;
    
    segPts[seg]->push_back(Pt);
  }

  for(int iseg=0; iseg<NSEGS; iseg++){
    averPt[iseg].x = averPt[iseg].x / segPts[iseg]->size();
    averPt[iseg].y = averPt[iseg].y / segPts[iseg]->size();
    averPt[iseg].z = averPt[iseg].z / segPts[iseg]->size();
  }

  cout<<"\r load "<<npoint<<" points from basis"<<endl;
  fbasis->Close();
  return;
}


void SignalBasis::MakeCoarseFineList(){
  cout<<"Generating coarse-fine list:";
  for(int iseg=0; iseg<NSEGS; iseg++)
    MakeCoarseFineList(iseg);
  cout<<endl;

  int ccount = 0;
  for(int iseg=0; iseg<NSEGS; iseg++)
    ccount += Ptlist[iseg]->size();
  cout<<"Generate "<<ccount<<" coarse points in total"<<endl;
    
  return;
}


void SignalBasis::MakeCoarseFineList(int seg){
  vector<pointPsa> *ppts = segPts[seg];
  int               npts = ppts->size();

  // the (precalculated) center of the segment
  float xaver = averPt[seg].x;
  float yaver = averPt[seg].y;
  float zaver = averPt[seg].z;

  // precision of the search for the nearest grid point
  float dprec = fstep*0.49f; // slightly less than half fine step

  float cdist = dprec;
  int   cnear = NearestPoint(xaver, yaver, zaver, seg, cdist);
  float ccx = ppts->at(cnear).x;
  float ccy = ppts->at(cnear).y;
  float ccz = ppts->at(cnear).z;
  
  // find the bounding box of the segment
  float xmin = ccx; float xmax = ccx;
  float ymin = ccy; float ymax = ccy;
  float zmin = ccz; float zmax = ccz;
  for(int jj = 0; jj < npts; jj++){
    xmin = min(xmin, ppts->at(jj).x);  xmax = max(xmax, ppts->at(jj).x);
    ymin = min(ymin, ppts->at(jj).y);  ymax = max(ymax, ppts->at(jj).y);
    zmin = min(zmin, ppts->at(jj).z);  zmax = max(zmax, ppts->at(jj).z);
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
	      vector<int> clist;
	      clist.push_back(jnear);
	      Ptlist[seg]->push_back(clist);	      
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

    for(int kk = 0; kk<Ptlist[seg]->size(); kk++){
      // count fine grid points, assigning to this coarse point a blob of surrounding points
      int jcc = Ptlist[seg]->at(kk)[0];
      float cpx = ppts->at(jcc).x;
      float cpy = ppts->at(jcc).y;
      float cpz = ppts->at(jcc).z;
      for(int jj = 0; jj < npts; jj++){
	if( AssignIt(ppts->at(jj).x, ppts->at(jj).y, ppts->at(jj).z, cpx, cpy, cpz) ){
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
	float bx = ppts->at(jj).x;
	float by = ppts->at(jj).y;
	float bz = ppts->at(jj).z;
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
      vector<int> clist;
      clist.push_back(jjnew);
      Ptlist[seg]->push_back(clist);	      
      repeat = true;
    }
    
  } // while(repeat)

  // prepare the Coarse-Fine structure for this segment


  // now fill it

  // multiple assignments as calculated above
  for(int indC = 0; indC < Ptlist[seg]->size(); indC++){
    int jcc = Ptlist[seg]->at(indC)[0];
    float cpx = ppts->at(jcc).x;
    float cpy = ppts->at(jcc).y;
    float cpz = ppts->at(jcc).z;
    for(int jj = 0; jj < npts; jj++){
      if(jj == jcc){
	continue; // pivot already inserted
      }
      if( AssignIt(ppts->at(jj).x, ppts->at(jj).y, ppts->at(jj).z, cpx, cpy, cpz) ){
	if( !IsCoarse(jj,seg) ){
	  Ptlist[seg]->at(indC).push_back(jj);
	}
      }
    }
  }

  //cout<<"seg "<<seg<<": "<<Ptlist[seg]->size()<<" coarse points"<<endl;
  
  delete [] mAss;
  delete [] nAss;
  return;
}


void SignalBasis::MakeSegmentMap(int neighbours){
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

// input value of dist used to stop search when a point closer than this
inline int SignalBasis::NearestPoint(float px, float py, float pz, int iseg, float& dist){
  vector<pointPsa> *ppts = segPts[iseg];
  int               npts = ppts->size();
  
  int cnear = -1;
  float cprec = dist*dist;

  float cdist = (float)1.e20;
  for(int jj = 0; jj < npts; jj++){
    float dx = ppts->at(jj).x - px;
    float dy = ppts->at(jj).y - py;
    float dz = ppts->at(jj).z - pz;
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
inline int SignalBasis::NearestCoarse(float px, float py, float pz, int iseg, float& dist){
  vector<pointPsa> *ppts = segPts[iseg];
  int               numC = Ptlist[iseg]->size();
  
  int cnear = -1;
  float cprec = dist*dist;

  float cdist = (float)1.e20;
  for(int kk = 0; kk < numC; kk++){
    int jj = Ptlist[iseg]->at(kk)[0];
    float dx = ppts->at(jj).x - px;
    float dy = ppts->at(jj).y - py;
    float dz = ppts->at(jj).z - pz;
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


inline bool SignalBasis::AssignIt(float px, float py, float pz, float cx, float cy, float cz){
  const float surface = 0.1f;  // to include the surface

  // all points inside a "cube" centered on this coarse grid point
  float dx = px - cx; if(dx < 0) dx = -dx;
  float dy = py - cy; if(dy < 0) dy = -dy;
  float dz = pz - cz; if(dz < 0) dz = -dz;
  bool inside = ( dx<(cstepx+surface) && dy<(cstepy+surface) && dz<(cstepz+surface));
  
  return inside;
}


inline bool SignalBasis::IsCoarse(int jj, int seg){
  for(int cc=0; cc<Ptlist[seg]->size(); cc++){
    if(Ptlist[seg]->at(cc)[0] == jj)
      return true;
  }
  return false;
}
