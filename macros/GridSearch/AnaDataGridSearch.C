//./macros/AnaData inputfile outputfile dbfile

#include "TRint.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TInterpreter.h"
#include "TRandom.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>

#include "Global.hh"
#include "SignalBasis.hh"
#include "PSAFilterGridSearch.hh"

using namespace std;

void AnaData(string inputfile, string outputfile, string basisfile){

  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");
  gInterpreter->GenerateDictionary("vector<vector<float>>","vector");

  // PSA
  PSAFilterGridSearch *fpsa = new PSAFilterGridSearch(basisfile);

  // read data
  TFile *fin = new TFile(inputfile.c_str());
  if(!fin->IsOpen()){
    cerr<<"cannot find inputfile "<<inputfile<<endl;
    return;
  }
  TTree *intree = (TTree *)fin->Get("tree");
  const int nsig = 121; // npoints in 5ns step
  // Declaration of leaf types                                                                     
  int                      ievent;     // event id
  vector<int>             *ndet = 0;   // detector id
  vector<int>             *g4seg = 0;  // segment id
  vector<float>           *energy = 0; // deposit energy
  vector<vector<float>>   *posa = 0;   // absolute/global interaction position vector<float(3)>
  vector<vector<float>>   *posr = 0;   // relative/local interaction position vector<float(3)>

  // pulse shape vector<>
  vector<int>             *pdet = 0;   // detector id for pulse shape
  vector<float>           *ecore = 0;  // core energy
  vector<vector<int>>     *inter = 0;  // interaction id in G4 vector<>
    
  vector<vector<int>>     *pseg = 0;   // segment id in pulsedb
  vector<vector<int>>     *vngrid = 0;  // number of grid found around ppos
  vector<vector<int>>     *vextrpl = 0;  // number of grid for extrapolation
  vector<vector<float>>   *core = 0;   // core pulse shape vector<float(121)>
  vector<vector<float>>   *spulse = 0; // segment pulse shape vector<float(4356)>
    
  int                      category; // 1: max 1 seg fired in a det, >1 det fired; 2: max >1 seg fired in a det

  // set branch addresses and branch pointers
  intree->SetBranchAddress("ievent",&ievent);
  intree->SetBranchAddress("ndet",&ndet);
  intree->SetBranchAddress("g4seg",&g4seg);
  intree->SetBranchAddress("energy",&energy);
  intree->SetBranchAddress("posa",&posa);
  intree->SetBranchAddress("posr",&posr);

  intree->SetBranchAddress("pdet",&pdet);
  intree->SetBranchAddress("ecore",&ecore);
  intree->SetBranchAddress("inter",&inter);

  intree->SetBranchAddress("pseg",&pseg);
  intree->SetBranchAddress("ngrid",&vngrid);
  intree->SetBranchAddress("extrpl",&vextrpl);
  intree->SetBranchAddress("core",&core);
  intree->SetBranchAddress("spulse",&spulse);

  //intree->SetBranchAddress("category",&category);


  int nentries = intree->GetEntriesFast();
  cout<<"read "<<nentries<<" events from "<<inputfile<<endl;


#ifdef NOISE
  cout<<"AddNoise: create noise base size "<<NOISE<<endl;
  // get random noise
  float noise[NOISE];
  float noise0[NOISE];
  float tmpnoise;
  for(int i=0; i<NOISE; ){
    tmpnoise = gRandom->Uniform(-1,1);
    for(int isig=0; isig<nsig; isig++){
      if(i>=NOISE) break;
      tmpnoise += -0.2*(tmpnoise+gRandom->Uniform(-10,10));
      noise[i] = noise0[i] = tmpnoise;
      i++;
    }
  }

  // smooth noise
  int nsmooth = 3;
  for(int i=nsmooth; i<NOISE-nsmooth; i++){
    for(int j=1; j<=nsmooth; j++)
      noise[i] += noise0[i-j] + noise0[i+j];

    noise[i] = noise[i] / (2*nsmooth+1) * 2; // scale noise
  }
#endif

  
  // output
  TFile *fout = new TFile(outputfile.c_str(),"RECREATE");
  TTree *anatree = new TTree("tree","analyzed tree");
  Int_t   numNetCharges;
  Int_t   nhits;
  Int_t   seg;
  Float_t NetCh;
  Int_t   ngrid;
  Int_t   extrpl;
  Float_t simpos[3];
  Float_t anapos[3];
  Float_t dist;
  anatree->Branch("numNetCharges", &numNetCharges, "numNetCharges/I");
  anatree->Branch("nhits", &nhits, "nhits/I");
  anatree->Branch("seg", &seg, "seg/I");
  anatree->Branch("NetCh", &NetCh, "NetCh/F");
  //anatree->Branch("ngrid", &ngrid, "ngrid/I");
  //anatree->Branch("extrpl", &extrpl, "extrpl/I");
  anatree->Branch("simpos", simpos, "simpos[3]/F");
  anatree->Branch("anapos", anapos, "anapos[3]/F");
  anatree->Branch("dist", &dist, "dist/F");


  cout<<"gUseAdaptive = "<<gUseAdaptive<<"; "
      <<"gCoarseOnly  = "<<gCoarseOnly<<"; "
      <<endl;
  
  // clock
  time_t start, stop;
  time(&start);

  int ievt;
  for(ievt=0; ievt<nentries; ievt++){
    if(ievt%1000==0)
      cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<flush;

    intree->GetEntry(ievt);
    for(int idet=0; idet<pdet->size(); idet++){
      if(pdet->at(idet)!=0) continue; // only work on det 0

      //*************************************
      // calc average position
      //*************************************
      vector<int>           simseg;
      vector<int>           simnhits;
      vector<float>         simeng;
      vector<vector<float>> simdetpos;

      // get individual hits in one det
      for(int i=0; i<inter->at(idet).size(); i++){
	int interid = inter->at(idet)[i];

	simseg.push_back(pseg->at(idet)[i]-1); // pseg start from 1
	simnhits.push_back(1);
	simeng.push_back(energy->at(interid));
	vector<float> tmpdetpos;
	for(int ix=0; ix<3; ix++) tmpdetpos.push_back(posr->at(interid)[ix]); // pos in det frame
	simdetpos.push_back(tmpdetpos);
      }

      // sum up in a segment
      for(int i=0; i<simeng.size(); i++){
	for(int j=i+1; j<simeng.size(); j++){
	  if(simseg[i]==simseg[j]){
	    simnhits[i] = simnhits[i] + simnhits[j];
	    float tmpe = simeng[i] + simeng[j];
	    for(int ix=0; ix<3; ix++)
	      simdetpos[i][ix] = simeng[i]/tmpe*simdetpos[i][ix] + simeng[j]/tmpe*simdetpos[j][ix];
	    simeng[i] = tmpe;

	    simseg.erase(simseg.begin()+j);
	    simnhits.erase(simnhits.begin()+j);
	    simeng.erase(simeng.begin()+j);
	    simdetpos.erase(simdetpos.begin()+j);
	    j--;
	  }
	}
      }

      //*************************************
      // look at pulse shape
      //*************************************
      pointExp pE;
      for(int iseg=0; iseg<NSEGS; iseg++){
	for(int isig=0; isig<BSIZE; isig++){
	  pE.tAmp[iseg][isig] = (spulse->at(idet)[iseg*nsig+basis_tzero+isig*2] + spulse->at(idet)[iseg*nsig+basis_tzero+isig*2+1])/2;
	}
      }
      for(int isig=0; isig<BSIZE; isig++){
	pE.tAmp[NCHAN-1][isig] = (core->at(idet)[basis_tzero+isig*2] + core->at(idet)[basis_tzero+isig*2+1])/2;
      }

#ifdef NOISE
      for(int iseg=0; iseg<NCHAN; iseg++){
	int nidx = (int)gRandom->Uniform(0,NOISE);
	for(int isig=0; isig<BSIZE; isig++){
	  nidx = nidx%NOISE;
	  pE.tAmp[iseg][isig] += noise[nidx];
	  nidx+=2;
	}
      }
#endif


      // find hit segments
      int numsegs = 0;
      for(int iseg=0; iseg<NSEGS; iseg++){
	float segE = 0;
	for(int isig=BSIZE-2; isig>BSIZE-7; isig--) segE += pE.tAmp[iseg][isig];
	segE = segE/5.;
	if(segE > minSegEnergy){
	  pE.netChargeSegnum[numsegs] = iseg;
	  pE.netChargeEnergy[numsegs] = segE;
	  pE.nhits[numsegs] = -1;
	  for(int ii=0; ii<simseg.size(); ii++){
	    if(simseg[ii]==iseg){
	      pE.nhits[numsegs] = simnhits[ii];
	      pE.x[numsegs] = simdetpos[ii][0];
	      pE.y[numsegs] = simdetpos[ii][1];
	      pE.z[numsegs] = simdetpos[ii][2];
	      if(fabs(segE-simeng[ii])>2*minSegEnergy){
		cout<<"seg "<<iseg<<" nhits "<<simnhits[ii]
		    <<" : segE = "<<segE<<" ; simeng = "<<simeng[ii]<<endl;
		cout<<"tAmp: ";
		for(int isig=BSIZE-2; isig>BSIZE-7; isig--)
		  cout<<Form("%.0f ",pE.tAmp[iseg][isig]);
		cout<<endl;

	      }
	    }
	  }
	  numsegs++;
	}
      }

      pE.numNetCharges = numsegs;
      pE.netChSeg = -1;
      pE.bestPt = -1;
      pE.chi2min = float(1.e30);
      pE.isValid = false;
      pE.isInitialized = false;
      

      //*************************************
      // Coarse-Fine GridSearch    
      //*************************************
      fpsa->ProcessOneEvent(pE);


      //*************************************
      // output results
      //*************************************
      numNetCharges = pE.numNetCharges;
      for(int snum=0; snum<numsegs; snum++){

	nhits = pE.nhits[snum];
	seg   = pE.netChargeSegnum[snum];
	NetCh = pE.netChargeEnergy[snum];
	simpos[0] = pE.x[snum]; simpos[1] = pE.y[snum]; simpos[2] = pE.z[snum];
	
	int bestPt = pE.resPt[snum];
	fpsa->GetPtPos(seg, bestPt, anapos);

	dist = 0;
	for(int ix=0; ix<3; ix++){
	  dist += pow(anapos[ix]-simpos[ix] , 2);
	}
	dist = sqrt(dist);

	anatree->Fill();
      }
      
    } // end of loop idet

  } // end of loop ievt
  cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<endl;

  time(&stop);
  printf("============ Elapsed time: %.1f seconds =============\n",difftime(stop,start));

  cout<<"write to "<<outputfile<<" ..."<<endl;
  fout->cd();
  anatree->Write();
  fout->Close();
  fin->Close();

  return;
}

#ifndef __CINT__
int main(int argc, char *argv[]){

  if(argc>3){
    AnaData( string(argv[1]), string(argv[2]), string(argv[3]) );
  }else if(argc>2){
    AnaData( string(argv[1]), string(argv[2]), "pulsedb/pulseA.root" );
  }else if(argc>1){
    AnaData( string(argv[1]), "rootfiles/AnaData.root", "pulsedb/pulseA.root" );
  }else{
    AnaData( "rootfiles/G4SimData.root", "rootfiles/G4AnaData.root", "pulsedb/pulseA.root" );
  }

  return 0;
}
#endif
