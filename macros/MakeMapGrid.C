// ./macros/MakeMapGrid dettype PSAfile mapname

#include <TRint.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TProfile.h>
#include "TSystem.h"

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

double maxdist = 1; // 1mm
double mindist = 0.1;
// chain
int MaxRun = -1;
float minNetCh = 900;

// segmentation
const Int_t NCryType = 3;
double griddist = 2; // 2mm grid
double range[NCryType][3][2]; // range of XYZ
const Int_t MaxSteps = 50;

double Map[NCryType][MaxSteps][MaxSteps][MaxSteps][3]; // position res map (sigma_phi,sigma_r,sigma_z)
double MapPos[NCryType][MaxSteps][MaxSteps][MaxSteps][3]; // position of map (x,y,z)
int    Seg[NCryType][MaxSteps][MaxSteps][MaxSteps]; // segment at grid point

// read map
void LoadMapGrid(string mapfilename){
  ifstream fin(mapfilename.c_str());
  if(!fin){
    gROOT->ProcessLine(Form(".!cp Map/MapPointsGrid.dat %s",mapfilename.c_str()));
    fin.open(mapfilename.c_str());
  }

  cout<<"\e[1;32m read position resolution map \e[0m"<<endl;

  // init
  for(int itype=0; itype<NCryType; itype++){
    for(int iaxis=0; iaxis<3; iaxis++){
      range[itype][iaxis][0] = 1;
      range[itype][iaxis][1] = 0;
    }
    for(int ix=0; ix<MaxSteps; ix++)
      for(int iy=0; iy<MaxSteps; iy++)
	for(int iz=0; iz<MaxSteps; iz++){
	  Map[itype][ix][iy][iz][0] = -1;
	  Map[itype][ix][iy][iz][1] = -1;
	  Map[itype][ix][iy][iz][2] = -1;
	  Seg[itype][ix][iy][iz] = -1;
	}
  }
  
  // read Map
  const int kMaxBufLen = 500;
  char *buffer = new char[kMaxBufLen];
  while(!fin.eof()){
    fin.getline(buffer,kMaxBufLen);

    if(strncmp(buffer,"#range",6)==0){
      cout<<"reading range for each crystal type"<<endl;
      fin.getline(buffer,kMaxBufLen);
      int itype;
      while(1){
        fin >> itype;
        if(itype==-1) break;
	int tmp;
        for(int iaxis=0; iaxis<3; iaxis++){
	  fin >> range[itype][iaxis][0] >> range[itype][iaxis][1] >> tmp;
	  if(tmp>MaxSteps) {cerr<<"change MaxSteps to >= "<<tmp<<endl; return;}
	}
      }
    }

    if(strncmp(buffer,"#Map",4)==0){
      cout<<"reading Map"<<endl;
      fin.getline(buffer,kMaxBufLen);
      int itype, iseg, ipos[3];
      double pos[3],res[3];
      while(1){
        fin >> itype >> iseg >> pos[0] >> pos[1] >> pos[2];
        if(itype==-1) break;

	for(int iaxis=0; iaxis<3; iaxis++){
	  if(pos[iaxis]-range[itype][iaxis][0]<-mindist ||
	     pos[iaxis]-range[itype][iaxis][1]>mindist){
	    cout<<Form("axis%d: %.3f outside range %.3f ~ %.3f",iaxis,pos[iaxis],range[itype][iaxis][0],range[itype][iaxis][1])<<endl;
	    ipos[iaxis] = -1;
	  }else{
	    ipos[iaxis] = (int)((pos[iaxis] - range[itype][iaxis][0]) / griddist + 0.5);
	    if(fabs(pos[iaxis]-(range[itype][iaxis][0]+griddist*ipos[iaxis]))>mindist){
	      cout<<Form("axis%d: %.3f not a grid point",iaxis,pos[iaxis])<<endl;
	      ipos[iaxis] = -1;
	    }
	  }
	}
	
        for(int i=0; i<3; i++) fin>>res[i];
	if(ipos[0]<0 || ipos[1]<0 || ipos[2]<0) continue;
        for(int i=0; i<3; i++){
	  Map[itype][ipos[0]][ipos[1]][ipos[2]][i] = res[i];
	  MapPos[itype][ipos[0]][ipos[1]][ipos[2]][i] = pos[i];
	}
	Seg[itype][ipos[0]][ipos[1]][ipos[2]] = iseg;
      }
    }

  }
  cout<<"finish load Map"<<endl;

  fin.close();
  return;
}

// write map
void WriteMapGrid(string mapfilename){

  cout<<"Write map..."<<endl;
  ofstream fout("Map/tmp.dat",ios::out);
  
  //output range
  fout<<"#range #####################"<<endl;
  fout<<"# itype  xmin  xmax  xsteps  ymin  ymax  ysteps  zmin  zmax  zsteps #####################"<<endl;
  for(int itype=0; itype<NCryType; itype++){
    fout<<Form("  %d",itype);
    for(int iaxis=0; iaxis<3; iaxis++){
      int steps = (int)((range[itype][iaxis][1]-range[itype][iaxis][0])/griddist+0.5)+1;
      fout<<Form("  %.3f  %.3f  %d",range[itype][iaxis][0],range[itype][iaxis][1],steps);
    }
    fout<<endl;
  }
  fout<<" -1  -1  -1  -1  -1  -1  -1  -1  -1  -1"<<endl;

  // output map
  fout<<"#Map #####################"<<endl;
  fout<<"# itype seg  x  y  z  sigma_phi sigma_r sigma_z #####################"<<endl;
  for(int itype=0; itype<NCryType; itype++)
    for(int iz=0; iz<MaxSteps; iz++)
      for(int iy=0; iy<MaxSteps; iy++)
	for(int ix=0; ix<MaxSteps; ix++){

            if(Seg[itype][ix][iy][iz]<0) continue;

	    fout<<Form("  %d  %d   %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
		       itype,
		       Seg[itype][ix][iy][iz],
		       MapPos[itype][ix][iy][iz][0],
		       MapPos[itype][ix][iy][iz][1],
		       MapPos[itype][ix][iy][iz][2],
		       Map[itype][ix][iy][iz][0],
		       Map[itype][ix][iy][iz][1],
		       Map[itype][ix][iy][iz][2])<<endl;
          }
  fout<<" -1  -1  -1  -1  -1  -1  -1  -1"<<endl;

  fout.close();  
  gROOT->ProcessLine(Form(".!mv -f Map/tmp.dat %s",mapfilename.c_str()));
  return;
}


// make Map
void MakeMapGrid(int itype, string PSAfile0, string mapfilename){
  if(itype<0 || itype>2) return;
  
  LoadMapGrid(mapfilename);

  // get resolution from PSA data
  TChain *tree = new TChain();
  for(int irun=0; ; irun++){
    if(MaxRun>0 && !(irun<MaxRun)) break;
    
    string PSAfile = (string)Form("%s%04d.root",PSAfile0.c_str(),irun);
    if(gSystem->AccessPathName(PSAfile.c_str())) break;
    
    tree->AddFile(PSAfile.c_str(),0,"tree");
  }

  
  Int_t numNetCharges;
  Int_t nhits;
  Int_t seg;
  Float_t NetCh;
  Float_t simpos[3];
  Float_t anapos[3];
  Float_t dist;
  tree->SetBranchAddress("numNetCharges",&numNetCharges);
  tree->SetBranchAddress("nhits",&nhits);
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("NetCh",&NetCh);
  tree->SetBranchAddress("simpos",simpos);
  tree->SetBranchAddress("anapos",anapos);
  tree->SetBranchAddress("dist",&dist);

  int nentries = tree->GetEntriesFast();
  cout<<"read "<<nentries<<" events from tree"<<endl;

  int nbins = MaxSteps*MaxSteps*MaxSteps;
  TProfile *pres[3];
  pres[0] = new TProfile("pres0","position resolution [phi]",nbins,0.5,nbins+0.5);
  pres[1] = new TProfile("pres1","position resolution [r]",nbins,0.5,nbins+0.5);
  pres[2] = new TProfile("pres2","position resolution [z]",nbins,0.5,nbins+0.5);
  TH1D *hres[3];
  hres[0] = new TH1D("hres0","position resolution [phi]",nbins,0.5,nbins+0.5);
  hres[1] = new TH1D("hres1","position resolution [r]",nbins,0.5,nbins+0.5);
  hres[2] = new TH1D("hres2","position resolution [z]",nbins,0.5,nbins+0.5);
  TH1D *htmp = new TH1D("htmp","ibin",nbins,0.5,nbins+0.5);
  
  int ievt;
  for(ievt=0; ievt<nentries; ievt++){
    if(ievt%1000==0)
      cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<flush;
    tree->GetEntry(ievt);
    
    if(NetCh<minNetCh) continue;
    
    TVector3 vecsim(simpos[0],simpos[1],0);  
    TVector3 vecana(anapos[0],anapos[1],0);  
    double simphi = vecsim.Phi();
    double simr = vecsim.Mag();
    double anaphi = vecana.Phi();
    double anar = vecana.Mag();
    double diffphi = anaphi-simphi;
    if(diffphi>TMath::Pi())  diffphi-=2*TMath::Pi();
    if(diffphi<-TMath::Pi()) diffphi+=2*TMath::Pi();

    for(int iz=0; iz<MaxSteps; iz++){
      int nextz = 0;
      for(int iy=0; iy<MaxSteps; iy++){
	if(nextz){ nextz=0; break;}
	for(int ix=0; ix<MaxSteps; ix++){

	  if(Map[itype][ix][iy][iz][0]<0) continue;

	  if(fabs(simpos[2]-MapPos[itype][ix][iy][iz][2])>maxdist){ nextz=1; break;} //z
	  if(fabs(simpos[1]-MapPos[itype][ix][iy][iz][1])>maxdist) break; //y
	  if(fabs(simpos[0]-MapPos[itype][ix][iy][iz][0])>maxdist) continue; //x

	  int ibin =
	    iz*MaxSteps*MaxSteps +
	    iy*MaxSteps +
	    ix;

	  htmp->Fill(ibin);
	  pres[0]->Fill(ibin,simr*diffphi);
	  pres[1]->Fill(ibin,anar-simr);
	  pres[2]->Fill(ibin,anapos[2]-simpos[2]);
	} // end of loop x
      } // end of loop y
    } // end of loop z
    
  }// end of loop evts

  // calc RMS
  cout<<"calculate resolution..."<<endl;
  for(int ibin=1; ibin<nbins+1; ibin++){
    for(int iaxis=0; iaxis<3; iaxis++){
      double binerr = pres[iaxis]->GetBinError(ibin)*sqrt(pres[iaxis]->GetBinEntries(ibin));
      if(binerr>0) hres[iaxis]->SetBinContent(ibin,binerr);
    }
  }

  // update map for itype
  cout<<"update map for type "<<itype<<" ..."<<endl;
  for(int iz=0; iz<MaxSteps; iz++)
    for(int iy=0; iy<MaxSteps; iy++)
      for(int ix=0; ix<MaxSteps; ix++){

	if(Seg[itype][ix][iy][iz]<0) continue;

	int ibin =
	  iz*MaxSteps*MaxSteps +
	  iy*MaxSteps +
	  ix;

	if(hres[0]->GetBinContent(ibin)>0)
	  for(int iaxis=0; iaxis<3; iaxis++)
	    Map[itype][ix][iy][iz][iaxis] = hres[iaxis]->GetBinContent(ibin);
	    
      }

  TFile *ftmp = new TFile("htmp.root","RECREATE");
  htmp->Write();
  for(int iaxis=0; iaxis<3; iaxis++) pres[iaxis]->Write();
  for(int iaxis=0; iaxis<3; iaxis++) hres[iaxis]->Write();
  ftmp->Close();

  WriteMapGrid(mapfilename);

  return;
}


#ifndef __CINT__
int main(int argc, char *argv[]){
  // clock
  time_t start, stop;
  time(&start);

  if(argc>3){
    MakeMapGrid( atoi(argv[1]), string(argv[2]), string(argv[3]));
  }else if(argc>2){
    MakeMapGrid( atoi(argv[1]), string(argv[2]), "Map/MapGrid.dat");
  }else{
    MakeMapGrid( 0, "rootfiles/typeA/AnaData_run", "Map/MapGrid.dat");
  }

  time(&stop);
  cout<<Form("============= Elapsed time: %.1f seconds =============",difftime(stop,start))<<endl;
  
  return 0;
}
#endif
