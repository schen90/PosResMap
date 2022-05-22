#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TVector3.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

const int nbinsxy = 100;
const double xymin = -50;
const double xymax = 50;
const int nbinsz = 100;
const double zmin = -5;
const double zmax = 95;

const int nseg = 36;
const int nsli = 6;
const int nsec = 6;

void MakePlotChain(){
  TChain *tree = new TChain();
  for(int irun=1; irun<=16; irun++)
    tree->AddFile(Form("rootfiles/AnaData_db2.0/AnaData_run%04d.root",irun),0,"tree");

  Int_t seg;
  Int_t ngrid;
  Int_t extrpl;
  Double_t simpos[3];
  Double_t anapos[3];
  Double_t dist;
  Double_t minsum2;
  tree->SetBranchAddress("simseg",&seg);
  tree->SetBranchAddress("ngrid",&ngrid);
  tree->SetBranchAddress("extrpl",&extrpl);
  tree->SetBranchAddress("simpos",simpos);
  tree->SetBranchAddress("anapos",anapos);
  tree->SetBranchAddress("dist",&dist);
  tree->SetBranchAddress("minsum2",&minsum2);

  int nentries = tree->GetEntriesFast();
  cout<<"read "<<nentries<<" events from tree"<<endl;

  // define hist
  TFile *fout = new TFile("draw/plot/plotchain.root","RECREATE");
  TH3F *hxyz = new TH3F("hxyz","simulated xyz position",nbinsxy/2,xymin,xymax,nbinsxy/2,xymin,xymax,nbinsz/2,zmin,zmax);
  TH2F *hxy = new TH2F("hxy","x vs. y",nbinsxy/2,xymin,xymax,nbinsxy/2,xymin,xymax);
  TH2F *hxz = new TH2F("hxz","x vs. z",nbinsxy/2,xymin,xymax,nbinsz/2,zmin,zmax);
  TH2F *hyz = new TH2F("hyz","y vs. z",nbinsxy/2,xymin,xymax,nbinsz/2,zmin,zmax);

  TProfile3D *pxyz0[3];
  pxyz0[0] = new TProfile3D("pxyz0_0","position resolution [rphi]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
  pxyz0[1] = new TProfile3D("pxyz1_0","position resolution [r]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
  pxyz0[2] = new TProfile3D("pxyz2_0","position resolution [z]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);

  TProfile2D *pxy0[3];
  pxy0[0] = new TProfile2D("pxy0_0","position resolution [rphi]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
  pxy0[1] = new TProfile2D("pxy0_1","position resolution [r]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
  pxy0[2] = new TProfile2D("pxy0_2","position resolution [z]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
  TProfile2D *pxz0[3];
  pxz0[0] = new TProfile2D("pxz0_0","position resolution [rphi]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
  pxz0[1] = new TProfile2D("pxz0_1","position resolution [r]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
  pxz0[2] = new TProfile2D("pxz0_2","position resolution [z]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);

  
  TH3D *pxyz[3];
  pxyz[0] = new TH3D("pxyz0","position resolution [rphi]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
  pxyz[1] = new TH3D("pxyz1","position resolution [r]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
  pxyz[2] = new TH3D("pxyz2","position resolution [z]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);

  TH2D *pxy[3];
  pxy[0] = new TH2D("pxy0","position resolution [rphi]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
  pxy[1] = new TH2D("pxy1","position resolution [r]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
  pxy[2] = new TH2D("pxy2","position resolution [z]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
  TH2D *pxz[3];
  pxz[0] = new TH2D("pxz0","position resolution [rphi]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
  pxz[1] = new TH2D("pxz1","position resolution [r]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
  pxz[2] = new TH2D("pxz2","position resolution [z]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);

  TProfile3D *psli0[nsli][3];
  TH3D *psli[nsli][3];

  TProfile2D *pxysli0[nsli][3];
  TH2D *pxysli[nsli][3];
  for(int i=0; i<nsli; i++){
    psli0[i][0] = new TProfile3D(Form("psli%d_0_0",i),"position resolution [rphi]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
    psli0[i][1] = new TProfile3D(Form("psli%d_1_0",i),"position resolution [r]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
    psli0[i][2] = new TProfile3D(Form("psli%d_2_0",i),"position resolution [z]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);

    psli[i][0] = new TH3D(Form("psli%d_0",i),"position resolution [rphi]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
    psli[i][1] = new TH3D(Form("psli%d_1",i),"position resolution [r]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
    psli[i][2] = new TH3D(Form("psli%d_2",i),"position resolution [z]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);

    pxysli0[i][0] = new TProfile2D(Form("pxysli%d_0_0",i),"position resolution [rphi]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
    pxysli0[i][1] = new TProfile2D(Form("pxysli%d_1_0",i),"position resolution [r]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
    pxysli0[i][2] = new TProfile2D(Form("pxysli%d_2_0",i),"position resolution [z]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);

    pxysli[i][0] = new TH2D(Form("pxysli%d_0",i),"position resolution [rphi]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
    pxysli[i][1] = new TH2D(Form("pxysli%d_1",i),"position resolution [r]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
    pxysli[i][2] = new TH2D(Form("pxysli%d_2",i),"position resolution [z]",nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);
  }
  
  TProfile2D *pxzsec0[nsec][3];
  TH2D *pxzsec[nsec][3];
  for(int i=0; i<nsec; i++){
    pxzsec0[i][0] = new TProfile2D(Form("pxzsec%d_0_0",i),"position resolution [rphi]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
    pxzsec0[i][1] = new TProfile2D(Form("pxzsec%d_1_0",i),"position resolution [r]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
    pxzsec0[i][2] = new TProfile2D(Form("pxzsec%d_2_0",i),"position resolution [z]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);

    pxzsec[i][0] = new TH2D(Form("pxzsec%d_0",i),"position resolution [rphi]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
    pxzsec[i][1] = new TH2D(Form("pxzsec%d_1",i),"position resolution [r]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
    pxzsec[i][2] = new TH2D(Form("pxzsec%d_2",i),"position resolution [z]",nbinsxy,xymin,xymax,nbinsz,zmin,zmax);
  }
  
  int ievt;
  for(ievt=0; ievt<nentries; ievt++){
    if(ievt%1000==0)
      cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<flush;
    tree->GetEntry(ievt);
    seg = seg-1;
    int isli = seg%6;
    int isec = (int)(seg/6);

    TVector3 vecsim(simpos[0],simpos[1],0);
    TVector3 vecana(anapos[0],anapos[1],0);
    double simphi = vecsim.Phi();
    double simr = vecsim.Mag();
    double anaphi = vecana.Phi();
    double anar = vecana.Mag();
    double diffphi = anaphi-simphi;
    if(diffphi>TMath::Pi()) diffphi-= 2*TMath::Pi();
    if(diffphi<-TMath::Pi()) diffphi+= 2*TMath::Pi();
    
    hxyz->Fill(simpos[0],simpos[1],simpos[2]);
    hxy->Fill(simpos[0],simpos[1]);
    hxz->Fill(simpos[0],simpos[2]);
    hyz->Fill(simpos[1],simpos[2]);

    pxyz0[0]->Fill(simpos[0],simpos[1],simpos[2],simr*diffphi);
    pxyz0[1]->Fill(simpos[0],simpos[1],simpos[2],anar-simr);
    pxyz0[2]->Fill(simpos[0],simpos[1],simpos[2],anapos[2]-simpos[2]);

    pxy0[0]->Fill(simpos[0],simpos[1],simr*diffphi);
    pxy0[1]->Fill(simpos[0],simpos[1],anar-simr);
    pxy0[2]->Fill(simpos[0],simpos[1],anapos[2]-simpos[2]);

    pxz0[0]->Fill(simpos[0],simpos[2],simr*diffphi);
    pxz0[1]->Fill(simpos[0],simpos[2],anar-simr);
    pxz0[2]->Fill(simpos[0],simpos[2],anapos[2]-simpos[2]);

    // slices
    psli0[isli][0]->Fill(simpos[0],simpos[1],simpos[2],simr*diffphi);
    psli0[isli][1]->Fill(simpos[0],simpos[1],simpos[2],anar-simr);
    psli0[isli][2]->Fill(simpos[0],simpos[1],simpos[2],anapos[2]-simpos[2]);

    pxysli0[isli][0]->Fill(simpos[0],simpos[1],simr*diffphi);
    pxysli0[isli][1]->Fill(simpos[0],simpos[1],anar-simr);
    pxysli0[isli][2]->Fill(simpos[0],simpos[1],anapos[2]-simpos[2]);

    // sectors
    pxzsec0[isec][0]->Fill(simpos[0],simpos[2],simr*diffphi);
    pxzsec0[isec][1]->Fill(simpos[0],simpos[2],anar-simr);
    pxzsec0[isec][2]->Fill(simpos[0],simpos[2],anapos[2]-simpos[2]);
  }
  cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<endl;

  for(int ix=1; ix<nbinsxy+1; ix++)
    for(int iy=1; iy<nbinsxy+1; iy++)
      for(int iz=1; iz<nbinsz+1; iz++){
	for(int iaxis=0; iaxis<3; iaxis++){
	  int nbin = pxyz0[iaxis]->GetBin(ix,iy,iz);
	  double binerr = pxyz0[iaxis]->GetBinError(ix,iy,iz)*sqrt(pxyz0[iaxis]->GetBinEntries(nbin));
	  if(binerr>0) pxyz[iaxis]->SetBinContent(ix,iy,iz,binerr);

	  for(int i=0; i<nsli; i++){
	    nbin = psli0[i][iaxis]->GetBin(ix,iy,iz);
	    binerr = psli0[i][iaxis]->GetBinError(ix,iy,iz)*sqrt(psli0[i][iaxis]->GetBinEntries(nbin));
	    if(binerr>0) psli[i][iaxis]->SetBinContent(ix,iy,iz,binerr);
	  }
	}
      }

  for(int ix=1; ix<nbinsxy+1; ix++)
    for(int iy=1; iy<nbinsxy+1; iy++)
      for(int iaxis=0; iaxis<3; iaxis++){
	int nbin = pxy0[iaxis]->GetBin(ix,iy);
	double binerr = pxy0[iaxis]->GetBinError(ix,iy)*sqrt(pxy0[iaxis]->GetBinEntries(nbin));
	if(binerr>0) pxy[iaxis]->SetBinContent(ix,iy,binerr);

	for(int i=0; i<nsli; i++){
	  nbin = pxysli0[i][iaxis]->GetBin(ix,iy);
	  binerr = pxysli0[i][iaxis]->GetBinError(ix,iy)*sqrt(pxysli0[i][iaxis]->GetBinEntries(nbin));
	  if(binerr>0) pxysli[i][iaxis]->SetBinContent(ix,iy,binerr);
	}
      }

  for(int ixy=1; ixy<nbinsxy+1; ixy++)
    for(int iz=1; iz<nbinsz+1; iz++)
      for(int iaxis=0; iaxis<3; iaxis++){
	int nbin = pxz0[iaxis]->GetBin(ixy,iz);
	double binerr = pxz0[iaxis]->GetBinError(ixy,iz)*sqrt(pxz0[iaxis]->GetBinEntries(nbin));
	if(binerr>0) pxz[iaxis]->SetBinContent(ixy,iz,binerr);

	for(int i=0; i<nsec; i++){
	  nbin = pxzsec0[i][iaxis]->GetBin(ixy,iz);
	  binerr = pxzsec0[i][iaxis]->GetBinError(ixy,iz)*sqrt(pxzsec0[i][iaxis]->GetBinEntries(nbin));
	  if(binerr>0) pxzsec[i][iaxis]->SetBinContent(ixy,iz,binerr);
	}
      }

  
  cout<<"Write output to file..."<<endl;
  fout->cd();
  hxyz->Write();
  hxy->Write();
  hxz->Write();
  hyz->Write();

  for(int iaxis=0; iaxis<3; iaxis++){
    pxyz[iaxis]->Write();
    for(int i=0; i<nsli; i++){
      psli[i][iaxis]->Write();
    }
    pxy[iaxis]->Write();
    for(int i=0; i<nsli; i++){
      pxysli[i][iaxis]->Write();
    }
    pxz[iaxis]->Write();
    for(int i=0; i<nsec; i++){
      pxzsec[i][iaxis]->Write();
    }
  }
  
  fout->Close();
  tree->Reset();
  return;
}
