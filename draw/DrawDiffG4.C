#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TVector3.h"
#include "TLegend.h"

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

TCanvas *c;

void DrawDiffG4(){
  TChain *tree = new TChain();
  for(int irun=0; irun<1; irun++){
    //tree->AddFile(Form("rootfiles/G4Ana1n/G4AnaData%04d.root",irun),0,"tree");
    //tree->AddFile(Form("rootfiles/G4Ana/G4AnaData%04d_n90keV.root",irun),0,"tree");
    tree->AddFile(Form("rootfiles/G4Ana/det2/G4AnaData%04d.root",irun),0,"tree");
    //tree->AddFile(Form("rootfiles/G4Anadev/G4AnaData%04d_0_1.root",irun),0,"tree");
    //tree->AddFile(Form("rootfiles/G4Anadev/G4AnaData%04d_0_10.root",irun),0,"tree");
    //tree->AddFile(Form("rootfiles/G4Anadev/G4AnaData%04d_0_30.root",irun),0,"tree");
    //tree->AddFile(Form("rootfiles/G4Anadev/G4AnaData%04d_1.root",irun),0,"tree");
    //tree->AddFile(Form("rootfiles/G4Anadev/G4AnaData%04d_10.root",irun),0,"tree");
    //tree->AddFile(Form("rootfiles/G4Anadev/G4AnaData%04d_30.root",irun),0,"tree");
  }
  
  Int_t seg;
  Float_t NetCh;
  Float_t simpos[3];
  Float_t anapos[3];
  Float_t dist;
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("NetCh",&NetCh);
  tree->SetBranchAddress("simpos",simpos);
  tree->SetBranchAddress("anapos",anapos);
  tree->SetBranchAddress("dist",&dist);

  int nentries = tree->GetEntriesFast();
  cout<<"read "<<nentries<<" events from tree"<<endl;

  // define hist
  TH1D *hdist = new TH1D("hdist","dist: anapos - simpos",100,0,10);
  TH1D *hdiff[3];
  hdiff[0] = new TH1D("hdiff0","rdiffphi",100,-10,10);
  hdiff[1] = new TH1D("hdiff1","diffr",100,-10,10);
  hdiff[2] = new TH1D("hdiff2","diffz",100,-10,10);
  
  int ievt;
  for(ievt=0; ievt<nentries; ievt++){
    if(ievt%1000==0)
      cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<flush;
    tree->GetEntry(ievt);
    if(NetCh<900) continue;

    TVector3 vecsim(simpos[0],simpos[1],0);
    TVector3 vecana(anapos[0],anapos[1],0);
    double simphi = vecsim.Phi();
    double simr = vecsim.Mag();
    double anaphi = vecana.Phi();
    double anar = vecana.Mag();
    double diffphi = anaphi-simphi;
    if(diffphi>TMath::Pi()) diffphi-= 2*TMath::Pi();
    if(diffphi<-TMath::Pi()) diffphi+= 2*TMath::Pi();
    
    hdist->Fill(dist);
    hdiff[0]->Fill(simr*diffphi);
    hdiff[1]->Fill(anar-simr);
    hdiff[2]->Fill(anapos[2]-simpos[2]);
  }
  cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<endl;

  c = new TCanvas("c","c",400,300);

  hdiff[0]->SetTitle("");
  hdiff[0]->SetStats(0);
  hdiff[0]->SetLineColor(1);
  hdiff[1]->SetLineColor(2);
  hdiff[2]->SetLineColor(4);

  hdiff[0]->GetYaxis()->SetRangeUser(0,0.45e5);
  hdiff[0]->GetYaxis()->SetMaxDigits(4);
  hdiff[0]->GetXaxis()->SetLabelSize(0.05);
  hdiff[0]->GetYaxis()->SetLabelSize(0.05);
  hdiff[0]->Draw();
  hdiff[1]->Draw("same");
  hdiff[2]->Draw("same");

  TLegend *label = new TLegend(0.65,0.65,0.9,0.9);
  label->AddEntry(hdiff[0],"diff in Phi","l");
  label->AddEntry(hdiff[1],"diff in R","l");
  label->AddEntry(hdiff[2],"diff in Z","l");
  label->Draw("same");

  return;
}
