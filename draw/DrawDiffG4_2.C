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

void DrawDiffG4_2(){
  TChain *tree = new TChain();
  for(int irun=0; irun<1; irun++){
    tree->AddFile(Form("rootfiles/G4Ana1n/G4AnaData%04d.root",irun),0,"tree");
  }

  Int_t numNetCharges;
  Int_t nhits;
  Int_t seg;
  Float_t simpos[3];
  Float_t anapos[3];
  Float_t dist;
  tree->SetBranchAddress("numNetCharges",&numNetCharges);
  tree->SetBranchAddress("nhits",&nhits);
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("simpos",simpos);
  tree->SetBranchAddress("anapos",anapos);
  tree->SetBranchAddress("dist",&dist);

  int nentries = tree->GetEntriesFast();
  cout<<"read "<<nentries<<" events from tree"<<endl;

  // define hist
  TH1D *hdist[4];
  hdist[0] = new TH1D("hdist0","netQ=1 && nhits=1",100,0,10);
  hdist[1] = new TH1D("hdist1","netQ>1 && nhits=1",100,0,10);
  hdist[2] = new TH1D("hdist2","netQ=1 && nhits>1",100,0,10);
  hdist[3] = new TH1D("hdist3","netQ>1 && nhits>1",100,0,10);
  
  int ievt;
  for(ievt=0; ievt<nentries; ievt++){
    if(ievt%1000==0)
      cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<flush;
    tree->GetEntry(ievt);

    if(numNetCharges==1 && nhits==1)  hdist[0]->Fill(dist);
    if(numNetCharges>1  && nhits==1)  hdist[1]->Fill(dist);
    if(numNetCharges==1 && nhits>1)   hdist[2]->Fill(dist);
    if(numNetCharges>1  && nhits>1)   hdist[3]->Fill(dist);
  }
  cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<endl;

  double max[4], MAX=0;
  for(int i=0; i<4; i++){
    max[i] = hdist[i]->GetMaximum();
    if(max[i]>MAX) MAX=max[i];
  }

  for(int i=0; i<4; i++){
    hdist[i]->Scale(MAX/max[i]);
  }
  
  
  c = new TCanvas("c","c",400,300);

  hdist[0]->SetTitle("");
  hdist[0]->SetStats(0);
  hdist[0]->SetLineColor(1);
  hdist[1]->SetLineColor(2);
  hdist[2]->SetLineColor(4);
  hdist[3]->SetLineColor(6);

  hdist[0]->GetYaxis()->SetRangeUser(0,1.6e4);
  hdist[0]->GetYaxis()->SetMaxDigits(4);
  hdist[0]->GetXaxis()->SetLabelSize(0.05);
  hdist[0]->GetYaxis()->SetLabelSize(0.05);
  hdist[0]->Draw("h");
  hdist[1]->Draw("hsame");
  hdist[2]->Draw("hsame");
  hdist[3]->Draw("hsame");

  TLegend *label = new TLegend(0.45,0.55,0.9,0.9);
  label->AddEntry(hdist[0],"netQ=1 && nhits=1","l");
  label->AddEntry(hdist[1],"netQ>1 && nhits=1","l");
  label->AddEntry(hdist[2],"netQ=1 && nhits>1","l");
  label->AddEntry(hdist[3],"netQ>1 && nhits>1","l");
  label->Draw("same");

  return;
}
