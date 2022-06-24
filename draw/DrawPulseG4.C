#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TInterpreter.h"

#include <vector>

TGraph *gr0;
TGraph *gr1;

void DrawPulseG4(){
  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");
  gInterpreter->GenerateDictionary("vector<vector<float>>","vector");

  TCanvas *c = new TCanvas("c","c",600,400);
  c->SetMargin(0.12,0.01,0.12,0.01);
  int seg1,seg2;
  Double_t spulse[37][56];
  Double_t core1[56];
  Double_t spulse1[36][56];
  Float_t spulse2[37][56];

  vector<vector<int>>   *pseg   = 0;
  vector<vector<float>> *vcore   = 0;
  vector<vector<float>> *vspulse = 0;

  TFile *f1 = new TFile("rootfiles/G4Sim/G4SimData0000.root");
  TTree *tree1 = (TTree *)f1->Get("tree");
  tree1->SetBranchAddress("pseg",&pseg);
  tree1->SetBranchAddress("core",&vcore);
  tree1->SetBranchAddress("spulse",&vspulse);

  TFile *f2 = new TFile("pulsedb/LibTrap_A001.root");
  TTree *tree2 = (TTree *)f2->Get("tree");
  tree2->SetBranchAddress("seg",&seg2);
  tree2->SetBranchAddress("spulse",spulse2);

  Float_t x[2072];
  Float_t y[2072];
  for(int iseg=0; iseg<37; iseg++)
    for(int i=0; i<56; i++)
      x[iseg*56+i] = iseg+i*1./56-0.5;

  int nentries = tree1->GetEntriesFast();
  int ientry;
  bool kfound = false;
  for(ientry = 0; ientry<nentries; ientry++){
    tree1->GetEntry(ientry);

    for(int i=0; i<pseg->at(0).size(); i++){
      if(pseg->at(0)[i]==23){
	cout<<"entry1 "<<ientry<<" seg"<<i<<" = "<<pseg->at(0)[i]<<endl;
	kfound = true;
	break;
      }
    }
    if(kfound){
      for(int iseg=0; iseg<36; iseg++){
	for(int isig=0; isig<56; isig++){
	  spulse1[iseg][isig] = vspulse->at(0)[iseg*56+isig];
	}
      }
      for(int isig=0; isig<56; isig++){
	core1[isig] = vcore->at(0)[isig];
      }
	
      break;
    }
  }

  for(int iseg=0; iseg<37; iseg++){
    for(int isig=0; isig<56; isig++){
      if(iseg==36) spulse[iseg][isig] = core1[isig];
      else         spulse[iseg][isig] = spulse1[iseg][isig];
    }
  }
  
  for(int iseg=0; iseg<37; iseg++)
    for(int isig=0; isig<56; isig++)
      y[iseg*56+isig] = spulse[iseg][isig];
  
  gr0 = new TGraph(2072,x,y);
  gr0->GetXaxis()->SetRangeUser(-0.5,36.5);
  //gr0->GetYaxis()->SetRangeUser(-0.3,1.1);
  gr0->SetLineWidth(2);
  gr0->SetTitle("");
  
  nentries = tree2->GetEntriesFast();
  for(ientry = 31300; ientry<nentries; ientry++){
    tree2->GetEntry(ientry);
    if(seg2==21){
      cout<<"entry2 "<<ientry<<" seg = "<<seg2<<endl;
      break;
    }
  }


  for(int iseg=0; iseg<37; iseg++){
    for(int isig=0; isig<56; isig++){
      spulse[iseg][isig] = spulse2[iseg][isig];
    }
  }
  
  for(int iseg=0; iseg<37; iseg++)
    for(int isig=0; isig<56; isig++)
      y[iseg*56+isig] = spulse[iseg][isig];

  gr1 = new TGraph(2072,x,y);
  gr1->SetLineColor(2);
  gr1->SetLineWidth(2);

  
  gr0->GetXaxis()->SetTitle("Segment Number");
  gr0->GetYaxis()->SetTitle("Charge");
  gr0->GetXaxis()->CenterTitle();
  gr0->GetYaxis()->CenterTitle();
  gr0->GetXaxis()->SetTitleSize(0.06);
  gr0->GetXaxis()->SetLabelSize(0.05);
  gr0->GetXaxis()->SetTitleOffset(0.9);
  gr0->GetYaxis()->SetTitleOffset(1);
  gr0->GetYaxis()->SetTitleSize(0.06);
  gr0->GetYaxis()->SetLabelSize(0.05);
  gr0->Draw("APL");
  //gr1->Draw("Lsame");
  
}
