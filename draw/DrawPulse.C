#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"

#define NOISE

TGraph *gr0;
TGraph *gr1;

void DrawPulse(){
  TCanvas *c = new TCanvas("c","c",700,400);
  c->SetMargin(0.08,0.01,0.12,0.01);
  int seg1,seg2;
  Double_t pos[3];
  Double_t spulse[37][121];
  Double_t core1[121];
  Double_t spulse1[36][121];
  Double_t core2[121];
  Double_t spulse2[36][121];

  TFile *f1 = new TFile("rootfiles/SimData0.root");
  TTree *tree1 = (TTree *)f1->Get("tree");
  tree1->SetBranchAddress("seg",&seg1);
  tree1->SetBranchAddress("pos",pos);
  tree1->SetBranchAddress("core",core1);
  tree1->SetBranchAddress("spulse",spulse1);

  TFile *f2 = new TFile("pulsedb/pulseA.root");
  TTree *tree2 = (TTree *)f2->Get("tree");
  tree2->SetBranchAddress("seg",&seg2);
  tree2->SetBranchAddress("core",core2);
  tree2->SetBranchAddress("spulse",spulse2);

  Float_t x[4477];
  Float_t y[4477];
  for(int iseg=0; iseg<37; iseg++)
    for(int i=0; i<121; i++)
      x[iseg*121+i] = iseg+i*1./121-0.5;

  int nentries = tree1->GetEntriesFast();
  int ientry;
  for(ientry = 398; ientry<nentries; ientry++){
    tree1->GetEntry(ientry);
    if(seg1==22){
      cout<<"entry1 "<<ientry<<" seg = "<<seg1
	  <<" pos = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<endl;
      break;
    }
  }

  for(int iseg=0; iseg<37; iseg++){
    for(int isig=0; isig<121; isig++){
      if(iseg==36) spulse[iseg][isig] = core1[isig];
      else         spulse[iseg][isig] = spulse1[iseg][isig];
    }
  }

#ifdef NOISE
  // add noise here
  double noise[4477],noise2[4477];
  double tmpnoise;
  for(int iseg=0; iseg<37; iseg++){
    tmpnoise = gRandom->Uniform(-0.01,0.01);
    for(int i=0; i<121; i++){
      tmpnoise += -0.2*(tmpnoise+gRandom->Uniform(-0.1,0.1));
      noise[iseg*121+i]=tmpnoise;
    }
  }

  // smooth noise
  int nsmooth=3;
  for(int i=0; i<4477; i++){
    noise2[i] = noise[i];
    if(i<nsmooth || i>=4477-nsmooth) continue;
    for(int j=1; j<=nsmooth; j++){
      noise2[i]+=noise[i+j]+noise[i-j];
    }
    noise2[i]=noise2[i]/(2*nsmooth+1);
  }

  for(int iseg=0; iseg<37; iseg++)
    for(int isig=0; isig<121; isig++)
      spulse[iseg][isig]+=noise2[isig+121*iseg];
#endif

  
  for(int iseg=0; iseg<37; iseg++)
    for(int isig=0; isig<121; isig++)
      y[iseg*121+isig] = spulse[iseg][isig];

  gr0 = new TGraph(4477,x,y);
  gr0->GetXaxis()->SetRangeUser(-0.5,36.5);
  gr0->GetYaxis()->SetRangeUser(-0.3,1.1);
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
    for(int isig=0; isig<121; isig++){
      if(iseg==36) spulse[iseg][isig] = core2[isig];
      else         spulse[iseg][isig] = spulse2[iseg][isig];
    }
  }
  
  for(int iseg=0; iseg<37; iseg++)
    for(int isig=0; isig<121; isig++)
      y[iseg*121+isig] = spulse[iseg][isig];

  gr1 = new TGraph(4477,x,y);
  gr1->SetLineColor(2);
  gr1->SetLineWidth(2);

  gr0->GetXaxis()->SetTitle("Segment Number");
  gr0->GetYaxis()->SetTitle("Normalised Charge");
  gr0->GetXaxis()->CenterTitle();
  gr0->GetYaxis()->CenterTitle();
  gr0->GetXaxis()->SetTitleSize(0.06);
  gr0->GetXaxis()->SetLabelSize(0.05);
  gr0->GetXaxis()->SetTitleOffset(0.9);
  gr0->GetYaxis()->SetTitleOffset(0.65);
  gr0->GetYaxis()->SetTitleSize(0.06);
  gr0->GetYaxis()->SetLabelSize(0.05);
  gr0->Draw("APL");
  //gr1->Draw("Lsame");
  
}
