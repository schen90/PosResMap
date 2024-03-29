#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"

TGraph *gr;
TGraph *gr2;

void CompPulse(int ientry = 26915){
  //---------------------------
  Int_t simseg = 0;
  Double_t simpos[3] = {0,0,0};
  //---------------------------

  TCanvas *c = new TCanvas("c","c",1265,400);
  c->SetMargin(0.06,0.01,0.12,0.01);

  // simulation
  TFile *f = new TFile("pulseA.root");
  TTree *tree = (TTree *)f->Get("tree");
  Int_t seg;
  Double_t pos[3];
  Double_t core[121];
  Double_t spulse[4356];
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("pos",pos);
  tree->SetBranchAddress("core",core);
  tree->SetBranchAddress("spulse",spulse);
  int nentries = tree->GetEntriesFast();

  // with Dev
  TFile *f2 = new TFile("Dev/pulseA.root");
  TTree *tree2 = (TTree *)f2->Get("tree");
  Int_t seg2;
  Double_t pos2[3];
  Double_t core2[121];
  Double_t spulse2[4356];
  tree2->SetBranchAddress("seg",&seg2);
  tree2->SetBranchAddress("pos",pos2);
  tree2->SetBranchAddress("core",core2);
  tree2->SetBranchAddress("spulse",spulse2);
  
 
  tree->GetEntry(ientry);
  tree2->GetEntry(ientry);
  
  Double_t x[4477];
  for(int iseg=0; iseg<37; iseg++)
    for(int i=0; i<121; i++)
      x[iseg*121+i] = iseg+i*1./121-0.5;

  double tmpspulse[4477];
  for(int i=0; i<4356; i++) tmpspulse[i] = spulse[i];
  for(int i=0; i<121; i++) tmpspulse[i+4356] = core[i];
  
  double tmpspulse2[4477];
  for(int i=0; i<4356; i++) tmpspulse2[i] = spulse2[i];
  for(int i=0; i<121; i++) tmpspulse2[i+4356] = core2[i];

  gr = new TGraph(4477,x,tmpspulse);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitleOffset(0.9);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetTitle("Segment Number");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleOffset(0.5);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetTitle("Normalised Charge");
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->SetRangeUser(-0.5,36.5);
  gr->GetYaxis()->SetRangeUser(-0.1,1.1);

  gr2 = new TGraph(4477,x,tmpspulse2);
  gr2->SetLineColor(2);

  cout<<Form("sim seg: %d  pos: %f  %f  %f",seg,pos[0],pos[1],pos[2])<<endl;
  cout<<Form("dev seg: %d  pos: %f  %f  %f",seg2,pos2[0],pos2[1],pos2[2])<<endl;
  
  gr->Draw("APL");
  gr2->Draw("Lsame");
  
  f->Close();
  return;
}
