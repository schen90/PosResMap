#include "TStyle.h"
#include "TList.h"
#include "TRandom.h"
#include "TH2.h"
#include "TProfile3D.h"

void Plot2D(){
  gStyle->SetOptStat(0);
  // get hist 2
  //TFile *f = new TFile("plotchain_db0.root");
  //TFile *f = new TFile("plotchain_coarse.root");
  //TFile *f = new TFile("plotchain_cf.root");
  //TFile *f = new TFile("plotchain_full.root");
  TFile *f = new TFile("plotchain_n90keV.root");

  TH2D *h[6];
  h[0] = (TH2D *)f->Get("pxy0");
  h[1] = (TH2D *)f->Get("pxy1");
  h[2] = (TH2D *)f->Get("pxy2");
  h[3] = (TH2D *)f->Get("pxz0");
  h[4] = (TH2D *)f->Get("pxz1");
  h[5] = (TH2D *)f->Get("pxz2");
  
  TCanvas *c = new TCanvas("c","abc",900,600);
  c->Divide(3,2);

  for(int i=0; i<6; i++){
    c->cd(i+1);
    //h[i]->GetZaxis()->SetRangeUser(0,3);
    h[i]->Draw("colz");
  }

  
}
