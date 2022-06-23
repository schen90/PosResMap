#include "TCanvas.h"
#include "TFile.h"
#include "TH3.h"
#include "TH2.h"

using namespace std;

TCanvas *c = new TCanvas("c","c",900,300);

void CompareRes(){
  //TFile *f1 = new TFile("det0_Noise/plotresol.root");
  TFile *f1 = new TFile("det0_Dev/plotresol.root");
  TH3D *h1[3];
  for(int i=0; i<3; i++) h1[i] = (TH3D *)f1->Get(Form("pxyz%d",i));  

  //TFile *f2 = new TFile("det0_Noise/plotboot.root");
  TFile *f2 = new TFile("det0_Dev/plotboot.root");
  TH3D *h2[3];
  for(int i=0; i<3; i++) h2[i] = (TH3D *)f2->Get(Form("pxyz%d",i));  

  TH2D *h[3];
  h[0] = new TH2D("h0","boot resolution : resolution [phi]",100,0,12,100,0,12);
  h[1] = new TH2D("h1","boot resolution : resolution [r]",100,0,3,100,0,3);
  h[2] = new TH2D("h2","boot resolution : resolution [z]",100,0,8,100,0,8);
  
  int Nx = h1[0]->GetNbinsX();
  int Ny = h1[0]->GetNbinsY();
  int Nz = h1[0]->GetNbinsZ();

  for(int ix=0; ix<Nx; ix++){
    for(int iy=0; iy<Ny; iy++){
      for(int iz=0; iz<Nz; iz++){
	for(int dir=0; dir<3; dir++){
	  double res1 = h1[dir]->GetBinContent(ix,iy,iz);
	  double res2 = h2[dir]->GetBinContent(ix,iy,iz);
	  if(res1>0 && res2>0){
	    h[dir]->Fill(res1,res2);
	  }
	}
      }
    }
  }

  c->Divide(3,1);
  c->cd(1);
  h[0]->Draw("colz");
  c->cd(2);
  h[1]->Draw("colz");
  c->cd(3);
  h[2]->Draw("colz");
  
  /*
  TFile *f = new TFile("CompareRes.root","RECREATE");
  f->cd();
  for(int dir=0; dir<3; dir++) h[dir]->Write();
  f->Close();
  */
  return;
}
