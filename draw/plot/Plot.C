#include "TStyle.h"
#include "TList.h"
#include "TRandom.h"
#include "TH3.h"
#include "TProfile3D.h"

void Plot(){
  // get hist 3
  //TFile *f = new TFile("plotdata.root");
  TFile *f = new TFile("plotchain.root");
  TH3 *h3 = (TH3F *)f->Get("pxyz0");
  //TH3 *h3 = (TH3F *)f->Get("psegr3_0");

  TList *lf = h3->GetListOfFunctions();
  if (lf){
    TF1 *tf = new TF1("TransferFunction","x*x/100",0,10.);
    lf->Add(tf);
  }
  gStyle->SetCanvasPreferGL(1);
  gStyle->SetPalette(1);
  TCanvas *c = new TCanvas("c","abc",700,700);
  h3->Draw("glcolz");

  
}
