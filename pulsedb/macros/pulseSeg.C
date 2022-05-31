#include "TCanvas.h"
#include "TRandom.h"
#include "iostream"
#include "TVirtualFFT.h"
#include "TH1.h"
#include "TGraph.h"
#include "time.h"

using namespace std;

const int NSegCore = 37;
const int NSig = 121;
const int npoint = 4477;

float ymin = -0.05;
float ymax = 1.1; 

void pulseSeg(){
  time_t t;
  gRandom->SetSeed(time(&t));
  
  const int nevts = 1;
  int seg = 20;
  
  TH1D *h = new TH1D("h","h",NSig-1,0,120);
  TFile *f = new TFile("pulse.root");
  TGraph *gr = (TGraph *)f->Get("Graph");
  
  int istart = NSig * seg;
  Double_t x[NSig], y[NSig];

  for(int isig=0; isig<NSig; isig++){
    x[isig] = isig;
    y[isig] = gr->GetPointY(istart+isig);
  }

  TGraph *gseg = new TGraph(NSig, x, y);
  gseg->SetNameTitle("gseg","");
  gseg->GetXaxis()->SetRangeUser(0,120);
  gseg->GetXaxis()->SetTitle("t / ns");
  gseg->GetYaxis()->SetTitle("amp.");
  //gseg->Draw("APL");

  
  for(int isig=0; isig<NSig; isig++)
    h->SetBinContent(isig,y[isig]);


  TCanvas *c = new TCanvas("c","c",600,900);
  c->Divide(1,3);
  c->cd(1);
  h->GetYaxis()->SetRangeUser(ymin,ymax);
  h->Draw("h");

  c->cd(2);
  TH1 *hm =0;
  TVirtualFFT::SetTransform(0);
  hm = h->FFT(hm, "MAG");
  hm->SetTitle("Magnitude of the 1st transform");
  hm->Draw();

  c->cd(3);
  TH1 *hp =0;
  hp = h->FFT(hp, "PH");
  hp->SetTitle("Phase of the 1st transform");
  hp->Draw();

  //return;
  const int NSig2 = 80;
  int nsig = NSig2;
  Double_t amp[NSig2];
  for(int i=0; i<NSig2; i++)
    amp[i] = y[10+i];
  //TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  TVirtualFFT *fft = TVirtualFFT::FFT(1, &nsig, "R2C M K");
  fft->SetPoints(amp);
  fft->Transform();
  Double_t *re_full = new Double_t[NSig];
  Double_t *im_full = new Double_t[NSig];
  fft->GetPointsComplex(re_full,im_full);
  for(int i=0; i<NSig2; i++){
    re_full[i] = re_full[i]/(NSig2);
    im_full[i] = im_full[i]/(NSig2);

    if(seg==20){
      if(i==0){
	re_full[i] = im_full[i] = 0;
      }else if(i>0 && i<4){
	re_full[i] = re_full[i] * gRandom->Gaus(1,0.2);
	im_full[i] = im_full[i] * gRandom->Gaus(1,0.2);
      }
      re_full[0] += -2*re_full[i];
      im_full[0] += -2*im_full[i];

    }else{
      if(i==0){
	re_full[i] = im_full[i] = 0;
      }else if(i>0 && i<10){
	re_full[i] = re_full[i] * gRandom->Gaus(1,0.5);
	im_full[i] = im_full[i] * gRandom->Gaus(1,0.5);
      }
      re_full[0] += -2*re_full[i];
      im_full[0] += -2*im_full[i];

    }
  }

  // filter ---------------------------
  /*
  for(int isig=0; isig<NSig; isig++){
    if(isig<20) continue;
    re_full[isig] = 0;
  }
  for(int isig=0; isig<NSig; isig++){
    if(isig<20) continue;
    im_full[isig] = 0;
  }

  //re_full[1] = im_full[1] = 0;
  //re_full[3] = im_full[3] = 0;
  //re_full[4] = im_full[4] = 0;
  //re_full[5] = im_full[5] = 0;
  */
  //-----------------------------------

  
  c->cd(2);
  TH1D *h2 = new TH1D("h2","re",NSig,0,NSig);
  for(int isig=0; isig<NSig; isig++) h2->SetBinContent(isig,re_full[isig]);
  h2->Draw();

  c->cd(3);
  TH1D *h3 = new TH1D("h3","im",NSig,0,NSig);
  for(int isig=0; isig<NSig; isig++) h3->SetBinContent(isig,im_full[isig]);
  h3->Draw();

  //return;
  
  // backward transform
  c->cd(3);
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &nsig, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1 *hb = 0;
  hb = TH1::TransformHisto(fft_back,hb,"Re");
  hb->SetTitle("The backward transform");
  TH1D *hb2 = new TH1D("hb2","The backward transform",NSig-1,0,120);
  for(int i=1; i<NSig; i++){
    hb2->SetBinContent(i,h->GetBinContent(i));
    if(i>8 && i<1+NSig2)
      hb2->SetBinContent(i,hb->GetBinContent(i-9));

    if(i>NSig2 && i<10+NSig2)
      hb2->SetBinContent(i,((i-NSig2)*h->GetBinContent(i)+(10+NSig2-i)*hb->GetBinContent(i-9))/10);
      
  }
  c->cd(3);
  hb2->SetLineColor(2);
  hb2->GetYaxis()->SetRangeUser(ymin,ymax);
  hb2->Draw("h");
  h->Draw("hsame");
  
  return;
}
