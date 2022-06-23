const int npoints = 56;

void DrawBase(){
  TH1F *h = new TH1F("h","h",100,-50,50);
  float x[npoints], y[npoints];
  ifstream fin("NoiseBase.txt");
  for(int i=0; i<npoints; i++){
    x[i] = i;
    fin>>y[i];
    y[i] = y[i]/100;//<-------
    h->Fill(y[i]);
  }

  TGraph *g = new TGraph(npoints, x, y);

  g->Draw("APL");
  //h->Draw();
}
