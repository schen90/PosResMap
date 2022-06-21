{
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","pos",300,600);
  c->Divide(1,2);

  c->cd(1);
  XYrels->Draw("col");
  c->cd(2);
  XZrels->Draw("col");

}
