{
  TCanvas *c = new TCanvas("c","c",500,800);
  c->Divide(1,2);
  c->cd(1);
  tree->Draw("dist:NetCh>>h(100,0,1500,100,-1,11)","","colz");
  h->SetTitle("dist : NetCh");
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("colz");

  c->cd(2);
  TProfile *p = new TProfile("p","dist : NetCh",100,0,2000,0,100);
  p->GetXaxis()->SetRangeUser(0,1500);
  p->GetYaxis()->SetRangeUser(-1,11);
  p->GetXaxis()->SetLabelSize(0.05);
  p->GetYaxis()->SetLabelSize(0.05);
  tree->Draw("dist:NetCh>>p");

}
