{
TFile *fhist=new TFile("/home/dim2/urqmd_1036_115.mc.root");
TChain *chain = new TChain("mctree");
chain->Add("/home/dim2/urqmd_1036_1**.mctree.root");
TFile *f0=new TFile ("./ALLTREE.root","recreate");
f0->cd();
chain->Write();
f0->Close();
}
