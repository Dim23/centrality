#define NewTree_cxx
#include "NewTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void NewTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L NewTree.C
//      root> NewTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

    TH1F *Hpart=new TH1F("Hpart","Hpart",500,0,500); 
    TH1F *Hnh=new TH1F("Hh","Hh",2000,0,2000); 
    TH1F *Hbimp=new TH1F("Hbimp","Hbimp",100,0,20); 

    double p,pt,eta;
    int s; 

TFile *f0=new TFile ("./HIST.root","recreate");
f0->cd();

TTree *mytree = new TTree("mytree","My Tree"); 
int snh;
float snbimp;
mytree->Branch("bimp", &snbimp, "bimp/F");
mytree->Branch("nh", &snh, "nh/I");

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        s=0;

        for (int i=0;i<nh;i++)
        {
            pt=sqrt(momx[i]*momx[i]+momy[i]*momy[i]);
            p=sqrt(momx[i]*momx[i]+momy[i]*momy[i]+momz[i]*momz[i]);
            eta=log((p+momz[i])/(p-momz[i]));

            if (pt>0.15 && charge[i]!=0 && eta*eta<1 ) 
            {   
                s+=1;
            };
        }

        snbimp=bimp;
        snh=s;
        mytree->Fill();
         Hpart->Fill(npart);
         Hnh->Fill(s);
         Hbimp->Fill(bimp);

// if (Cut(ientry) < 0) continue;
    }


Hpart->Write();
Hnh->Write();
Hbimp->Write();
mytree->Write();
f0->Close();

}
