


#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>

#include <TRatioPlot.h>
using namespace std;
//gROOT->ForceStyle();
//gROOT->SetStyle("Pub");


//Задаем постоянные;
double sigma=685,pi=TMath::Pi();

double stepx=0.05; int Nx=(1/stepx),Nn=450;
//Уссловная вероятность;
double f_k(double x,double p,double q,double r,double s,double nmax) 
{
    return nmax*exp(p*x+q*pow(x,2)+r*pow(x,3)+s*pow(x,4));
};

double Pnb(double nn, double x,double p,double q,double r,double s,double nmax,double teta) 
{
    double kk=f_k(x,p,q,r,s,nmax)/teta;
    return ROOT::Math::gamma_pdf(nn,kk,teta);
};
//Функция для фита;
double Pn(double nn,double p,double q,double r,double s,double nmax,double teta) 
{
    double Integr=0.0,x=0.0;
    for (int n=0;n<Nx;n++) 
    {
        Integr+=Pnb(nn,x,p,q,r,s,nmax,teta);
        x+=stepx;
    }
    return (Integr-0.5*(Pnb(nn,0,p,q,r,s,nmax,teta)+Pnb(nn,1,p,q,r,s,nmax,teta)))*stepx;
};
double FitF(double *w, double *par) 
{
    double p=par[0],q=par[1],r=par[2],s=par[3],nmax=par[4],teta=par[5];
    return Pn(w[0],p,q,r,s,nmax,teta);
};
//
double cb(double b) 
{
    return pi*b*b/sigma;
};


vector <double> NormVect(vector <double> a,double stepb){
int N=a.size();
double I=0;
for(int i=0;i<N;i++){
I+=a[i];
}
I=(I-0.5*(a[0]+a[N-1]))*stepb;
for(int i=0;i<N;i++){
a[i]=a[i]/I;
}
return a;
}

class Fit {
    public:
      
double p,q,r,s,nmax,teta;

//Считываем данные с файла;
    TFile *file=new TFile("/home/dim2/HIST.root");
    TTree *t=(TTree*)file->Get("mytree");

   Int_t           npart;
   Float_t         bimp;
    TH1F *Gev=new TH1F("N_{ch}"," ",500,0,500);


void Start()
{
    //Задаем функцию для фитирования
    TF1 *fc=new TF1("fc",FitF,0,Nn,6);
    

    //Заполняем граф и создаем канву для графика;
    t->SetBranchAddress("bimp", &bimp);
    t->SetBranchAddress("nh", &npart);

    Gev->Sumw2();
    Int_t nentries = (Int_t)t->GetEntries();

    for (Int_t i=0; i<nentries; i++) 
    {
        t->GetEntry(i); Gev->Fill(npart);
    }
    Gev->Scale(1.0/Gev->GetEntries());
    TH1F *GevC=(TH1F*)Gev->Clone();

            // Начальные значения параметров;
            fc->SetParameter(0,-5); fc->SetParameter(1,3);
            fc->SetParameter(2,-6); fc->SetParameter(3,2.35);
            fc->SetParameter(4,400); fc->SetParameter(5,2); 

             
   TCanvas *c=new TCanvas("c","Graph Draw Options",800,600);
// Upper plot will be in pad1
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
   pad1->SetBottomMargin(0.05);
  // Upper and lower plot are joined
         // Vertical grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();       
        pad1->cd()->SetLogy();
        Gev->GetYaxis()->SetRangeUser(1e-6,0.5);
        Gev->GetXaxis()->SetRangeUser(0,Nn);
        Gev->Fit(fc,"RM");
        Gev->Draw();

         p=(fc->GetParameter(0));q=(fc->GetParameter(1));r=(fc->GetParameter(2));
         s=(fc->GetParameter(3));nmax=(fc->GetParameter(4));teta=(fc->GetParameter(5));
         

   
   c->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
   pad2->SetBottomMargin(0.3);
  pad2->SetTopMargin(0.01);
   // vertical grid
   pad2->Draw();
   pad2->cd();       // pad2 becomes the current pad
   // Define the ratio plot

TF1 *NEWfc=new TF1;
NEWfc=Gev->GetFunction("fc");

   GevC->Divide(NEWfc);
   GevC->SetMinimum(0.4);  // Define Y ..
   GevC->SetMaximum(1.6); // .. range
   GevC->Sumw2();
   GevC->SetStats(0);
   Gev->SetStats(0);      // No statistics on lower plot
  
// Y axis ratio plot settings
    Gev->GetYaxis()->SetTitle("1/N dN_{ch}/dN");
    Gev->GetYaxis()->SetNdivisions(505);
    Gev->GetYaxis()->SetTitleSize(20);
    Gev->GetYaxis()->SetTitleFont(43);
    Gev->GetYaxis()->SetTitleOffset(1.1);
    Gev->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Gev->GetYaxis()->SetLabelSize(15);
// X axis ratio plot settings
    Gev->GetXaxis()->SetRangeUser(0,Nn);
    Gev->GetXaxis()->SetTitleSize(20);
    Gev->GetXaxis()->SetTitleFont(43);
    Gev->GetXaxis()->SetTitleOffset(1.5);
    Gev->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Gev->GetXaxis()->SetLabelSize(15);


// Y axis ratio plot settings
    GevC->GetYaxis()->SetTitle("Fit/Data");
    GevC->GetYaxis()->SetNdivisions(505);
    GevC->GetYaxis()->SetTitleSize(20);
    GevC->GetYaxis()->SetTitleFont(43);
    GevC->GetYaxis()->SetTitleOffset(1.1);
    GevC->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    GevC->GetYaxis()->SetLabelSize(15);
// X axis ratio plot settings
    GevC->GetXaxis()->SetRangeUser(0,Nn);
    GevC->GetXaxis()->SetTitle("N_{ch}");
    GevC->GetXaxis()->SetTitleSize(20);
    GevC->GetXaxis()->SetTitleFont(43);
    GevC->GetXaxis()->SetTitleOffset(3);
    GevC->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    GevC->GetXaxis()->SetLabelSize(15);

    GevC->Draw();
TLine *line=new TLine(0,1,Nn,1);
line->Draw("SAME");

        };


 //Метод для нахождения интеграла P(n);
        double IntPn(double n0,double nn,int M)
        {
            double Integr=0.0,n=n0,step=(nn-n0)/M;
            for (int k=0;k<=M;k++)
            {
                Integr+=Pn(n,p,q,r,s,nmax,teta); 
                n+=step;       
            }
        return (Integr-0.5*(Pn(n0,p,q,r,s,nmax,teta)+Pn(nn,p,q,r,s,nmax,teta)))*step;
        };

//Метод для нахождения центральности;
        double c(double n0,int M)
        {
            double I=0.0;
            I=IntPn(n0,Nn,M);
        return I;
        };

//Метод для нахождения n по заданной центральности;
        double cFind(double c0,double n0) {
            double n=n0,cen=0,stepn=1,step=1;
            int M=round((Nn-n0)/step);
            while (cen<c0)
            {
                n=n-stepn; cen=c(n,M); 
            }
        return n;}

 //Метод для нахождения интеграла P(n\b);
        double IntgPnb(double n0,double nn,int M,double b)
        {
            double Integr=0.0,n=n0,step=(nn-n0)/M;
            for (int k=0;k<=M;k++)
            {
                Integr+=Pnb(n,cb(b),p,q,r,s,nmax,teta);
                 n+=step; 
            }
        return (Integr-0.5*(Pnb(n0,cb(b),p,q,r,s,nmax,teta)+Pnb(nn,cb(b),p,q,r,s,nmax,teta)))*step;
   
        };


//Метод для нахождения распределения b для промежутка n0<n<nn;
        double IntPb(double n0,double nn,int M,double b)
        {
            double I1,I2;
            I1=IntgPnb(n0,nn,M,b);
            I2=IntPn(n0,nn,M);
        return 2*pi*b*I1/(sigma*I2);
        } ;




//Запись и вывод полученного распределения;
       vector <double> PlotIntb(double n0,double nn,int M,double bn,int N)
        {
            vector<double> bi,Pbi;
            bi.clear();
            Pbi.clear();
            double stepb=bn/N;
            for (int n=0;n<=N;n++)
            {
                bi.push_back(stepb*n);
                Pbi.push_back(IntPb(n0,nn,M,bi[n]));
            }
            return NormVect(Pbi,stepb); 
         };


TH1F* Histb(double n0,double nn,double bn,int N)
{
TH1F *hbH=new TH1F("hbH","impact",N,0,bn); 
double k=N/bn;
int nentries = (int)t->GetEntries();
t->SetBranchAddress("bimp", &bimp);
for (int i=0; i<nentries; i++) 
{
    t->GetEntry(i);
    if (npart <= nn && npart >= n0) 
    {  
        hbH->Fill(bimp);
     } 

}
        hbH->Scale(k/hbH->GetEntries());
        return hbH;
}


void Plotbn(double n0,double nn,int M,double bn,int N)
        {
         vector<double> b,Nc;
         vector<double> d1;
         d1.clear();
         
         double stepb=bn/N;
         for (int n=0;n<=N;n++)
            {
                b.push_back(stepb*n);
            }

        d1=PlotIntb(n0,nn,M,bn,N);
        TGraph* grb1=new TGraph(N,b.data(),d1.data());grb1->SetName("grBimp_1");

        grb1->GetYaxis()->SetRangeUser(0,1);
        grb1->GetXaxis()->SetRangeUser(0,bn);
        
        TCanvas *c3=new TCanvas("c3","Graph Draw Options",10,0,800,800);
        c3->cd();
        grb1->SetLineColor(1);grb1->Draw(); 

        TH1F *hb=new TH1F("hb","impact",N,0,bn); 
        hb=Histb(n0,nn,bn,N);
        hb->Draw("SAME");

}


double meanb(vector<double> d,vector<double> b,double stepb,int N){
double mb=0,id=0;
 for (int n=0;n<=N;n++)
 {
  mb+=d[n]*b[n];
  id+=d[n];
 }
return mb/id;
 }


double RMSb(vector<double> d,vector<double> b,double stepb,int N){

  double mb=0,id=0,bm=meanb(d,b,stepb,N);
 for (int n=0;n<=N;n++)
 {
  mb+=d[n]*pow((b[n]-bm),2);
    id+=d[n];
  }
return sqrt(mb/id);
 }

      
void PlotAllb(int M,double bn,int N)
        {
         vector<double> b,Nc(11);
         vector<double> d[10];
         
         b.clear();
         
         double stepb=bn/N;
         for (int n=0;n<=N;n++)
            {
                b.push_back(stepb*n);
            }
          
        Nc[0]=Nn;float c0=0.1;Nc[10]=0;
        for (int n=1;n<10;n++)
            {   
                Nc[n]=cFind(c0,Nc[n-1]);
                c0+=0.1;
            }

        for (int n=0;n<10;n++)
            {   
                d[n]=PlotIntb(Nc[n+1],Nc[n],M,bn,N);
            }
    
    
double x[10],ex[10],x1[10],ex1[10],y[10],ey[10];       
TH1F *hb[10];
        for (int n=0;n<10;n++)
            {   
                hb[n]=Histb(Nc[n+1],Nc[n],bn,N);
                x[n]=meanb(d[n],b,stepb,N);
                ex[n]=RMSb(d[n],b,stepb,N);
                x1[n]=hb[n]->GetMean();
                ex1[n]=hb[n]->GetRMS();
                y[n]=(2*n+1)*5;ey[n]=0;
            }

        
TCanvas *c3=new TCanvas("c3","Graph Draw Options",10,0,800,800); c3->cd();

for (int k=0;k<10;k++) {
std::cout <<x[k]<<" " <<ex[k] <<" Fit "<<x1[k]<<" " <<ex1[k]<< " Model "<< y[k]<< " nn "<<Nc[k]<<std::endl;
                            }

TGraphErrors *gr = new TGraphErrors(10,y,x,ey,ex);
TGraphErrors *gr1 = new TGraphErrors(10,y,x1,ey,ex1);
        
gr->GetYaxis()->SetRangeUser(0,bn);
gr->GetXaxis()->SetRangeUser(0,100);
gr->SetLineColor(1);
gr1->SetLineColor(2);
gr->SetTitle(" Reconstructed ");
gr1->SetTitle(" UrQMD ");
gr->SetName("FitImpact");
gr1->SetName("ModalImpact");

gr->GetYaxis()->SetTitle("<b>");
gr->GetXaxis()->SetTitle("Centrality,%");

gr->SetMarkerStyle(21);
gr1->SetMarkerStyle(20);
gr->SetMarkerSize(1.6);
gr1->SetMarkerSize(1.6);
gr->SetMarkerColor(2);
gr1->SetMarkerColor(1);

gr->Draw("AP");
gr1->Draw("P SAME"); 
gPad->BuildLegend();


};
 
//Метод для нахождения распределения прецельного параметра;
        double Pb(double n0,double b)
        {
            double I1,I2;
            
            I1=Pnb(n0,cb(b),p,q,r,s,nmax,teta);
            I2=Pn(n0,p,q,r,s,nmax,teta);
            
        return 2*pi*b*I1/(sigma*I2);
   
        };

//Запись и вывод полученного распределения;
        void Plotb(double n0,double bn,int N)
        {
            double bi[N+1],Pbi[N+1],stepb=bn/N;

            for (int n=0;n<=N;n++)
            {
                bi[n]=stepb*n;
                Pbi[n]=Pb(n0,bi[n]);
            }
            TGraph* grb=new TGraph(N,bi,Pbi);
            grb->SetTitle(" Probability distribution of impact parametr");
            TCanvas *c4=new TCanvas("c4","Graph Draw Options",10,0,400,400);
            c4->cd();
            grb->Draw("AL"); 
         };

};
