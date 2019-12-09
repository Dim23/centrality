


#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <TRatioPlot.h>

using namespace std;

//Задаем постоянные;
double sigma=685,pi=TMath::Pi(),bmax=14.76;
double stepx=0.001; int Nx=(1/stepx),Nn=450,binN=250;

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


//Функция для распределения n при cb=xn;
vector<double> Pnx(double p,double q,double r,double s,double nmax,double teta,double x0,double xn,int Nc) 
{
    double Integr,x,stepc=(xn-x0)/Nc;
    vector<double> VPn;
   int N=500;
    for (int n=0;n<=N;n++)
    {
       
        Integr=0.0;x=x0;
        for (int i=0;i<=Nc;i++) {
        if(i%2==0)
                    {
                        Integr+=2*Pnb(n,x,p,q,r,s,nmax,teta);
                    }
                else
                    {
                        Integr+=4*Pnb(n,x,p,q,r,s,nmax,teta);
                    }
                x+=stepc;       
            }
            Integr=(Integr-(Pnb(n,x0,p,q,r,s,nmax,teta)+Pnb(n,xn,p,q,r,s,nmax,teta)))*stepx/3;
            VPn.push_back(Integr);
        };
        
        
    return VPn;
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



vector <double> NormVect(vector <double> a,double stepb)
{ 
int M=a.size();double Integ=0;
for (int k=0;k<M;k++)
            {
               /* if(k%2==0)
                    {
                        Integ+=2*a[k];
                    }
                else
                    {
                        Integ+=4*a[k];
                    }  */
                Integ+=a[k];
            }
Integ=Integ*stepb;

for(int i=0;i<M;i++){
a[i]=a[i]/Integ;
}
            return a; 
}


double meanb(vector<double> d,vector<double> b,double stepb){
double mb=0,id=0;
 for (int n=0;n<d.size();n++)
 {
  mb+=d[n]*b[n];
  id+=d[n];
 }
return mb/id;
 }


double RMSb(vector<double> d,vector<double> b,double stepb){

  double mb=0,id=0,bm=meanb(d,b,stepb);
 for (int n=0;n<d.size();n++)
 {
  mb+=d[n]*pow((b[n]-bm),2);
    id+=d[n];
  }
return sqrt(mb/id);
 }



class Fit {
    public:
      
double p,q,r,s,nmax,teta;
   Int_t           npart;
   Float_t         bimp;
//Считываем данные с файла;
    TFile *file=new TFile("/home/dim2/HIST.root");
    TTree *t=(TTree*)file->Get("mytree");
    TH1F *Gev=new TH1F("N_{ch}"," ",500,0,500);
    
    
void Start()
{
    //Задаем функцию для фитирования
    TF1 *fc=new TF1("fc",FitF,0,450,6);
    t->SetBranchAddress("nh", &npart);
    t->SetBranchAddress("bimp", &bimp);
    Gev->Sumw2();
    int nentries = (int)t->GetEntries();

    for (Int_t i=0; i<nentries; i++) 
    {
        t->GetEntry(i); Gev->Fill(npart);
    }
    Gev->Scale(1/Gev->GetEntries());
    TH1F *GevC=(TH1F*)Gev->Clone();
            // Начальные значения параметров;
            fc->SetParameter(0,-5); fc->SetParameter(1,3);
            fc->SetParameter(2,-6); fc->SetParameter(3,2.35);
            fc->SetParameter(4,400); fc->SetParameter(5,2); 

             
   TCanvas *c=new TCanvas("c","Graph Draw Options",800,500);
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
        Gev->Fit(fc,"R");
        Gev->Draw();

         p=(fc->GetParameter(0));q=(fc->GetParameter(1));r=(fc->GetParameter(2));
         s=(fc->GetParameter(3));nmax=(fc->GetParameter(4));teta=(fc->GetParameter(5));
         
TLine *line1=new TLine(nmax,1e-6,nmax,0.5);
line1->Draw("SAME");
   
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
   //Gev->SetStats(0);      // No statistics on lower plot
  
    Gev->GetYaxis()->SetTitle("1/N dN_{ch}/dN");
    Gev->GetXaxis()->SetRangeUser(0,Nn);
    GevC->GetYaxis()->SetTitle("Fit/Data");
    GevC->GetXaxis()->SetRangeUser(0,Nn);
    GevC->GetXaxis()->SetTitle("N_{ch}");
// Y axis ratio plot settings

    Gev->GetYaxis()->SetNdivisions(505);
    Gev->GetYaxis()->SetTitleSize(20);
    Gev->GetYaxis()->SetTitleFont(43);
    Gev->GetYaxis()->SetTitleOffset(1.1);
    Gev->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Gev->GetYaxis()->SetLabelSize(15);
// X axis ratio plot settings
    
    Gev->GetXaxis()->SetTitleSize(20);
    Gev->GetXaxis()->SetTitleFont(43);
    Gev->GetXaxis()->SetTitleOffset(1.5);
    Gev->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Gev->GetXaxis()->SetLabelSize(15);


// Y axis ratio plot settings
    
    GevC->GetYaxis()->SetNdivisions(505);
    GevC->GetYaxis()->SetTitleSize(20);
    GevC->GetYaxis()->SetTitleFont(43);
    GevC->GetYaxis()->SetTitleOffset(1.1);
    GevC->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    GevC->GetYaxis()->SetLabelSize(15);
// X axis ratio plot settings

    GevC->GetXaxis()->SetTitleSize(20);
    GevC->GetXaxis()->SetTitleFont(43);
    GevC->GetXaxis()->SetTitleOffset(3);
    GevC->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    GevC->GetXaxis()->SetLabelSize(15);

    GevC->Draw();
TLine *line=new TLine(0,1,Nn,1);
line->Draw("SAME");


        };


double Norm(double step)
        {
            double Integ=0.0,n=0; int M=round(Nn/step);
            for (int k=0;k<M;k++)
            {
                if(k%2==0)
                    {
                        Integ+=2*Pn(n,p,q,r,s,nmax,teta);
                    }
                else
                    {
                        Integ+=4*Pn(n,p,q,r,s,nmax,teta);
                    }
                n+=step;       
            }
            return (Integ-(Pn(0,p,q,r,s,nmax,teta)+Pn(Nn,p,q,r,s,nmax,teta)))*step/3;
        };

 //Метод для нахождения интеграла P(n);
        double IntPn(double n0,double nn,double stepn)
        {
            double Integ=0.0,n=n0,norm=Norm(stepn); int M=round((nn-n0)/stepn);
            for (int k=0;k<M;k++)
            {
                if(k%2==0)
                    {
                        Integ+=2*Pn(n,p,q,r,s,nmax,teta);
                    }
                else
                    {
                        Integ+=4*Pn(n,p,q,r,s,nmax,teta);
                    }
                n+=stepn;       
            }
            return (Integ-(Pn(n0,p,q,r,s,nmax,teta)+Pn(nn,p,q,r,s,nmax,teta)))*stepn/(3*norm);
        };

//Метод для нахождения центральности;
        double c(double n0,double stepn)
        {
            double I=0.0;
            I=IntPn(n0,Nn+1,stepn);
        return I;
        };

//Метод для нахождения n по заданной центральности;
        double cFind(double c0,double n0,double stepn) {
            double n=n0,cen=0,step=1;
            while (cen<c0)
            {
                n=n-step; cen=c(n,stepn); 
            }
        return n;}

//Метод для нахождения распределения прецельного параметра для диапозона центральностей;

        double IntPb(double n0,double nn,double stepn,double b)
        {
            double I1=0,I2=0;
            double n=n0;
            int M=round((nn-n0)/stepn);
            for (int k=0;k<M;k++)
            {
                if(k%2==0)
                    {
                        I1+=2*Pnb(n,cb(b),p,q,r,s,nmax,teta);I2+=2*Pn(n,p,q,r,s,nmax,teta);
                    }
                else
                    {
                        I1+=4*Pnb(n,cb(b),p,q,r,s,nmax,teta);I2+=4*Pn(n,p,q,r,s,nmax,teta);
                    }
                n+=stepn;       
            }
            I1=(I1-(Pnb(n0,cb(b),p,q,r,s,nmax,teta)+Pnb(nn,cb(b),p,q,r,s,nmax,teta)));
            I2=(I2-(Pn(n0,p,q,r,s,nmax,teta)+Pn(nn,p,q,r,s,nmax,teta)));
       
        return 2*pi*b*I1/(sigma*I2);
        } ;



//Запись и вывод полученного распределения;
       vector <double> PlotIntb(double n0,double nn,double stepn,double bn,int N)
        {
            vector<double> Pbi;
            Pbi.clear();
            double bi=0,stepb=bn/N;
                        for (int n=0;n<=N;n++)
            {
                bi=stepb*n;
                Pbi.push_back(IntPb(n0,nn,stepn,bi));
            }
            return NormVect(Pbi,stepb); 
         };


TH1F* Histb(double n0,double nn,double bn,int N)
{ 
TH1F *hbH=new TH1F("hbH","impact",N,0,bn);
double k=N/bn;
int nentries = (int)t->GetEntries();


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



TH1F* Histnch(double b0,double bb)
{ 
TH1F *hnH=new TH1F("N_{ch}","UrQMD",250,0,500);

int nentries = (int)t->GetEntries();
for (int i=0; i<nentries; i++) 
{
    t->GetEntry(i);
    if (bimp >b0 && bimp <bb) 
    {  
        hnH->Fill(npart);
     } 

}
        hnH->Scale(0.5/hnH->GetEntries());
        return hnH;
}


void PlotPn(double b0,double bb,int Nc){

        TH1F *hn=new TH1F("","UrQMD",250,0,500); 
        hn=Histnch(b0,bb);

        vector<double> d1,Vn;
         for (int n=0;n<=500;n++)
            {
                Vn.push_back(n);
            }
        d1=Pnx(p,q,r,s,nmax,teta,cb(b0),cb(bb),Nc);
        d1=NormVect(d1,1);
               
            

TGraph* grb1=new TGraph(Vn.size(),Vn.data(),d1.data());
grb1->SetLineWidth(3);
grb1->SetTitle("Reconstracted");
grb1->GetYaxis()->SetTitle("P(N_{ch})");
grb1->GetXaxis()->SetTitle("N_{ch}");

grb1->GetYaxis()->SetTitleOffset(1.4);
grb1->Draw();
hn->Draw("SAME");
//gPad->BuildLegend();
}


void PlotMeann(int N,int Nc)
{
        vector<double> VPn,Vn,Cb,meanFit,meanData,ErorFit,ErorData;
        double stepb=bmax/N,mFit,mData,erFit,erData;
        for (int n=0;n<=500;n++)
            {
                Vn.push_back(n);
            }
    for(int i=0;i<=N;i++)
    {
            VPn.clear();
            Cb.push_back(cb(stepb*i));
            VPn=Pnx(p,q,r,s,nmax,teta,0,Cb[i],Nc);
            mFit=meanb(VPn,Vn,1);
            meanFit.push_back(mFit);
            //ErorFit.push_back(RMSb(VPn,Vn,1));
            mData=Histnch(0,stepb*i)->GetMean();
            meanData.push_back(mData);
            //erData=Histnch(stepb*i)->GetRMS();
            //ErorData.push_back(erData);
     }

//vector<double> ex(N+1,0);
//TGraphErrors *grb1 = new TGraphErrors(N+1,Cb.data(),meanFit.data(),ex.data(),ErorFit.data());
//TGraphErrors *grb2 = new TGraphErrors(N+1,Cb.data(),meanData.data(),ex.data(),ErorData.data());
TCanvas *c3=new TCanvas("c3","Graph Draw Options",10,0,800,800); c3->cd();
    TGraph* grb1=new TGraph(Cb.size(),Cb.data(),meanFit.data());
    grb1->GetYaxis()->SetTitle("<N_{ch}(c_{b})>");
    grb1->GetXaxis()->SetTitle("c_{b}");
    grb1->GetXaxis()->SetRangeUser(0,1);
    grb1->SetLineColor(1);
    grb1->SetMarkerStyle(21);
    grb1->SetMarkerSize(1);
    grb1->SetTitle("Reconstracted");
    grb1->Draw("AP");
    grb1->SetMarkerColor(1);
TGraph* grb2=new TGraph(Cb.size(),Cb.data(),meanData.data());
grb2->SetLineColor(2);
grb2->SetTitle("UrQMD");
grb2->SetMarkerStyle(20);
grb2->SetMarkerSize(1);
grb2->SetMarkerColor(2);
grb2->Draw("P SAME");
gPad->BuildLegend();
}

void Plotbn(double n0,double nn,double stepn,double bn,int N)
        {
         vector<double> b,Nc;
         vector<double> d1;
         d1.clear();
         double stepb=bn/N;
         for (int n=0;n<=N;n++)
            {
                b.push_back(stepb*n);
            }

        d1=PlotIntb(n0,nn,stepn,bn,N);
        d1=NormVect(d1,stepb);
        TGraph* grb1=new TGraph(N,b.data(),d1.data());grb1->SetTitle(" ");

        grb1->GetYaxis()->SetRangeUser(0,1);
        grb1->GetXaxis()->SetRangeUser(0,bn);
        
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1);
   pad2->SetBottomMargin(0.2);
        grb1->SetLineColor(1);
        grb1->GetXaxis()->SetTitle("b");
        grb1->GetYaxis()->SetTitle("P(b)");
  
    grb1->GetYaxis()->SetTitleSize(20);
    grb1->GetYaxis()->SetTitleFont(43);
    grb1->GetXaxis()->SetTitleOffset(1.1);
grb1->GetYaxis()->SetTitleOffset(1.1);
    grb1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    grb1->GetYaxis()->SetLabelSize(16);
// X axis ratio plot settings
    
    grb1->GetXaxis()->SetTitleSize(20);
    grb1->GetXaxis()->SetTitleFont(43);
 
    grb1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    grb1->GetXaxis()->SetLabelSize(16);
   
        grb1->Draw(); 

        TH1F *hb=new TH1F("hb","impact",N,0,bn); 
        hb=Histb(n0,nn,bn,N);

        hb->Draw("SAME");

}




      
void PlotAllb(double stepn,double bn,int N)
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
                Nc[n]=cFind(c0,Nc[n-1],stepn);
                c0+=0.1;
            }

        for (int n=0;n<10;n++)
            {   
                d[n]=PlotIntb(Nc[n+1],Nc[n],stepn,bn,N);
            }
    
    
double x[10],ex[10],x1[10],ex1[10],y[10],ey[10];       
TH1F *hb[10];
        for (int n=0;n<10;n++)
            {   
                hb[n]=Histb(Nc[n+1],Nc[n],bn,N);
                x[n]=meanb(d[n],b,stepb);
                ex[n]=RMSb(d[n],b,stepb);
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
gr->SetLineColor(2);
gr1->SetLineColor(1);
gr->SetTitle("Reconstracted");
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



TGraph* grb1=new TGraph(N,b.data(),d[0].data());grb1->SetName("grBimp_1");grb1->SetTitle("0-10%");
        
        TGraph* grb2=new TGraph(N,b.data(),d[1].data()); grb2->SetName("grBimp_2"); grb2->SetTitle("10-20%");
       
        TGraph* grb3=new TGraph(N,b.data(),d[2].data());grb3->SetName("grBimp_3"); grb3->SetTitle("20-30%");
        
       TGraph* grb4=new TGraph(N,b.data(),d[3].data());grb4->SetName("grBimp_4"); grb4->SetTitle("30-40%");
        TGraph* grb5=new TGraph(N,b.data(),d[4].data());grb5->SetName("grBimp_5"); grb5->SetTitle("40-50%");
       TGraph* grb6=new TGraph(N,b.data(),d[5].data());grb6->SetName("grBimp_6"); grb6->SetTitle("50-60%");
       TGraph* grb7=new TGraph(N,b.data(),d[6].data());grb7->SetName("grBimp_7"); grb7->SetTitle("60-70%");
        TGraph* grb8=new TGraph(N,b.data(),d[7].data());grb8->SetName("grBimp_8"); grb8->SetTitle("70-80%");
        TGraph* grb9=new TGraph(N,b.data(),d[8].data());grb9->SetName("grBimp_9"); grb9->SetTitle("80-90%");
     TGraph* grb10=new TGraph(N,b.data(),d[9].data());grb10->SetName("grBimp_10"); grb10->SetTitle("90-100%");
       
        grb1->GetYaxis()->SetRangeUser(0,1);
        grb1->GetXaxis()->SetRangeUser(0,bn);
  
        TCanvas *c5=new TCanvas("c5","Graph Draw Options",10,0,800,800);
      

        c5->cd();
        grb1->SetLineColor(1);
grb1->GetYaxis()->SetTitle("P(b/c)");
grb1->GetXaxis()->SetTitle("b");
        grb2->SetLineColor(2);
        grb3->SetLineColor(8);
        grb4->SetLineColor(4);
        grb5->SetLineColor(5);

        grb6->SetLineColor(1);grb6->SetLineStyle(9);
        
        grb7->SetLineColor(2);grb7->SetLineStyle(9);
        
        grb8->SetLineColor(8);grb8->SetLineStyle(9);
        grb9->SetLineColor(4);grb9->SetLineStyle(9);
        grb10->SetLineColor(5);grb10->SetLineStyle(9);


grb1->Draw(); 
        grb2->Draw("SAME"); 
        grb3->Draw("SAME"); 
        grb4->Draw("SAME"); 
       grb5->Draw("SAME"); 
        grb6->Draw("SAME"); 
        grb7->Draw("SAME"); 
        grb8->Draw("SAME"); 
        grb9->Draw("SAME"); 
        grb10->Draw("SAME"); 

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
