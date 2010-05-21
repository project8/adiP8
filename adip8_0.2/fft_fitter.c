#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include <stdio.h>
#include "TROOT.h"
using namespace std;

void fit_fft_to_gaussian(double* pars,double* f,double* out,int N)
{
  //find max (initial guess at Gaussian center)
  int imax = 0; double omax=-1;
  for (int i=0;i<N;i++)
    {
      if (out[i] > omax) { omax = out[i]; imax = i; } 
    }
  //find 10% height (initial guess at Gaussian width)
  int i10 = N-1;
  for (int i=imax;i<N;i++)
    {
      if (out[i] < omax/10 && out[i+1] < omax/10) { i10 = i; break;}	
    }
  int iwid = i10 - imax;
  //  cout << imax << " " << i10 << " " << iwid << endl;

  //to do the fit we will generate a TGraph and fit that
  //OK, let's generate a pointer to a subarray of OUT
  const int nP = 12*iwid;
  double* subout = &out[imax-nP/2];
  double* subx = &f[imax-nP/2];
  
  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;
 
  TGraph* data = new TGraph(nP,subx,subout);
  data->SetName("gausFit");
  TF1* fgaus = new TF1("fgaus","[2]*TMath::Gaus(x,[0],[1])",subx[0],subx[nP-1]); 
  fgaus->SetParameters(f[imax], f[i10]-f[imax], out[imax]*1e6);
  fgaus->SetNpx(1000);
  gStyle->SetOptFit(1111);
  data->Fit(fgaus, "", "", subx[0], subx[nP-1]);
  data->Draw("A*");
  
//OK, great, let's try to return the fit parameters
  pars[0] = fgaus->GetParameter(0);
  pars[1] = fgaus->GetParameter(1);
  pars[2] = fgaus->GetParameter(2);
  pars[3] = 0;
  pars[4] = fgaus->GetParError(0);
  pars[5] = fgaus->GetParError(1);
  pars[6] = fgaus->GetParError(2);
  pars[7] = 0;
  pars[8] = fgaus->GetChisquare();

  data->Write();
}

void fit_fft_to_lorentzian(double* pars,double* f,double* out,int N)
{
  //find max (initial guess at Gaussian center)
  int imax = 0; double omax=-1;
  for (int i=0;i<N;i++)
    {
      if (out[i] > omax) { omax = out[i]; imax = i; } 
    }
  //find 10% height (initial guess at Gaussian width)
  int i10 = N-1;
  for (int i=imax;i<N;i++)
    {
      if (out[i] < omax/10 && out[i+1] < omax/10) { i10 = i; break;}	
    }
  int iwid = i10 - imax;
    //cout << imax << " " << i10 << " " << iwid << endl;

  //to do the fit we will generate a TGraph and fit that
  //OK, let's generate a pointer to a subarray of OUT
  const int nP = iwid*12;
  double* subout = &out[imax-nP/2];
  double* subx = &f[imax-nP/2];
  //create error bars, affects chi2 but not fit
  double erY[nP];
  for (int i=0;i<nP;i++)
    {
      erY[i] = out[imax]*.01;
    }
  
  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;
  
  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;
 
  TGraphErrors* data = new TGraphErrors(nP,subx,subout, 0, erY);
  data->SetName("lorentzFit");

  TF1* lor = new TF1("lorentz","[2]*TMath::CauchyDist(x, [0], [1])",subx[0],subx[nP-1]);
  lor->SetParameters(f[imax], f[i10]-f[imax], out[imax]*1e6);
  lor->SetParNames("mean", "width", "amplitutde");
  //TF1* lor = new TF1("lorentz","[2]*TMath::Min(1.0/(1.0-TMath::Power((x-[0])/[1],2)),0.0)",subx[0],subx[nP-1]);
  //lor->SetParameters(subx[3*iwid],subx[3*iwid]-subx[2*iwid],1);
  lor->SetLineColor(2);
  lor->SetNpx(1000);
  gStyle->SetOptFit(1111);
  data->Fit(lor, "", "", subx[0], subx[nP-1]);
  data->Draw("Ap");
  double integral = lor->Integral(subx[0], subx[nP-1]);
  cout << "Lorentzian function integral : " << integral << endl;
  cout << "Lorentzian time : " << 1/f[1] << endl;
  cout << "Lorentzian function time average: " << integral*f[1] << endl;
  //make sure integral and frequency have same units.
  //OK, great, let's try to return the fit parameters
  pars[0] = lor->GetParameter(0);
  pars[1] = lor->GetParameter(1);
  pars[2] = lor->GetParameter(2);
  pars[3] = 0;
  pars[4] = lor->GetParError(0);
  pars[5] = lor->GetParError(1);
  pars[6] = lor->GetParError(2);
  pars[7] = 0;
  pars[8] = lor->GetChisquare();

  lor->Write();
  data->Write();
}
void fit_fft_to_line_broadening(double* pars,double* f,double* out,int N)
{
  //find max (initial guess at Gaussian center)
  int imax = 0; double omax=-1;
  for (int i=0;i<N;i++)
    {
      if (out[i] > omax) { omax = out[i]; imax = i; } 
    }
  //find 10% height (initial guess at Gaussian width)
  int i10 = N-1;
  for (int i=imax;i<N;i++)
    {
      if (out[i] < omax/10 && out[i+1] < omax/10) { i10 = i; break;}	
    }
  int iwid = i10 - imax;
    //cout << imax << " " << i10 << " " << iwid << endl;

  //to do the fit we will generate a TGraph and fit that
  //OK, let's generate a pointer to a subarray of OUT
  const int nP = iwid*24; 
  double* subout = &out[imax-nP/2];
  double* subx = &f[imax-nP/2];
  //create error bars
  double erY[nP];
  for (int i=0;i<nP;i++)
    {
      erY[i] = out[imax]*.01;
    }
  
  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;
 
  TGraphErrors* data = new TGraphErrors(nP,subx,subout, 0, erY);
  data->SetName("broadFit");

  TF1* linebroad = new TF1("linebroadening","[2]*TMath::Sin(2*TMath::Pi()*(x-[0])*[1]/2)^2/(2*TMath::Pi()*(x-[0])*[1]/2)^2",subx[0],subx[nP-1]);
  linebroad->SetParameters(f[imax],1/f[1],out[imax]);
  linebroad->SetParNames("mean", "duration", "amplitutde");
  linebroad->SetLineColor(2);
  linebroad->SetNpx(1000);
  gStyle->SetOptFit(1111);
  data->Fit(linebroad, "", "", subx[0], subx[nP-1]);
  data->Draw("Ap");
  double integral = linebroad->Integral(subx[0], subx[nP-1]);
  cout << "Line Broadening function integral : " << integral << endl;
  cout << "Line Broadening time : " << 1/f[1] << endl;
  cout << "Line Broadening function time average: " << integral*f[1] << endl;
  //make sure integral and frequency have same units.
  
  //OK, great, let's try to return the fit parameters
  pars[0] = linebroad->GetParameter(0);
  pars[1] = linebroad->GetParameter(1);
  pars[2] = linebroad->GetParameter(2);
  pars[3] = 0;
  pars[4] = linebroad->GetParError(0);
  pars[5] = linebroad->GetParError(1);
  pars[6] = linebroad->GetParError(2);
  pars[7] = 0;
  pars[8] = linebroad->GetChisquare();

  linebroad->Write();
  data->Write();
}

