#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include <stdio.h>
#include "TROOT.h"
using namespace std;

void fit_fft_to_gaussian(double *pars, double *f, double *out, int N)
{
  //fit to gaussian, appropriate for non-rectangular windowing
  //find max (initial guess at Gaussian center)
  int imax = 0;
  double omax = -1;
  for (int i = 0; i < N; i++) {
    if (out[i] > omax) {
      omax = out[i];
      imax = i;
    }
  }
  //find 10% height (initial guess at Gaussian width)
  int i10 = N - 1;
  for (int i = imax; i < N; i++) {
    if (out[i] < omax / 10 && out[i + 1] < omax / 10) {
      i10 = i;
      break;
    }
  }
  int iwid = i10 - imax;
  //  cout << imax << " " << i10 << " " << iwid << endl;

  //to do the fit we will generate a TGraph and fit that
  //OK, let's generate a pointer to a subarray of OUT
  const int nP = 12 * iwid;
  double *subout = &out[imax - nP / 2];
  double *subx = &f[imax - nP / 2];

  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;

  TGraph *data = new TGraph(nP, subx, subout);
  data->SetName("gausFit");
  TF1 *fgaus = new TF1("fgaus", "[2]*TMath::Gaus(x,[0],[1])", subx[0], subx[nP - 1]);
  fgaus->SetParameters(f[imax], f[i10] - f[imax], out[imax] * 1e6);
  fgaus->SetNpx(1000);
  gStyle->SetOptFit(1111);
  int status = data->Fit(fgaus, "q", "", subx[0], subx[nP - 1]);
  //data->Draw("A*");

//OK, great, let's try to return the fit parameters
  pars[0] = fgaus->GetParameter(0);
  pars[1] = fgaus->GetParameter(1);
  pars[2] = fgaus->GetParameter(2);
  pars[3] = 0;
  pars[4] = fgaus->GetParError(0);
  pars[5] = fgaus->GetParError(1);
  pars[6] = fgaus->GetParError(2);
  pars[7] = status;

  data->Write();
}

void fit_fft_to_lorentian(double *pars, double *f, double *pow, int N, int j)
{
  cout << "****************************************************** " << endl;
  cout << "Fitting Power Spectrum P(f) to Lorentzian Function: " << endl;
  //fit to lorentian, appropriate for non-rectangular windowing
  //find max (initial guess at mean freq)
  int imax = 0;
  double omax = -1;
  for (int i = 0; i < N; i++) {
    if (pow[i] > omax) {
      omax = pow[i];
      imax = i;
    }
  }
  //find 10% height (initial guess at width)
  int i10 = N - 1;
  for (int i = imax; i < N; i++) {
    if (pow[i] < omax / 10 && pow[i + 1] < omax / 10) {
      i10 = i;
      break;
    }
  }
  int iwid = i10 - imax;
  //cout << imax << " " << i10 << " " << iwid << endl;

  //to do the fit we will generate a TGraph and fit that
  //OK, let's generate a pointer to a subarray of OUT
  const int nP = iwid * 24;
  double *subout = &pow[imax - nP / 2];
  double *subx = &f[imax - nP / 2];
  //create error bars, affects chi2 but not fit
  double erY[nP];
  for (int i = 0; i < nP; i++) {
    erY[i] = pow[imax] * .01;
  }

  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;

  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;

  TGraphErrors *data = new TGraphErrors(nP, subx, subout, 0, erY);
  data->SetName(Form("lorentzFit%i", j));

  TF1 *lor = new TF1("lorentz", "[2]*TMath::CauchyDist(x, [0], [1])", subx[0], subx[nP - 1]);
  lor->SetParameters(f[imax], f[i10] - f[imax], pow[imax]);
  lor->SetParNames("mean", "width", "amplitude");
  //TF1* lor = new TF1("lorentz","[2]*TMath::Min(1.0/(1.0-TMath::Power((x-[0])/[1],2)),0.0)",subx[0],subx[nP-1]);
  //lor->SetParameters(subx[3*iwid],subx[3*iwid]-subx[2*iwid],1);
  lor->SetLineColor(2);
  lor->SetNpx(1000);
  gStyle->SetOptFit(1111);
  int status = data->Fit(lor, "q", "", subx[0], subx[nP - 1]);
  //data->Draw("Ap");
  double integral = lor->Integral(subx[0], subx[nP - 1]);
  double ta_power = integral * f[1];
  cout << "Lorentzian function integral : " << integral << " fJ " << endl;
  cout << "Lorentzian time : " << 1 / f[1] << " s " << endl;
  cout << "Lorentzian function time average power: " << integral * f[1] << " fW " << endl;
  cout << "****************************************************** " << endl;
  //make sure integral and frequency have same units.
  //OK, great, let's try to return the fit parameters
  pars[0] = lor->GetParameter(0);
  pars[1] = lor->GetParameter(1);
  pars[2] = lor->GetParameter(2);
  pars[3] = ta_power;
  pars[4] = lor->GetParError(0);
  pars[5] = lor->GetParError(1);
  pars[6] = lor->GetParError(2);
  pars[7] = status;

  data->Write();

}

void fit_fft_to_sinc(double *pars, double *f, double *pow, double f0, int j)
{
  cout << "****************************************************** " << endl;
  cout << "Fitting Power Spectrum P(f) to Sinc2(f) Function: " << endl;
  //appropriate for rectangular windowing
  //find max (initial guess at mean freq)
  int imax = 0;
  double delf = f[1];
  imax = f0/delf; 
  cout << "Peak found at freq: " << f0 << endl;
  //to do the fit we will generate a TGraph and fit that
  //OK, let's generate a pointer to a subarray of OUT
  const int nP = 18;
  double *subout = &pow[imax - nP / 2];
  double *subx = &f[imax - nP / 2];
  //create error bars
  double erY[nP];
  for (int i = 0; i < nP; i++) {
    erY[i] = pow[imax] * .01;
  }

  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;
  //cout << "Nsamples " << N << endl;

  TGraphErrors *data = new TGraphErrors(nP, subx, subout, 0, erY);
  data->SetName(Form("sincFit%i", j));

  //for large sampling rate
  TF1 *sinc = new TF1("sinc", "[2]*TMath::Sin(2*TMath::Pi()*(x-[0])*[1]/2)^2/(2*TMath::Pi()*(x-[0])*[1]/2)^2", subx[0], subx[nP - 1]);
  sinc->SetParameters(f0, 1 / f[1], pow[imax]);
  //cout << "Initializing to mean: " << f[imax] << " duration: " << 1 / f[1] << " amp: " << pow[imax] << endl;
  sinc->SetParLimits(0, subx[0], subx[nP - 1]);
  sinc->SetParLimits(1, 0, 1000 / f[1]);
  sinc->SetParNames("mean", "duration", "amplitude");

  sinc->SetLineColor(2);
  sinc->SetNpx(1000);
  gStyle->SetOptFit(1111);
  int status = data->Fit(sinc, "q", "", subx[0], subx[nP - 1]);
  if (!status == 0) {
    status = data->Fit(sinc, "Mq", "", subx[0], subx[nP - 1]);
  }
  //data->Draw("Ap");
  double integral = sinc->Integral(subx[0], subx[nP - 1]);
  //make sure integral and frequency have same units.

  //OK, great, let's try to return the fit parameters
  pars[0] = data->GetFunction("sinc")->GetParameter(0);
  pars[1] = data->GetFunction("sinc")->GetParameter(1);
  pars[2] = data->GetFunction("sinc")->GetParameter(2);
  pars[3] = integral;
  pars[4] = sinc->GetParError(0);
  pars[5] = sinc->GetParError(1);
  pars[6] = sinc->GetParError(2);
  pars[7] = status;
  cout << "Sinc duration: " << pars[1] << " s " << endl;
  cout << "Sinc function TA power : " << integral << " fW " << endl;

  data->Write();
  //sinc->Write();
  cout << "****************************************************** " << endl;
}

void fit_fft_to_sinc_2nd(double *pars, double *f, double *pow, int N, int j)
{
  //fit 2nd harmonic to sinc^2, appropriate for rectangular windowing
  //find max (initial guess at mean freq)
  int imax = 0;
  double omax = -1;
  for (int i = 0; i < N; i++) {
    if (pow[i] > omax) {
      omax = pow[i];
      imax = i;
    }
  }
  //find 10% height (initial guess at width)
  int i10 = N - 1;
  for (int i = imax; i < N; i++) {
    if (pow[i] < omax / 10 && pow[i + 1] < omax / 10) {
      i10 = i;
      break;
    }
  }
  int iwid = i10 - imax;
  //cout << imax << " " << i10 << " " << iwid << endl;

  //to do the fit we will generate a TGraph and fit that
  //OK, let's generate a pointer to a subarray of OUT
  const int nP = iwid * 24;
  double *subout = &pow[2 * imax - nP / 2];
  double *subx = &f[2 * imax - nP / 2];
  //create error bars
  double erY[nP];
  for (int i = 0; i < nP; i++) {
    erY[i] = pow[2 * imax] * .01;
  }

  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;

  TGraphErrors *data = new TGraphErrors(nP, subx, subout, 0, erY);
  data->SetName(Form("sincFit_2nd_%i", j));

  TF1 *sinc = new TF1("sinc", "[2]*TMath::Sin(2*TMath::Pi()*(x-[0])*[1]/2)^2/(2*TMath::Pi()*(x-[0])*[1]/2)^2", subx[0], subx[nP - 1]);
  sinc->SetParameters(2 * f[2 * imax], 1 / f[1], pow[2 * imax]);
  sinc->SetParNames("mean", "duration", "amplitude");
  sinc->SetLineColor(2);
  sinc->SetNpx(1000);
  gStyle->SetOptFit(1111);
  int status = data->Fit(sinc, "q", "", subx[0], subx[nP - 1]);
  //data->Draw("Ap");
  double integral = sinc->Integral(subx[0], subx[nP - 1]);
  double ta_power = integral * f[1];
  cout << "Sinc function integral : " << integral << " fJ " << endl;
  cout << "Sinc duration : " << 1 / f[1] << " s " << endl;
  cout << "Sinc function time average power: " << integral * f[1] << " fW " << endl;
  //make sure integral and frequency have same units.

  //OK, great, let's try to return the fit parameters
  pars[0] = sinc->GetParameter(0);
  pars[1] = sinc->GetParameter(1);
  pars[2] = sinc->GetParameter(2);
  pars[3] = ta_power;
  pars[4] = sinc->GetParError(0);
  pars[5] = sinc->GetParError(1);
  pars[6] = sinc->GetParError(2);
  pars[7] = status;

  data->Write();
}

void fit_pow_to_cos(double *pars, double *t, double *pow_t, double f0, int N, int j)
{
  cout << "****************************************************** " << endl;
  cout << "Fitting Instantaneous Power P(t) to cos2(t)" << endl;
  //appropriate for single harmonic
  //find max in p(t) (initial guess at amp)
  double duration = t[0]*N;
  int imax_t = 0;
  double omax_t = -1;
  int nP = 50;
  for (int i = 0; i < nP; i++) {
    if (pow_t[i] > omax_t) {
      omax_t = pow_t[i];
      imax_t = i;
    }
  }
  //create error bars
  double erY[nP];
  for (int i = 0; i < nP; i++) {
    erY[i] = pow_t[imax_t] * .01;
  }
  //create TGraph
  TGraphErrors *data = new TGraphErrors(nP, t, pow_t, 0, erY);
  data->SetName(Form("cosFit%i", j));


  TF1 *cos = new TF1("cosine", "[2]*TMath::Cos(2*TMath::Pi()*x*[0]+[1])^2", t[0], t[nP - 1]);
  cos->SetParameters(f0, 0, pow_t[imax_t]);
  cos->SetParNames("frequency", "offset", "amplitude");
  cos->SetLineColor(2);
  cos->SetNpx(1000);
  gStyle->SetOptFit(1111);
  int status = data->Fit(cos, "q", "", t[0], t[nP - 1]);
  //data->Draw("Ap");
  double integral = cos->Integral(t[0], duration );
  double ta_power = integral / duration;
  cout << "cos function integral : " << integral << " fJ " << endl;
  cout << "cos function duration : " << duration << " s " << endl;
  cout << "cos function time average power: " << integral / t[N - 1] << " fW " << endl;
  //make sure integral and time have same units.

  pars[0] = cos->GetParameter(0);
  pars[1] = cos->GetParameter(1);
  pars[2] = cos->GetParameter(2);
  pars[3] = ta_power;
  pars[4] = cos->GetParError(0);
  pars[5] = cos->GetParError(1);
  pars[6] = cos->GetParError(2);
  pars[7] = status;

  data->Write();
  cout << "****************************************************** " << endl;
}

void add_noise(double *f, double *pow, int N, int j)
{
  double thermalNoise = 1.4e-6; //fJ 
  double sigLevel = 5.6 * 1.4e-6;  //fJ 
  //double thermalNoise = 1.4e-6*f[1];//fW 
  //double sigLevel= 5.6*1.4e-6*f[1];//fW 
  cout << "Power sig: " << sigLevel << endl;
  //find max (initial guess at mean)
  int imax = 0;
  double omax = -1;
  for (int i = 0; i < N; i++) {
    if (pow[i] > omax) {
      omax = pow[i];
      imax = i;
    }
  }
  //OK, let's generate a pointer to a subarray of OUT
  const int nP = 100;
  double *subout = &pow[imax - nP / 2];
  double *subx = &f[imax - nP / 2];
  //add noised
  for (int i = 0; i < nP; i++) {
    subout[i] += thermalNoise;
  }

  //  cout << iwid << " " << subout[0] << " " << subout[3*iwid] << endl;
  cout << "Nsamples " << N << endl;

  TGraphErrors *data = new TGraphErrors(nP, subx, subout, 0, 0);
  data->SetName(Form("noise%i", j));
  //noise level
  TF1 *noise = new TF1("noise", "[0]", subx[0], subx[nP - 1]);
  noise->SetTitle("Thermal Noise");
  noise->SetParameter(0, thermalNoise);
  noise->SetParNames("noise");
  noise->SetLineColor(4);
  noise->SetNpx(1000);

  //significance level
  TF1 *sig = new TF1("sig", "[0]", subx[0], subx[nP - 1]);
  sig->SetTitle("Required signal significance");
  sig->GetHistogram()->SetXTitle("Freq [Hz]");
  sig->GetHistogram()->SetYTitle("Time-averged power [fW]");
  sig->SetParameter(0, sigLevel);
  sig->SetParNames("sig");
  sig->SetLineColor(2);
  sig->SetNpx(1000);

  //data->Draw("AP");
  sig->Write();
  noise->Write();
  data->Write();
}
