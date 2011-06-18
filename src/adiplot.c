#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include "fftw3.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TMath.h"
#include "radiation.h"
#include "paramanage.h"

using namespace std;
struct t_parameter parameter;

int main(int argc, char *argv[])
{
  TApplication theApp("App", &argc, argv);
  TString cardname = theApp.Argv(1);

  //use Adipark parameter-reading routine to read control card
  char *parameter_filename = new char[255];
  strcpy(parameter_filename, cardname.Data());
  init_parameters();
  int parameter_file_found = 0;
  if (load_init_data(parameter_filename) == 0) {

    char *parameter_filename_sub = new char[255];
    cout << "Try parameter file w/o suffix: ";
    strncat(parameter_filename_sub, parameter_filename, strlen(parameter_filename) - 4);
    cout << parameter_filename_sub << endl;
    if (load_init_data(parameter_filename_sub) == 0) {
      cout << "Try parameter file w/o point: ";
      strcpy(parameter_filename_sub, "");
      strncat(parameter_filename_sub, parameter_filename, strlen(parameter_filename) - 1);
      cout << parameter_filename_sub << endl;
      if (load_init_data(parameter_filename_sub) == 0) {
        cout << endl << "ERROR: missing some settings in INI-file!";
        cout << endl << endl << flush;
      } else {
        parameter_file_found = 1;
        cardname.Resize(strlen(parameter_filename) - 1);
      }
    } else {
      parameter_file_found = 1;
      cardname.Resize(strlen(parameter_filename) - 4);
    }
  } else {
    parameter_file_found = 1;
  }

  init_data();                  //set default tl_data
  switch (parameter.rad_calc_mode) {
    case 2:                      //parallel wires
      init_tl_data(kFALSE);
      break;
    case 3:                      //square wg
      init_sq_wg_data(OMEGA0 / Clight / M2CM);
      //Warning!  Wave imp. for TE modes freq dep, not implemented properly
      break;
    case 4:                      //circular wg
      init_circ_wg_data(OMEGA0 / Clight / M2CM);
      //Warning!  Wave imp. for TE modes freq dep, not implemented properly
      break;
    case 5:                      //parallel plates/strips
      init_pp_data();
      break;
    case 6:                      //coaxial cable
      init_coax_data();
      break;
    case 7:                      //offset parallel wires
      init_tl_data(kTRUE);
      break;
  }
  double kT = (parameter.antenna_temp) * K_BOL;    //fJ, T=35K, amplifier Teff = 25 K
  if (kT == 0 ) kT = K_BOL;
  double impedance = parameter.impedance; //normalized load impedance
  double refCo = (impedance - 1) / (impedance + 1);
  if (std::isinf(impedance)) {
    refCo = 1;
  }

  TString rootname = cardname + TString(".root");
  TFile *tfin = new TFile(rootname, "read");

  TTree *runcard = (TTree *) tfin->Get("runcard");
  Long64_t nEvents = runcard->GetEntries();
  tfin->Close();

  TFile *tfout = new TFile(cardname + "_extended.root", "update");
  TNtuple *eloss = new TNtuple("eloss", "Parameters for correcting eloss", "i:z0:zcut:dfdt:foundStart:nfound");
  float ntarray[4];

  //loop over number of electrons in file 
  for (Int_t event = 1; event <= nEvents; event++) {

    TTree *wf_1 = (TTree *) tfout->Get(Form("wf_%d", event));
    Long64_t nEntries = wf_1->GetEntries();
    cout << cardname << " open with " << nEntries << endl;

    //initialize values
    wf_1->GetEntry(1);
    double delt = wf_1->FindLeaf("t")->GetValue();//in secs
    double tmax = parameter.max_tof_in_sec;
    double tSize = 5.0e-5;      //in seconds, < 50 us if free space energy loss
    if (refCo == 1) {
      tSize = 1.3e-5;           //in seconds, 4xfree space energy loss
    }
    if ( tmax < tSize) tSize = tmax; 
    Int_t Ntot = tSize / delt;  //max size of fft
    double delf = 1 / tSize;    //in Hz
    const Int_t nX = nEntries / Ntot;     //number of ffts
    cout << "Number of transforms: " << nX;
    cout << " " << tSize << " s long " << endl;
    Int_t nY = Ntot / 2 + 1;//number of freq bins
    double *wf;
    wf = (double *) fftw_malloc(sizeof(double) * Ntot);
    fftw_complex *out_fd;
    out_fd = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Ntot);
    fftw_plan p_flong;
    double *powerf = (double *) fftw_malloc(sizeof(double) * (Ntot / 2 + 1));

    p_flong = fftw_plan_dft_r2c_1d(Ntot, wf, out_fd, FFTW_ESTIMATE);
    //spectrogram 
    TH2F *h2 = new TH2F(Form("fvst_%d", event), "Spectrogram;Time (ms);Frequency (MHz); Power (fW)", nX, -0.5 * tSize * 1e3, (nX - 0.5) * tSize * 1e3, 
                                                                                                    nY, -0.5 * delf * 1e-6, (nY - 0.5) * delf * 1e-6);
    //noise in single PS
    TH1F *hn = new TH1F(Form("hfWn_%d", event), "Noise; Power (fW); Counts", nY, 0, 0.1);
    TH1F *hPSa = new TH1F(Form("hPS_avg_%d", event), "Averaged Power Spectra; Frequency (MHz); Power (fW)", nY, -0.5 * delf * 1e-6, (nY - 0.5) * delf * 1e-6);
    TH1F **hPS = new TH1F *[nX];
    TH1F *hna = new TH1F(Form("hfWn_avg_%d", event), "Averaged Noise; Power (fW); Counts", nY, 0, 1.1);
    cout << " max time : " << tSize * nX << endl;
    int bin;
    double noiseP = kT * delf;
    double zCut = 4.6 + log(nY * nX / 4);   //requirement power*noiseP to accept as signal
    double z0 = 4.6 + log(nY * nX);    //requirement power*noiseP to accept as signal
    Bool_t signalFT = kFALSE;
    cout << "Required signal: " << zCut * noiseP << endl;
    int NsigFT = 0;
    double slope = 3.6e8;       //delf/delt for energy loss
    double maxP = 0;
    double t0_i = 0;
    //loop over number of fourier transforms
    for (int i = 0; i < nX; i++) {
      //cout << "time " << (i+0.5)*tSize << endl;
      for (int k = 0; k < Ntot; k++) {
        wf_1->GetEntry(k + i * Ntot);
        wf[k] = wf_1->FindLeaf("nEf")->GetValue();
      }
      //from time to freq. 
      fftw_execute(p_flong);
      //calculate power and energy from long sample fft
      for (int j = 0; j <= Ntot / 2; j++) {
        bin = h2->GetBin(i + 1, j + 1);
        //mean square amplitude, time-averaged power at bin f is fW_f
        //should have units of fW 
        h2->SetBinContent(bin, 1e+15 / tl_data.Zw / Ntot / Ntot * (pow(out_fd[j][0], 2) + pow(out_fd[j][1], 2)));
        powerf[j] = 1e+15 / tl_data.Zw / Ntot / Ntot * (pow(out_fd[j][0], 2) + pow(out_fd[j][1], 2));
      }
      //now add in negative frequencies 
      for (int j = Ntot / 2 + 1; j < Ntot; j++) {
        bin = h2->GetBin(i + 1, Ntot - j + 1);
        //mean square amplitude, time-averaged power at bin f is fW_f
        //should have units of fW 
        h2->AddBinContent(bin, 1e+15 / tl_data.Zw / Ntot / Ntot * (pow(out_fd[Ntot - j][0], 2) + pow(out_fd[Ntot - j][1], 2)));
        powerf[Ntot - j] += 1e+15 / tl_data.Zw / Ntot / Ntot * (pow(out_fd[Ntot - j][0], 2) + pow(out_fd[Ntot - j][1], 2));
        if ((powerf[Ntot - j] > zCut * noiseP) || (powerf[Ntot - j] + powerf[Ntot - (j - 1)] > zCut * noiseP)) {
          signalFT = kTRUE;
          if (NsigFT == 0) {
            t0_i = i;
          }
        }
        if (powerf[Ntot - j] > maxP) {
          maxP = powerf[Ntot - j];
        }
      }
      if (signalFT) {
        //add this FT into signal plot, correct for energy loss 
        int shift = int ((i - t0_i) * slope / delf / delf + 0.5);    //round to nearest integer 
        hPS[NsigFT] = new TH1F(Form("hPS_%d_%d", NsigFT, event), Form("Power Spectra %d w/ signal; Frequency (MHz); Power (fW)", i), nY, -0.5 * delf * 1e-6, (nY - 0.5) * delf * 1e-6);
        for (int j = 0; j <= Ntot / 2; j++) {
          hPSa->AddBinContent(j, powerf[j + shift]);
          hPS[NsigFT]->SetBinContent(j, powerf[j]);
        }
        cout << "Averaging in PS: " << i << " at time " << i * tSize << " with shift" << shift << endl;
        hPS[NsigFT]->Write();
        NsigFT += 1;
      }
      if (!signalFT) {
        //check that noise is correct distribution
        for (int j = 0; j <= Ntot / 2; j++) {
          hn->Fill(powerf[j]);
        }
      }
      signalFT = kFALSE;
    }
    cout << "Max power: " << maxP << endl;
    cout << "Averaged PS: " << NsigFT << endl;
    hPSa->SetTitle(Form("%d Averaged Power Spectra", NsigFT));
    hPSa->Scale(1.0 / NsigFT);
    //check averaged plot has correct noise
    double tempP;
    for (int j = 0; j <= Ntot / 2; j++) {
      tempP = hPSa->GetBinContent(j);       //check that noise is correct distribution
      hna->Fill(tempP);
    }
    ntarray[0] = event;
    ntarray[1] = z0;
    ntarray[2] = zCut;
    ntarray[3] = slope;
    ntarray[4] = t0_i * tSize;
    ntarray[5] = NsigFT;
    eloss->Fill(ntarray);
    h2->Write();
    hn->Write();
    hna->Write();
    hPSa->Write();
  }
  eloss->Write();
  //tfout->Write();
  tfout->Close();
}
